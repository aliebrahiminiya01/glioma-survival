"""
Shared helpers for Steps 3-5: loading, a fixed train/test split, preprocessing,
and survival metrics. Every model (Cox, RSF, DeepSurv) is evaluated with the
same functions, so the comparison table is fair.

A "model" here is any sklearn-style estimator with
    fit(X, y)                     y = sksurv structured array (event, time)
    predict(X)                    -> risk score, higher = dies earlier
    predict_survival_function(X)  -> list of step functions S(t)   (sksurv)
or a LifelinesCox wrapper below.
"""
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sklearn.base import BaseEstimator, clone
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.metrics import (concordance_index_censored, concordance_index_ipcw,
                            cumulative_dynamic_auc, integrated_brier_score)
from sksurv.util import Surv

NUMERIC = ["age", "kps", "grade"]
CATEGORICAL = ["sex", "subtype", "mgmt", "chr7_10"]
# Reference level per categorical (dropped in one-hot) -> HR is relative to these.
# Same references as Step 2: the higher-risk group, so HR < 1 means protective.
REFERENCE = {"sex": "Female", "subtype": "IDHwt", "mgmt": "Unmethylated", "chr7_10": "no"}

# Median follow-up is ~12 months and only ~66 patients are still at risk at 60 months,
# so metrics beyond 36 months would rest on a handful of test patients.
EVAL_TIMES = np.array([6, 12, 24, 36])
SEED = 42


def load_data(path="data/processed/glioma.csv"):
    """-> X (DataFrame of features), y (sksurv structured array)."""
    df = pd.read_csv(path, index_col="patient_id")
    X = df[NUMERIC + CATEGORICAL].copy()
    X["grade"] = X["grade"].astype(float)
    y = Surv.from_arrays(event=df["event"].astype(bool), time=df["time"])
    return X, y


def split(X, y, test_size=0.2):
    """Fixed hold-out test set (same SEED in every step). Stratified on the event
    so train and test have the same death/censoring ratio."""
    return train_test_split(X, y, test_size=test_size, stratify=y["event"], random_state=SEED)


def make_preprocessor(numeric=NUMERIC, categorical=CATEGORICAL):
    """Fitted inside each CV fold, so test-fold statistics never leak into training.

    numeric:     median imputation + "was missing" flag (KPS is missing for ~1/3) + scaling
    categorical: missing becomes its own category ("missing" can carry information,
                 e.g. molecular tests not done) + one-hot, REFERENCE level dropped
                 so Cox is not collinear.
    """
    num = make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler())
    cat = make_pipeline(SimpleImputer(strategy="constant", fill_value="missing"),
                        OneHotEncoder(drop=[REFERENCE[c] for c in categorical],
                                      handle_unknown="ignore", sparse_output=False))
    pre = ColumnTransformer([("num", num, list(numeric)), ("cat", cat, list(categorical))],
                            verbose_feature_names_out=False)
    return pre.set_output(transform="pandas")


class LifelinesCox(BaseEstimator):
    """lifelines CoxPHFitter in sklearn form, so we can use stratification
    (scikit-survival's Cox does not support strata).

    strata_prefix: every column starting with it becomes a stratum variable,
    e.g. "subtype_" -> each subtype gets its own baseline hazard. Fixes the
    proportional-hazards violation found in Step 2.
    """

    def __init__(self, penalizer=0.01, strata_prefix=None):
        self.penalizer = penalizer
        self.strata_prefix = strata_prefix

    def _frame(self, X, y=None):
        df = pd.DataFrame(X).copy()
        if y is not None:
            df["time"], df["event"] = y["time"], y["event"].astype(int)
        return df

    def fit(self, X, y):
        df = self._frame(X, y)
        strata = [c for c in df.columns if self.strata_prefix and c.startswith(self.strata_prefix)]
        self.cph_ = CoxPHFitter(penalizer=self.penalizer)
        self.cph_.fit(df, "time", "event", strata=strata or None)
        return self

    def predict(self, X):
        return self.cph_.predict_log_partial_hazard(self._frame(X)).to_numpy()

    def predict_survival_at(self, X, times):
        # lifelines returns (times x patients); we want (patients x times).
        # With strata it also regroups patients by stratum -> restore the input order.
        df = self._frame(X)
        surv = self.cph_.predict_survival_function(df, times=times)
        return surv[df.index].T.to_numpy()


def survival_at(model, X, times):
    """S(t | x) for every patient at every time -> array (n_patients, n_times)."""
    model = getattr(model, "best_estimator_", model)   # GridSearchCV -> the refit best model
    if hasattr(model, "predict_survival_at"):
        return model.predict_survival_at(X, times)
    if hasattr(model, "steps") and hasattr(model.steps[-1][1], "predict_survival_at"):
        Xt = model[:-1].transform(X)
        return model.steps[-1][1].predict_survival_at(Xt, times)
    return np.vstack([fn(times) for fn in model.predict_survival_function(X)])


def evaluate(model, y_train, X_test, y_test, times=EVAL_TIMES):
    """All metrics for an already fitted model on one test set.

    harrell_c : classic C-index (biased when censoring is heavy)
    uno_c     : IPCW C-index, truncated at the last eval time (more honest here)
    auc_<t>   : time-dependent AUC - can the model tell who dies before t vs after t?
    ibs       : integrated Brier score - calibration + discrimination of S(t),
                lower is better (0.25 ~ a useless model that always predicts 0.5)
    """
    risk = model.predict(X_test)
    res = {
        "harrell_c": concordance_index_censored(y_test["event"], y_test["time"], risk)[0],
        "uno_c": concordance_index_ipcw(y_train, y_test, risk, tau=times.max())[0],
    }
    auc, mean_auc = cumulative_dynamic_auc(y_train, y_test, risk, times)
    res.update({f"auc_{t}m": a for t, a in zip(times, auc)})
    res["mean_auc"] = mean_auc
    res["ibs"] = integrated_brier_score(y_train, y_test, survival_at(model, X_test, times), times)
    return res


def cross_validate(model, X, y, n_splits=5, n_repeats=3, times=EVAL_TIMES):
    """Repeated stratified K-fold on the TRAINING set only -> one row per fold."""
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)
    rows = []
    for fold, (tr, va) in enumerate(cv.split(X, y["event"])):
        m = clone(model).fit(X.iloc[tr], y[tr])
        row = {"fold": fold, **evaluate(m, y[tr], X.iloc[va], y[va], times)}
        if hasattr(m, "best_params_"):   # nested CV: record what the inner search picked
            row["params"] = str(m.best_params_)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(cv_results: dict) -> pd.DataFrame:
    """{name: per-fold DataFrame} -> table of mean ± std per model."""
    out = {}
    for name, r in cv_results.items():
        r = r.drop(columns=["fold", "params"], errors="ignore")
        out[name] = r.mean().round(3).astype(str) + " ± " + r.std().round(3).astype(str)
    return pd.DataFrame(out).T


class KaplanMeierBaseline(BaseEstimator):
    """No-feature baseline: everyone gets the same KM curve and the same risk.
    C-index/AUC = 0.5 by construction; its IBS is the number a real model must beat."""

    def fit(self, X, y):
        from lifelines import KaplanMeierFitter
        self.km_ = KaplanMeierFitter().fit(y["time"], y["event"])
        return self

    def predict(self, X):
        return np.zeros(len(X))

    def predict_survival_at(self, X, times):
        s = self.km_.survival_function_at_times(times).to_numpy()
        return np.tile(s, (len(X), 1))
