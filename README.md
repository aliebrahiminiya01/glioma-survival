# Glioma Survival Analysis (TCGA LGG + GBM)

Predicting overall survival of diffuse glioma patients from clinical and molecular features.

Data: cBioPortal `lgggbm_tcga_pub` (Ceccarelli et al., Cell 2016). After cleaning, 1040 patients, 524 deaths.

## Setup
```bash
python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python src/prepare_data.py      # downloads + builds data/processed/glioma.csv
.venv/bin/python -m ipykernel install --user --name glioma-survival --display-name "Glioma Survival (.venv)"
```
Then open the notebooks with the **Glioma Survival (.venv)** kernel.

## Steps
- [x] 1. Data: download + clean → `src/prepare_data.py`
- [x] 2. EDA, Kaplan-Meier, log-rank, Cox PH → `notebooks/02_eda_km_cox.ipynb`
- [x] 3. Proper evaluation: train/test split + CV, C-index, time-dependent AUC, Brier score → `notebooks/03_evaluation.ipynb`, `src/survival_utils.py`
- [x] 4. ML models: Random Survival Forest, Gradient Boosting Survival (nested CV) → `notebooks/04_tree_models.ipynb`
- [ ] 5. DeepSurv: postponed until imaging features (UNET / radiomics) are available
- [x] 6. Interpretation (hazard ratios, permutation importance) → `notebooks/06_interpretation.ipynb`, final report → **[REPORT.md](REPORT.md)**
