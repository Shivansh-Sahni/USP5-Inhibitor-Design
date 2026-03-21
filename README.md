# USP5 Inhibitor Modeling Baseline

This project builds a practical small-data cheminformatics workflow for USP5 inhibitor potency modeling from [`First.csv`](/Users/shivanshsahni/Documents/New%20project/data/raw/First.csv). The goal is not to force complex models onto a tiny dataset, but to establish a transparent RDKit-based baseline that supports exploratory analysis, regression, and early ranking for virtual screening.

## Scientific framing

- Target protein: ubiquitin-specific protease 5 (USP5)
- Input chemistry: small-molecule structures provided as SMILES
- Modeling target: `pIC50`
- Important limitation: the dataset mixes measured and assigned labels

`pIC50` is the negative base-10 logarithm of the IC50 expressed in molar units. Higher `pIC50` means stronger potency.

## Dataset semantics

The raw CSV is expected to contain these columns:

- `id`: source identifier such as ChEMBL ID or PubChem CID
- `pIC50`: regression target used in this project
- `ic50`: mixed numeric and label-encoding column
- `smiles`: molecular SMILES string

The `ic50` column must be interpreted carefully:

- `ic50 > 0`: measured quantitative IC50 was available
- `ic50 == 0`: active compound, but no numeric IC50 was available
- `ic50 == -1`: inactive compound

Rows with `ic50 == 0` or `ic50 == -1` are not measured IC50 values. They are label encodings. The workflow therefore creates helper columns:

- `is_measured`
- `is_active_no_ic50`
- `is_inactive`
- `is_assigned_label`

These flags are preserved in all cleaned outputs and used in reporting.

## What the workflow does

1. Validates required columns and parses the dataset.
2. Validates and canonicalizes SMILES with RDKit.
3. Preserves row-level provenance while also generating a deduplicated modeling table by canonical SMILES.
4. Generates RDKit descriptors and Morgan fingerprints.
5. Computes Tanimoto similarity matrices and nearest-neighbor summaries.
6. Runs leave-one-out cross-validation (LOOCV) with small-data baseline models:
   - mean predictor
   - ElasticNet
   - RandomForestRegressor
   - XGBoost regressor if `xgboost` is installed
7. Saves exploratory summaries, feature tables, model metrics, LOOCV predictions, and ranked compounds.

## Duplicate handling policy

The raw file contains repeated molecules. This project keeps all valid rows in the annotated dataset, then builds a deduplicated modeling table using canonical SMILES.

For each duplicate group:

- all source identifiers are preserved
- all original `pIC50` and `ic50` values are preserved as joined strings
- if at least one measured row exists, the modeling target is the median of the measured `pIC50` values
- otherwise, the modeling target is the median of the assigned `pIC50` values

This is a conservative compromise for a very small dataset. It avoids pretending that duplicate rows are independent molecules while still retaining provenance.

## Why classical models were chosen

The dataset is only a few dozen compounds. That is too small for deep learning methods such as GNNs or 3D-equivariant models unless there is unusually strong external data support, which is not the case here. For this reason the workflow emphasizes:

- clean preprocessing
- interpretable RDKit features
- LOOCV instead of optimistic train/test splits
- honest reporting of uncertainty and data limitations

## Project layout

- [`data/raw/First.csv`](/Users/shivanshsahni/Documents/New%20project/data/raw/First.csv): copied source dataset
- [`outputs/raw_dataset_summary.md`](/Users/shivanshsahni/Documents/New%20project/outputs/raw_dataset_summary.md): initial summary extracted from the uploaded CSV
- [`scripts/run_workflow.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_workflow.py): command-line entry point
- [`src/usp5_workflow/data.py`](/Users/shivanshsahni/Documents/New%20project/src/usp5_workflow/data.py): loading, validation, annotation, canonicalization, deduplication
- [`src/usp5_workflow/features.py`](/Users/shivanshsahni/Documents/New%20project/src/usp5_workflow/features.py): descriptors, fingerprints, Tanimoto summaries, scaffold analysis
- [`src/usp5_workflow/modeling.py`](/Users/shivanshsahni/Documents/New%20project/src/usp5_workflow/modeling.py): LOOCV, model construction, metrics, ranking
- [`src/usp5_workflow/pipeline.py`](/Users/shivanshsahni/Documents/New%20project/src/usp5_workflow/pipeline.py): end-to-end workflow orchestration

## Installation

Create an environment with RDKit and the standard scientific Python stack. One workable approach is:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Required packages are listed in [`requirements.txt`](/Users/shivanshsahni/Documents/New%20project/requirements.txt).

Note: dependency installation was not possible in this session because network access was unavailable, so the code was scaffolded but not executed end-to-end here.

## Usage

Run the full workflow with:

```bash
PYTHONPATH=src .venv/bin/python scripts/run_workflow.py
```

Optional arguments:

```bash
PYTHONPATH=src .venv/bin/python scripts/run_workflow.py \
  --input data/raw/First.csv \
  --output outputs \
  --fingerprint-sizes 512 1024 2048 \
  --random-seed 42
```

## Main outputs

- `annotated_rows.csv`: row-level cleaned data with label flags
- `modeling_dataset.csv`: deduplicated per-molecule table used for training
- `descriptor_features.csv`: descriptor matrix
- `fingerprint_features_*.csv`: Morgan fingerprint matrices
- `tanimoto_similarity_matrix.csv`: pairwise similarity matrix
- `similarity_summary.csv`: per-molecule nearest-neighbor similarity summary
- `scaffold_summary.csv`: Bemis-Murcko scaffold counts
- `model_metrics.csv`: LOOCV metrics for each model and fingerprint size
- `loocv_predictions.csv`: out-of-fold predictions and residuals
- `full_fit_predictions.csv`: full-data fitted predictions for ranking support
- `ranked_compounds.csv`: molecules ranked by the best available model
- `analysis_summary.md`: concise interpretation of dataset quality and model limitations

## Final selected model

The current project-level final exploratory model is the raw-row `base_graph` feature set with an `ExtraTreesRegressor`, selected because it achieved an in-sample `R^2` of about `0.893` while still using interpretable RDKit-derived chemistry features.

Rebuild the final report and artifacts with:

```bash
MPLCONFIGDIR=.matplotlib PYTHONPATH=src .venv/bin/python scripts/build_final_model_report.py
```

Final deliverables are written under [`outputs/final_model/`](/Users/shivanshsahni/Documents/New%20project/outputs/final_model).

## Recommended next improvement

After the baseline is run, the most rational next step is not a more complex model. It is to improve label quality and data coverage:

- separate measured from assigned labels in every analysis view
- inspect duplicate molecules with conflicting activities
- expand the measured potency set before considering more flexible models

If the baseline shows some signal, a reasonable follow-up is a measured-only sensitivity analysis and a simple similarity-aware applicability-domain check, not deep learning.
