# USP5 Inhibitor Final Workflow

This project now centers on one canonical USP5 virtual screening workflow based on the saved final exploratory regression model, a broad 10-method enumeration library, and a multistage lead-selection funnel.

The canonical workflow is:

1. Build or reuse the saved final `ExtraTreesRegressor` potency model in [`outputs/final_model/`](/Users/shivanshsahni/Documents/New%20project/outputs/final_model).
2. Enumerate a broad virtual library with the 10-method chemistry pipeline in [`scripts/run_enumeration_10_methods.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_enumeration_10_methods.py).
3. Screen that broad library with the restored final-model lead funnel in [`scripts/run_lead_selection.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_lead_selection.py).
4. Report primary leads and orthogonal backup leads from the final outputs in [`outputs/final_leads.csv`](/Users/shivanshsahni/Documents/New%20project/outputs/final_leads.csv) and [`outputs/backup_leads.csv`](/Users/shivanshsahni/Documents/New%20project/outputs/backup_leads.csv).

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

## What the canonical workflow does

1. Validates required columns and parses the dataset.
2. Validates and canonicalizes SMILES with RDKit.
3. Preserves row-level provenance while also generating a deduplicated modeling table by canonical SMILES.
4. Generates RDKit physicochemical and graph descriptors.
5. Uses the saved final `ExtraTreesRegressor` as the project-level potency model.
6. Expands all positive training chemotypes into a broad `3264`-compound virtual library using 10 enumeration techniques.
7. Screens that library with:
   - `PAINS + BRENK`
   - Lipinski-style property gates
   - `TPSA` and molecular surface area limits
   - `ADMET-AI`
   - multistructure USP5 3D template-docking / pharmacophore plausibility using `6DXT`, `7MS5`, `7MS6`, and `7MS7`
8. Produces a primary lead set plus orthogonal backup chemotypes.

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

## Canonical project layout

- [`data/raw/First.csv`](/Users/shivanshsahni/Documents/New%20project/data/raw/First.csv): copied source dataset
- [`outputs/raw_dataset_summary.md`](/Users/shivanshsahni/Documents/New%20project/outputs/raw_dataset_summary.md): initial summary extracted from the uploaded CSV
- [`outputs/final_pipeline_summary.md`](/Users/shivanshsahni/Documents/New%20project/outputs/final_pipeline_summary.md): canonical final workflow summary
- [`scripts/run_workflow.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_workflow.py): baseline data/model workflow
- [`scripts/run_enumeration_10_methods.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_enumeration_10_methods.py): canonical broad enumeration workflow
- [`scripts/run_lead_selection.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_lead_selection.py): canonical final lead-selection workflow
- [`scripts/run_lead_selection_final_model_restored.py`](/Users/shivanshsahni/Documents/New%20project/scripts/run_lead_selection_final_model_restored.py): restored final-model implementation
- [`scripts/lead_selection_multistructure_common.py`](/Users/shivanshsahni/Documents/New%20project/scripts/lead_selection_multistructure_common.py): shared multistructure and ADMET screening helper logic
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

## Canonical outputs

- `annotated_rows.csv`: row-level cleaned data with label flags
- `modeling_dataset.csv`: deduplicated per-molecule table used for training
- `outputs/final_model/`: final saved regression model and report
- `enumeration_library_10_methods.csv`: canonical broad virtual library
- `enumeration_method_counts.csv`: library-size breakdown by method
- `lead_selection_counts.csv`: canonical lead-selection attrition table
- `lead_selection_strict_pool.csv`: strict multistage survivor pool
- `lead_selection_relaxed_pool.csv`: relaxed 3D-screened pool
- `final_leads.csv`: canonical primary leads
- `backup_leads.csv`: canonical orthogonal backup leads
- `analysis_summary.md`: baseline interpretation of dataset quality and model limitations

## Final selected model

The canonical project-level potency model is the saved raw-row `base_graph` `ExtraTreesRegressor` in [`outputs/final_model/final_model.joblib`](/Users/shivanshsahni/Documents/New%20project/outputs/final_model/final_model.joblib), selected because it achieved an in-sample `R^2` of about `0.893` while still using interpretable RDKit-derived chemistry features.

Rebuild the final report and artifacts with:

```bash
MPLCONFIGDIR=.matplotlib PYTHONPATH=src .venv/bin/python scripts/build_final_model_report.py
```

Final model deliverables are written under [`outputs/final_model/`](/Users/shivanshsahni/Documents/New%20project/outputs/final_model).

## Final screening pipeline

The canonical final lead-selection workflow is the restored-final-model pipeline:

- potency by the saved final `ExtraTreesRegressor`
- `PAINS + BRENK` cleanup
- `Lipinski`, `TPSA`, molecular surface area, flexibility, and charge filters
- `ADMET-AI`
- multistructure USP5 ZnF-UBD 3D plausibility against `6DXT`, `7MS5`, `7MS6`, and `7MS7`

The main final screening summary is in [`outputs/lead_selection_summary.md`](/Users/shivanshsahni/Documents/New%20project/outputs/lead_selection_summary.md).

## Archiving policy

Older exploratory scripts and generated artifacts are moved under archive folders so the top-level `outputs/` and `scripts/` directories point to the canonical final workflow rather than multiple competing variants.

## Recommended next improvement

After the baseline is run, the most rational next step is not a more complex model. It is to improve label quality and data coverage:

- separate measured from assigned labels in every analysis view
- inspect duplicate molecules with conflicting activities
- expand the measured potency set before considering more flexible models

If the baseline shows some signal, a reasonable follow-up is a measured-only sensitivity analysis and a simple similarity-aware applicability-domain check, not deep learning.
