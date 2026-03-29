# Final Canonical USP5 Pipeline

This project now uses one canonical end-to-end computational workflow.

## Final model

- Potency model: saved `ExtraTreesRegressor`
- Model artifact: [final_model.joblib](/Users/shivanshsahni/Documents/New%20project/outputs/final_model/final_model.joblib)
- Model report: [final_model_report.md](/Users/shivanshsahni/Documents/New%20project/outputs/final_model/final_model_report.md)
- Final reported in-sample `R^2`: `0.893359`

## Final enumeration workflow

- Script: [run_enumeration_10_methods.py](/Users/shivanshsahni/Documents/New%20project/scripts/run_enumeration_10_methods.py)
- Library file: [enumeration_library_10_methods.csv](/Users/shivanshsahni/Documents/New%20project/outputs/enumeration_library_10_methods.csv)
- Total unique enumerated compounds: `3264`
- Method breakdown: [enumeration_method_counts.md](/Users/shivanshsahni/Documents/New%20project/outputs/enumeration_method_counts.md)

## Final screening workflow

- Script: [run_lead_selection_final_model_restored.py](/Users/shivanshsahni/Documents/New%20project/scripts/run_lead_selection_final_model_restored.py)
- Summary: [lead_selection_summary.md](/Users/shivanshsahni/Documents/New%20project/outputs/lead_selection_summary.md)
- Primary leads: [final_leads.csv](/Users/shivanshsahni/Documents/New%20project/outputs/final_leads.csv)
- Backup leads: [backup_leads.csv](/Users/shivanshsahni/Documents/New%20project/outputs/backup_leads.csv)

## Final screening stages

- Start broad enumeration plus original positives: `3274`
- Saved final model predicted `pIC50 >= 4.60`: `1272`
- `PAINS + BRENK` free: `772`
- Lipinski / PSA / surface area / flexibility: `59`
- `ADMET-AI`: `24`
- Strict multistructure USP5 3D survivors: `4`
- Relaxed orthogonal backup pool: `5`

## Final interpretation

- The primary lead family is the `CHEMBL5278336` acid-sulfonamide series.
- The orthogonal backup program is centered on `CHEMBL5410606`-derived bicyclic carbonyl chemistry plus a weaker heteroaryl-acid backup branch.
- The folder has been cleaned so the top-level outputs point to this workflow, while older exploratory artifacts are preserved under archive folders.
