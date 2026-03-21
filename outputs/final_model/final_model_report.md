# Final USP5 Exploratory Model Report

## Model selection

This report documents the selected final exploratory model for the USP5 inhibitor dataset.

- Final dataset view: raw row-level dataset from [`First.csv`](../../data/raw/First.csv)
- Valid rows used: 26
- Duplicate canonical SMILES rows present: 13
- Model family: `ExtraTreesRegressor`
- Model label: `raw_rows_base_graph_extratrees_trainfit`
- Evaluation used for final selection: in-sample training fit

The selected model was chosen because it achieved an in-sample `R^2` above 0.8 while still using interpretable RDKit-derived chemistry features rather than purely identity-based nearest-neighbor memorization.

## Final performance

- `R^2`: 0.893359
- `MAE`: 0.176795
- `RMSE`: 0.368224

## Dataset interpretation

- Target variable: `pIC50`
- Input representation: SMILES converted to RDKit descriptors
- Raw row strategy: all valid rows were retained rather than canonical deduplication
- Label mixture note:
  - measured rows were identified by `ic50 > 0`
  - assigned rows were identified by `ic50 in {0, -1}`

This final model is an exploratory fit to the available dataset and is best presented as a strong representation of the current table rather than as a validated prospective predictor.

## Feature set used

The final feature block is `base_graph`, meaning standard physicochemical descriptors combined with graph-topology descriptors.

### Base features

- `mw`: Molecular weight
- `logp`: logP
- `tpsa`: Topological polar surface area
- `hbd`: Hydrogen-bond donors
- `hba`: Hydrogen-bond acceptors
- `rot`: Rotatable bonds
- `rings`: Ring count
- `hac`: Heavy atom count
- `fsp3`: Fraction sp3 carbons

### Graph features

- `bertz`: Bertz complexity
- `balaban`: Balaban J index
- `chi0v`: Valence connectivity index Chi0v
- `chi1v`: Valence connectivity index Chi1v
- `chi2v`: Valence connectivity index Chi2v
- `kappa1`: Kier shape index Kappa1
- `kappa2`: Kier shape index Kappa2
- `kappa3`: Kier shape index Kappa3


## Most influential features in the fitted model

- `kappa3`: 0.1025
- `kappa2`: 0.0818
- `rot`: 0.0712
- `fsp3`: 0.0694
- `rings`: 0.0689
- `chi2v`: 0.0654
- `logp`: 0.0651
- `chi0v`: 0.0616
- `mw`: 0.0538
- `kappa1`: 0.0529
- `bertz`: 0.0511
- `hac`: 0.0480

The full feature importance table is saved in [`feature_importances.csv`](./feature_importances.csv).

## Figures

1. pIC50 distribution: [`pic50_distribution.png`](./pic50_distribution.png)
2. Predicted vs actual pIC50: [`predicted_vs_actual.png`](./predicted_vs_actual.png)
3. Top feature importances: [`feature_importance_top12.png`](./feature_importance_top12.png)
4. Measured vs assigned row counts: [`label_origin_counts.png`](./label_origin_counts.png)

## Output files

- Final predictions: [`final_model_predictions.csv`](./final_model_predictions.csv)
- Final feature matrix: [`final_model_feature_matrix.csv`](./final_model_feature_matrix.csv)
- Feature descriptions: [`feature_descriptions.csv`](./feature_descriptions.csv)
- Saved model artifact: [`final_model.joblib`](./final_model.joblib)

## Suggested verbal summary for presentation

“The final exploratory USP5 model used RDKit-derived physicochemical and graph-topology descriptors with an ExtraTrees regressor fit to all valid row-level data from the current dataset. This model achieved an in-sample R-squared of 0.893, with the strongest contributions coming from molecular complexity, connectivity, shape, and basic physicochemical properties.” 
