# USP5 Baseline Analysis Summary

## Dataset
- Total input rows: 26
- Valid SMILES rows: 26
- Unique molecules after canonical deduplication: 19
- Measured rows (`ic50 > 0`): 16
- Assigned-label rows (`ic50 in {0, -1}`): 10
- Duplicate molecule groups: 6
- Duplicate groups with conflicting pIC50 values: 3

## Chemical diversity
- Unique Bemis-Murcko scaffolds: 18
- Mean nearest-neighbor Tanimoto: 0.430
- Max nearest-neighbor Tanimoto: 0.975

## Modeling
- Best LOOCV model: knn_k1_distance
- Best feature set: descriptors
- MAE: 0.752
- RMSE: 1.303
- R2: -0.286

## Interpretation
- This is a very small dataset, so LOOCV estimates can still be unstable.
- The best setup was selected by comparing many feature and model configurations on the same dataset, so the top score is optimistic.
- Assigned labels are mixed with measured potency values, so predictive performance should be treated as provisional.
- Duplicate molecules with inconsistent activities suggest assay or curation variability and should be reviewed before trusting fine-grained rankings.
- The next rational improvement is a measured-only sensitivity analysis and more experimentally measured compounds, not a more complex model.