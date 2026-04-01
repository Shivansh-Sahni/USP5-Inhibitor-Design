# Restored Final-Model Lead Selection

This run preserves the original saved ExtraTrees regression model but changes the final prioritization goal. The primary leads are now intentionally filtered away from the known ZnF-UBD-like chemistry and re-ranked for novelty relative to the existing molecules, while still preserving potency, ADMET, and lead-like property constraints.

## Stage counts

- `start_broad_enumeration_plus_original_positives`: 3274
- `stage1_original_final_model_predicted_pIC50_ge_4.60`: 1272
- `stage2_pains_brenk_free`: 772
- `stage3_lipinski_psa_surface_area_flexibility`: 59
- `stage4_admet_ai_multigate`: 24
- `stage5_non_znf_and_not_original`: 18
- `stage5b_non_znf_relaxed_pool`: 19

## Primary leads

- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(F)cccc2c1=O` | orthogonal score 0.7428 | predicted pIC50 5.019 (pred IC50 9.57 uM) | AMES 0.128 | hERG 0.161 | max existing similarity 0.812 | max ZnF-reference similarity 0.202 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` | parent `CHEMBL5410606`
- `Cc1cc(-n2c(C)cc(CCC(=O)O)c2C)c(C)n1CCC(=O)O` | orthogonal score 0.6583 | predicted pIC50 4.876 (pred IC50 13.31 uM) | AMES 0.068 | hERG 0.068 | max existing similarity 0.230 | max ZnF-reference similarity 0.186 | scaffold `c1ccn(-c2cc[nH]c2)c1` | parent `[9*]n1c(C)cc([16*])c1C`
- `COc1ccncc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | orthogonal score 0.6403 | predicted pIC50 4.876 (pred IC50 13.30 uM) | AMES 0.030 | hERG 0.178 | max existing similarity 0.818 | max ZnF-reference similarity 0.191 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1cccnc1` | parent `CHEMBL5410606`
- `Cc1cc(-c2cc(C)n(CCC(=O)O)c2C)c(C)n1CCC(=O)O` | orthogonal score 0.6396 | predicted pIC50 4.874 (pred IC50 13.36 uM) | AMES 0.093 | hERG 0.113 | max existing similarity 0.226 | max ZnF-reference similarity 0.217 | scaffold `c1cc(-c2cc[nH]c2)c[nH]1` | parent `[9*]n1c(C)cc([16*])c1C`
- `CC(=O)C1=C(O)C(=O)N(CCC(=O)O)C1c1ccc(Cl)cc1` | orthogonal score 0.5966 | predicted pIC50 4.908 (pred IC50 12.35 uM) | AMES 0.032 | hERG 0.030 | max existing similarity 0.527 | max ZnF-reference similarity 0.257 | scaffold `O=C1C=CC(c2ccccc2)N1` | parent `CHEMBL4129140`
- `CC(=O)C1=C(O)C(=O)N(c2ccc(O)cc2)C1c1ccc(Cl)cc1` | orthogonal score 0.5781 | predicted pIC50 4.878 (pred IC50 13.25 uM) | AMES 0.093 | hERG 0.224 | max existing similarity 0.407 | max ZnF-reference similarity 0.225 | scaffold `O=C1C=CC(c2ccccc2)N1c1ccccc1` | parent `CHEMBL4129140`

## Orthogonal backups

- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | orthogonal score 0.7123 | predicted pIC50 5.100 | AMES 0.128 | hERG 0.225 | max existing similarity 1.000 | max ZnF-reference similarity 0.187 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` | parent `CHEMBL5410606`
- `COc1ncccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | orthogonal score 0.5719 | predicted pIC50 4.876 | AMES 0.074 | hERG 0.073 | max existing similarity 0.788 | max ZnF-reference similarity 0.196 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1cccnc1` | parent `CHEMBL5410606`
- `CC(=O)C1=C(O)C(=O)N(c2ccccc2)C1c1ccc(Cl)cc1` | orthogonal score 0.5152 | predicted pIC50 4.864 | AMES 0.044 | hERG 0.176 | max existing similarity 0.407 | max ZnF-reference similarity 0.243 | scaffold `O=C1C=CC(c2ccccc2)N1c1ccccc1` | parent `CHEMBL4129140`
- `CC(=O)C1=C(O)C(=O)N(c2ccc(F)cc2)C1c1ccc(Cl)cc1` | orthogonal score 0.5054 | predicted pIC50 4.835 | AMES 0.058 | hERG 0.197 | max existing similarity 0.393 | max ZnF-reference similarity 0.219 | scaffold `O=C1C=CC(c2ccccc2)N1c1ccccc1` | parent `CHEMBL4129140`

## Interpretation

The original final model was preserved intact. The key change is post-model: known ZnF-like chemistry is explicitly deprioritized using a reference-similarity filter (`< 0.30` to ZnF templates), original known molecules are excluded from primary leads, and novelty versus the existing dataset is rewarded. This shifts the primary list toward orthogonal chemotypes with more distinct scaffolds and less dependence on the previously favored CHEMBL5278336-like series.