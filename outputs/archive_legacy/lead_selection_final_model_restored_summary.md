# Restored Final-Model Lead Selection

This run restores the original saved ExtraTrees regression model as the only potency model. The downstream workflow is still the stronger version: broad library input, PAINS/BRENK cleanup, lead-like property filters, ADMET-AI, and multistructure USP5 3D binding plausibility.

## Stage counts

- `start_broad_enumeration_plus_original_positives`: 3274
- `stage1_original_final_model_predicted_pIC50_ge_4.60`: 1272
- `stage2_pains_brenk_free`: 772
- `stage3_lipinski_psa_surface_area_flexibility`: 59
- `stage4_admet_ai_multigate`: 24
- `stage5_multistructure_template_docking_strict`: 4
- `stage5b_multistructure_template_docking_relaxed`: 5

## Primary leads

- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.6365 | predicted pIC50 5.197 (pred IC50 6.36 uM) | AMES 0.080 | hERG 0.213 | best binding 0.782 | scaffold `O=S(=O)(c1cccnc1)N1CCC(c2ccccc2)CC1` | parent `CHEMBL5278336`
- `O=C(O)CNC(=O)c1ncc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.4127 | predicted pIC50 4.895 (pred IC50 12.73 uM) | AMES 0.054 | hERG 0.106 | best binding 0.830 | scaffold `O=S(=O)(c1cncnc1)N1CCC(c2ccccc2)CC1` | parent `CHEMBL5278336`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(F)cc3)CC2)cn1` | score 0.4108 | predicted pIC50 5.011 (pred IC50 9.75 uM) | AMES 0.098 | hERG 0.123 | best binding 0.780 | scaffold `O=S(=O)(c1cccnc1)N1CCC(c2ccccc2)CC1` | parent `CHEMBL5278336`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cn3)CC2)cn1` | score 0.1993 | predicted pIC50 4.895 (pred IC50 12.74 uM) | AMES 0.036 | hERG 0.131 | best binding 0.762 | scaffold `O=S(=O)(c1cccnc1)N1CCC(c2ccccn2)CC1` | parent `CHEMBL5278336`

## Orthogonal backups

- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | score 0.5514 | predicted pIC50 5.100 | AMES 0.128 | hERG 0.225 | best binding 0.350 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` | parent `CHEMBL5410606`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(F)cccc2c1=O` | score 0.5056 | predicted pIC50 5.019 | AMES 0.128 | hERG 0.161 | best binding 0.397 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` | parent `CHEMBL5410606`
- `Cc1cc(-c2cc(C)n(CCC(=O)O)c2C)c(C)n1CCC(=O)O` | score 0.4640 | predicted pIC50 4.874 | AMES 0.093 | hERG 0.113 | best binding 0.556 | scaffold `c1cc(-c2cc[nH]c2)c[nH]1` | parent `[9*]n1c(C)cc([16*])c1C`
- `Cc1cc(-n2c(C)cc(CCC(=O)O)c2C)c(C)n1CCC(=O)O` | score 0.4548 | predicted pIC50 4.876 | AMES 0.068 | hERG 0.068 | best binding 0.500 | scaffold `c1ccn(-c2cc[nH]c2)c1` | parent `[9*]n1c(C)cc([16*])c1C`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cncc2c1=O` | score 0.4344 | predicted pIC50 4.877 | AMES 0.064 | hERG 0.179 | best binding 0.415 | scaffold `O=C(Cn1cnc2ccncc2c1=O)NCc1ccccc1` | parent `CHEMBL5410606`
- `COc1ncccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | score 0.4239 | predicted pIC50 4.876 | AMES 0.074 | hERG 0.073 | best binding 0.346 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1cccnc1` | parent `CHEMBL5410606`

## Interpretation

The original final model was preserved intact. The main scientific upgrades happen after potency scoring: better developability triage, ADMET-AI, and multistructure 3D evidence. This means the final list reflects the trusted model while still avoiding overreliance on a single chemotype in the final program recommendation.