# Model-First Lead Selection Summary

This funnel starts with the trained ExtraTrees potency model, then moves through PAINS/BRENK cleanup, lead-like property filters, ADMET-AI, and a structure-aware 3D binding plausibility screen against the USP5 `7MS7` co-crystal ligand pose.

## Stage counts

- `start_focused_library`: 2255
- `stage1_predicted_pIC50_ge_4.60`: 1050
- `stage2_pains_brenk_free`: 137
- `stage3_lipinski_psa_surface_area_flexibility`: 92
- `stage4_admet_ai_core_gates`: 59
- `stage5_template_docking_and_3d_pharmacophore`: 15

## Filter logic

- Stage 1: keep compounds with predicted `pIC50 >= 4.60` from the final ExtraTrees regression model.
- Stage 2: remove compounds with `PAINS` or `BRENK` alerts. `NIH` alerts were not used.
- Stage 3: require `Lipinski violations <= 1`, `TPSA 45-140`, `Labute ASA 140-235`, `rotatable bonds <= 10`, and `formal charge within -1 to +1`.
- Stage 4: run `admet-ai` and keep compounds with low `AMES`, low `hERG`, low `ClinTox`, plus acceptable `HIA` and `oral bioavailability`.
- Stage 5: align conformers into the bound USP5 `7MS7` ligand pose, require good 3D shape overlap, feature overlap, and low receptor clash count.

## Final leads

- `NC(=O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(F)cc3)CC2)cn1` | score 0.7123 | predicted pIC50 5.005 (pred IC50 9.88 uM) | AMES 0.172 | hERG 0.227 | HIA 0.994 | Bioavailability 0.901 | shape 0.780 | pharmacophore 0.778 | clashes 0 | transform `bioisostere_acid_to_amide`
- `NC(=O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.6965 | predicted pIC50 5.088 (pred IC50 8.17 uM) | AMES 0.135 | hERG 0.372 | HIA 0.995 | Bioavailability 0.886 | shape 0.757 | pharmacophore 0.722 | clashes 0 | transform `bioisostere_acid_to_amide`
- `O=C(O)CNC(=O)c1ncc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.6779 | predicted pIC50 4.895 (pred IC50 12.73 uM) | AMES 0.054 | hERG 0.106 | HIA 0.923 | Bioavailability 0.919 | shape 0.772 | pharmacophore 0.917 | clashes 0 | transform `heteroatom_CH_to_N`