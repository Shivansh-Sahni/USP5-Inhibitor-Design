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
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(F)cc3)CC2)cn1` | score 0.6600 | predicted pIC50 5.011 (pred IC50 9.75 uM) | AMES 0.098 | hERG 0.123 | HIA 0.900 | Bioavailability 0.884 | shape 0.619 | pharmacophore 0.972 | clashes 0 | transform `mmp_aryl_Cl_to_F`
- `NC(=O)CNC(=O)c1ncc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.6027 | predicted pIC50 4.932 (pred IC50 11.69 uM) | AMES 0.103 | hERG 0.219 | HIA 0.995 | Bioavailability 0.930 | shape 0.687 | pharmacophore 0.806 | clashes 0 | transform `bioisostere_acid_to_amide`
- `NC(=O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cn3)CC2)cn1` | score 0.5733 | predicted pIC50 4.933 (pred IC50 11.68 uM) | AMES 0.063 | hERG 0.243 | HIA 0.997 | Bioavailability 0.918 | shape 0.648 | pharmacophore 0.778 | clashes 0 | transform `bioisostere_acid_to_amide`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cn3)CC2)cn1` | score 0.5461 | predicted pIC50 4.895 (pred IC50 12.74 uM) | AMES 0.036 | hERG 0.131 | HIA 0.928 | Bioavailability 0.902 | shape 0.593 | pharmacophore 0.944 | clashes 0 | transform `heteroatom_CH_to_N`
- `O=C(O)CNC(=O)c1cnc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.4884 | predicted pIC50 4.895 (pred IC50 12.74 uM) | AMES 0.110 | hERG 0.147 | HIA 0.938 | Bioavailability 0.905 | shape 0.638 | pharmacophore 0.694 | clashes 0 | transform `heteroatom_CH_to_N`
- `CNC(=O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(F)cn3)CC2)cn1` | score 0.3929 | predicted pIC50 4.752 (pred IC50 17.71 uM) | AMES 0.057 | hERG 0.179 | HIA 0.999 | Bioavailability 0.930 | shape 0.592 | pharmacophore 0.611 | clashes 0 | transform `rxn_acid_to_amide_methylamine`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Br)cc3)CC2)cn1` | score 0.3253 | predicted pIC50 4.628 (pred IC50 23.55 uM) | AMES 0.086 | hERG 0.192 | HIA 0.896 | Bioavailability 0.880 | shape 0.610 | pharmacophore 1.000 | clashes 0 | transform `mmp_aryl_Cl_to_Br`