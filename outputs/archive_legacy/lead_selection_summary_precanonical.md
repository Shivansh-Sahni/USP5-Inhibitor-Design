# Lead Selection Summary

Starting library: similarity-window enumeration (`0.55-0.95` vs original positives).

## Stage counts

- `start_similarity_window_library`: 2255
- `stage1_alert_free_pains_brenk_nih`: 505
- `stage2_lead_like_core_properties`: 374
- `stage3_veber_surface_area_flexibility`: 367
- `stage4_admet_proxy_qed_charge_similarity_focus`: 274

## Notes

- Stage 1 removes PAINS, BRENK, and NIH alerts.
- Stage 2 applies core lead-like property limits.
- Stage 3 applies Veber/surface-area/flexibility limits.
- Stage 4 is an ADMET proxy stage using QED, charge sanity, and a tighter similarity focus.
- Dedicated hERG, Ames, or Kd models were not available locally, so these were not run as true predictive models.

## Final leads

- `O=C(NCC(=O)N1CCOCC1)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.8296 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.754 | MW 507.0 | cLogP 1.89 | TPSA 108.9 | QED 0.641 | SA 2.41 | transform `rxn_acid_to_amide_morpholine`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)nn1` | score 0.8244 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.737 | MW 438.9 | cLogP 1.51 | TPSA 129.6 | QED 0.699 | SA 2.45 | transform `heteroatom_CH_to_N`
- `O=C(NCC(=O)N1CCCCC1)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.8215 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.780 | MW 505.0 | cLogP 3.05 | TPSA 99.7 | QED 0.651 | SA 2.35 | transform `rxn_acid_to_amide_piperidine`
- `Cc1cc(C(=O)CN2CCC(O)CC2)c(C)n1-c1ccc(C#N)nc1` | score 0.8028 | nearest parent `CHEMBL5276142` (pIC50 4.67) | sim 0.732 | MW 338.4 | cLogP 2.00 | TPSA 82.2 | QED 0.862 | SA 2.54 | transform `heteroatom_CH_to_N`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)N2CCOCC2)nc2c(Cl)cccc2c1=O` | score 0.8014 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.761 | MW 499.0 | cLogP 2.17 | TPSA 102.8 | QED 0.510 | SA 2.42 | transform `rxn_acid_to_amide_morpholine`
- `COc1nnccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | score 0.7974 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.750 | MW 431.8 | cLogP 1.18 | TPSA 136.3 | QED 0.542 | SA 2.64 | transform `heteroatom_CH_to_N`
- `CCNC(=O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.7970 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.793 | MW 465.0 | cLogP 2.17 | TPSA 108.5 | QED 0.652 | SA 2.28 | transform `rxn_acid_to_amide_ethylamine`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)nncc2c1=O` | score 0.7946 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.750 | MW 431.8 | cLogP 1.18 | TPSA 136.3 | QED 0.542 | SA 2.57 | transform `heteroatom_CH_to_N`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cn3)CC2)cn1` | score 0.7916 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.789 | MW 438.9 | cLogP 1.51 | TPSA 129.6 | QED 0.699 | SA 2.45 | transform `heteroatom_CH_to_N`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)ccnc2c1=O` | score 0.7897 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.773 | MW 430.8 | cLogP 1.79 | TPSA 123.4 | QED 0.558 | SA 2.45 | transform `heteroatom_CH_to_N`