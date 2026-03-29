# Refined Lead Selection Summary

This refined funnel is stricter than the first pass and is designed to behave more like a real early lead-selection workflow.

## Stage counts

- `start_similarity_window_library`: 2255
- `stage1_alert_free_pains_brenk_nih`: 505
- `stage2_measured_parent_and_parent_pIC50_ge_4p5`: 344
- `stage3_lead_like_property_window`: 152
- `stage4_developability_admet_proxies`: 124
- `stage5_similarity_sweet_spot`: 69

## Stage logic

- Stage 1: remove PAINS, BRENK, and NIH alerts.
- Stage 2: require the nearest original parent to be a measured positive with pIC50 >= 4.5.
- Stage 3: tighter lead-like property window.
- Stage 4: developability/ADMET proxy filters using flexibility, charge, QED, and SA score.
- Stage 5: keep compounds in a focused similarity sweet spot (0.65 to 0.88).

## Final leads

- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)ccnc2c1=O` | score 0.7513 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.773 | MW 430.8 | cLogP 1.79 | TPSA 123.4 | QED 0.558 | SA 2.45 | transform `heteroatom_CH_to_N`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)N2CCOCC2)nc2c(Cl)cccc2c1=O` | score 0.7470 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.761 | MW 499.0 | cLogP 2.17 | TPSA 102.8 | QED 0.510 | SA 2.42 | transform `rxn_acid_to_amide_morpholine`
- `COc1cccnc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | score 0.7468 | nearest parent `CHEMBL5410606` (pIC50 5.10) | sim 0.776 | MW 430.8 | cLogP 1.79 | TPSA 123.4 | QED 0.558 | SA 2.46 | transform `heteroatom_CH_to_N`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)nn1` | score 0.7292 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.737 | MW 438.9 | cLogP 1.51 | TPSA 129.6 | QED 0.699 | SA 2.45 | transform `heteroatom_CH_to_N`
- `O=C(NCC(=O)N1CCCC1)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.7035 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.793 | MW 491.0 | cLogP 2.66 | TPSA 99.7 | QED 0.671 | SA 2.33 | transform `rxn_acid_to_amide_pyrrolidine`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cn3)CC2)cn1` | score 0.6848 | nearest parent `CHEMBL5278336` (pIC50 5.15) | sim 0.789 | MW 438.9 | cLogP 1.51 | TPSA 129.6 | QED 0.699 | SA 2.45 | transform `heteroatom_CH_to_N`
- `Cc1cc(C(=O)CN2CCC(O)CC2)c(C)n1-c1ccc(C#N)nc1` | score 0.6742 | nearest parent `CHEMBL5276142` (pIC50 4.67) | sim 0.732 | MW 338.4 | cLogP 2.00 | TPSA 82.2 | QED 0.862 | SA 2.54 | transform `heteroatom_CH_to_N`
- `Cc1cc(C(=O)CN2CCC(O)CC2)c(C)n1-c1ncc(C#N)cn1` | score 0.6403 | nearest parent `CHEMBL5276142` (pIC50 4.67) | sim 0.709 | MW 339.4 | cLogP 1.40 | TPSA 95.0 | QED 0.846 | SA 2.67 | transform `heteroatom_CH_to_N`
- `Cc1nc(C(=O)CN2CCC(O)CC2)c(C)n1-c1ccc(C#N)cc1` | score 0.5556 | nearest parent `CHEMBL5276142` (pIC50 4.67) | sim 0.661 | MW 338.4 | cLogP 2.00 | TPSA 82.2 | QED 0.862 | SA 2.57 | transform `heteroatom_CH_to_N`
- `Cc1cc(C(=O)CN2CCCCC2)c(C)n1-c1cnc(F)nc1` | score 0.5329 | nearest parent `CHEMBL1493046` (pIC50 4.54) | sim 0.708 | MW 316.4 | cLogP 2.69 | TPSA 51.0 | QED 0.643 | SA 2.54 | transform `ring_pyrrolidine_to_piperidine`