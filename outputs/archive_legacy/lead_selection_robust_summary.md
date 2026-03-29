# Robust Multi-Evidence USP5 Lead Selection

This report intentionally shifts away from over-reliance on any single signal. The workflow uses a broader starting library, a consensus potency score with uncertainty, PAINS/BRENK cleanup, lead-like properties, ADMET-AI, and multi-structure 3D binding plausibility against multiple USP5 ZnF-UBD co-crystals.

## Stage counts

- `start_broad_enumeration_plus_original_positives`: 3274
- `stage1_consensus_potency_and_uncertainty`: 2584
- `stage2_pains_brenk_free`: 1664
- `stage3_lipinski_psa_surface_area_flexibility`: 172
- `stage4_admet_ai_multigate`: 34
- `stage5_multistructure_template_docking_strict`: 4
- `stage5b_multistructure_template_docking_relaxed`: 5

## Potency modeling note

The potency model was treated conservatively. Leave-one-compound-out performance on the deduplicated dataset is modest, so consensus potency was used as one screen among several rather than as decisive proof of activity.

## Structural references

- `6DXT` small-molecule oxadiazole acid binder.
- `7MS5` difluorophenyl-piperidine sulfonamide keto-acid binder.
- `7MS6` fluorobenzoyl-glycine sulfonamide binder.
- `7MS7` chlorophenyl-piperidine sulfonamide glycine binder.

## Primary leads

- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.5224 | consensus pIC50 4.736 | potency std 0.453 | AMES 0.080 | hERG 0.213 | best binding 0.782 | best shape 0.692 | best pharm 0.917 | scaffold `O=S(=O)(c1cccnc1)N1CCC(c2ccccc2)CC1` | methods `original_positive`
- `O=C(O)CNC(=O)c1ncc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1` | score 0.5174 | consensus pIC50 4.481 | potency std 0.290 | AMES 0.054 | hERG 0.106 | best binding 0.830 | best shape 0.772 | best pharm 0.917 | scaffold `O=S(=O)(c1cncnc1)N1CCC(c2ccccc2)CC1` | methods `heteroatom_walks`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(F)cc3)CC2)cn1` | score 0.4467 | consensus pIC50 4.548 | potency std 0.318 | AMES 0.098 | hERG 0.123 | best binding 0.780 | best shape 0.681 | best pharm 0.972 | scaffold `O=S(=O)(c1cccnc1)N1CCC(c2ccccc2)CC1` | methods `matched_molecular_pair_expansion`
- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cn3)CC2)cn1` | score 0.3468 | consensus pIC50 4.481 | potency std 0.289 | AMES 0.036 | hERG 0.131 | best binding 0.762 | best shape 0.640 | best pharm 0.944 | scaffold `O=S(=O)(c1cccnc1)N1CCC(c2ccccn2)CC1` | methods `heteroatom_walks`

## Diverse backup leads

- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | score 0.5310 | consensus pIC50 4.778 | AMES 0.128 | hERG 0.225 | best binding 0.350 | scaffold `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` | primary parent `CHEMBL5410606`
- `CC(=O)C1=C(O)C(=O)N(c2ccc(O)cc2)C1c1ccc(Cl)cc1` | score 0.5282 | consensus pIC50 4.683 | AMES 0.093 | hERG 0.224 | best binding 0.355 | scaffold `O=C1C=CC(c2ccccc2)N1c1ccccc1` | primary parent `CHEMBL4129140`
- `CC(=O)C1=C(O)C(=O)N(c2ccccc2)C1c1ccc(Cl)cc1` | score 0.5142 | consensus pIC50 4.667 | AMES 0.044 | hERG 0.176 | best binding 0.364 | scaffold `O=C1C=CC(c2ccccc2)N1c1ccccc1` | primary parent `CHEMBL4129140`
- `COC1C(C(C)=O)=C(O)C(=O)N1C1C(C(C)=O)=C(O)C(=O)N1c1ccc(Cl)cc1` | score 0.4990 | consensus pIC50 4.521 | AMES 0.034 | hERG 0.035 | best binding 0.338 | scaffold `O=C1C=CCN1C1C=CC(=O)N1c1ccccc1` | primary parent `CHEMBL4129140`
- `CC(=O)C1=C(O)C(=O)N(c2ccc(Cl)cc2)C1N1CCCC1` | score 0.4844 | consensus pIC50 4.586 | AMES 0.101 | hERG 0.292 | best binding 0.475 | scaffold `O=C1C=CC(N2CCCC2)N1c1ccccc1` | primary parent `CHEMBL4129140`
- `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cncc2c1=O` | score 0.4780 | consensus pIC50 4.559 | AMES 0.064 | hERG 0.179 | best binding 0.415 | scaffold `O=C(Cn1cnc2ccncc2c1=O)NCc1ccccc1` | primary parent `CHEMBL5410606`

## Limitations

- The potency models remain low-data models and should not be overinterpreted.
- The 3D stage is a multistructure template-docking surrogate, not a production docking campaign with force-field or MD refinement.
- If one chemotype dominates strict survivors, that should be treated as current evidence concentration, not proof that other scaffolds are inactive.

## Literature anchors

- Mann et al. 2019, J Med Chem. Discovery of Small Molecule Antagonists of the USP5 Zinc Finger Ubiquitin-Binding Domain. PMID `31663737`.
- Wang et al. 2024, Comput Biol Med. Structure-based virtual screening of novel USP5 inhibitors targeting the zinc finger ubiquitin-binding domain. PMID `38603896`.
- RCSB structures used: `6DXT`, `7MS5`, `7MS6`, `7MS7`.