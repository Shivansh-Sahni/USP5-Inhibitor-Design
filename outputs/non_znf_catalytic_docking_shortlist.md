# Non-ZnF Catalytic Docking Shortlist

This shortlist was extracted from `outputs/lead_selection_relaxed_pool.csv` for use in a catalytic-pocket docking campaign against USP5 where the zinc-finger hypothesis is intentionally excluded.

Ranking logic:

- prioritize higher `orthogonal_composite_score`
- prefer lower `max_similarity_to_znf_reference`
- keep predicted potency as a secondary tiebreaker

Interpretation:

- `strict`: lower ZnF-reference similarity (`<= 0.225`)
- `borderline`: still potentially useful for catalytic-pocket docking, but less cleanly separated from ZnF-like chemistry

## Recommended first-pass docking set

If you do not want to dock everything, start with these 10:

1. `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O`
2. `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(F)cccc2c1=O`
3. `Cc1cc(-n2c(C)cc(CCC(=O)O)c2C)c(C)n1CCC(=O)O`
4. `Cc1cc(-c2cc(C)n(CCC(=O)O)c2C)c(C)n1CCC(=O)O`
5. `CC(=O)C1=C(O)C(=O)N(c2ccc(O)cc2)C1c1ccc(Cl)cc1`
6. `CC(=O)C1=C(O)C(=O)N(c2ccc(F)cc2)C1c1ccc(Cl)cc1`
7. `COC1C(C(C)=O)=C(O)C(=O)N1C1C(C(C)=O)=C(O)C(=O)N1c1ccccc1`
8. `COC1C(C(C)=O)=C(O)C(=O)N1C1C(C(C)=O)=C(O)C(=O)N1c1ccc(F)cc1`
9. `CC(=O)C1=C(O)C(=O)N(c2ccc(Cl)cc2)C1N1CCCC1`
10. `O=C(O)CCc1ccccc1-c1ccccc1N1CCC(O)CC1`

This set keeps multiple chemotypes instead of overloading one series.

## Strict shortlist

| Rank | Product SMILES | Parent | Orthogonal score | Max sim to ZnF ref | Pred pIC50 | Scaffold |
| --- | --- | --- | ---: | ---: | ---: | --- |
| 1 | `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | `CHEMBL5410606` | 0.7123 | 0.1868 | 5.1000 | `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` |
| 2 | `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(F)cccc2c1=O` | `CHEMBL5410606` | 0.6542 | 0.2022 | 5.0192 | `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccc1` |
| 3 | `Cc1cc(-n2c(C)cc(CCC(=O)O)c2C)c(C)n1CCC(=O)O` | `[9*]n1c(C)cc([16*])c1C` | 0.6082 | 0.1857 | 4.8758 | `c1ccn(-c2cc[nH]c2)c1` |
| 4 | `Cc1cc(-c2cc(C)n(CCC(=O)O)c2C)c(C)n1CCC(=O)O` | `[9*]n1c(C)cc([16*])c1C` | 0.5901 | 0.2167 | 4.8743 | `c1cc(-c2cc[nH]c2)c[nH]1` |
| 5 | `COc1ccncc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | `CHEMBL5410606` | 0.5724 | 0.1915 | 4.8761 | `O=C(Cn1cnc2ccccc2c1=O)NCc1cccnc1` |
| 6 | `COc1ncccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | `CHEMBL5410606` | 0.5719 | 0.1957 | 4.8761 | `O=C(Cn1cnc2ccccc2c1=O)NCc1cccnc1` |
| 7 | `COc1cccnc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | `CHEMBL5410606` | 0.5719 | 0.1935 | 4.8761 | `O=C(Cn1cnc2ccccc2c1=O)NCc1ccccn1` |
| 8 | `COc1cnccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cccc2c1=O` | `CHEMBL5410606` | 0.5594 | 0.1915 | 4.8761 | `O=C(Cn1cnc2ccccc2c1=O)NCc1ccncc1` |
| 9 | `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)cncc2c1=O` | `CHEMBL5410606` | 0.5564 | 0.1957 | 4.8769 | `O=C(Cn1cnc2ccncc2c1=O)NCc1ccccc1` |
| 10 | `CC(=O)C1=C(O)C(=O)N(c2ccc(O)cc2)C1c1ccc(Cl)cc1` | `CHEMBL4129140` | 0.5313 | 0.2254 | 4.8777 | `O=C1C=CC(c2ccccc2)N1c1ccccc1` |
| 11 | `COc1ccccc1CNC(=O)Cn1c(CCC(=O)O)nc2c(Cl)ccnc2c1=O` | `CHEMBL5410606` | 0.5256 | 0.1978 | 4.8761 | `O=C(Cn1cnc2cccnc2c1=O)NCc1ccccc1` |
| 12 | `CC(=O)C1=C(O)C(=O)N(c2ccc(F)cc2)C1c1ccc(Cl)cc1` | `CHEMBL4129140` | 0.5054 | 0.2192 | 4.8347 | `O=C1C=CC(c2ccccc2)N1c1ccccc1` |
| 13 | `CC(=O)C1=C(O)C(=O)N(c2ccc(Cl)cc2)C1N1CCCC1` | `CHEMBL4129140` | 0.4226 | 0.2055 | 4.7226 | `O=C1C=CC(N2CCCC2)N1c1ccccc1` |
| 14 | `COC1C(C(C)=O)=C(O)C(=O)N1C1C(C(C)=O)=C(O)C(=O)N1c1ccccc1` | `CHEMBL4129140` | 0.4211 | 0.1447 | 4.6387 | `O=C1C=CCN1C1C=CC(=O)N1c1ccccc1` |
| 15 | `COC1C(C(C)=O)=C(O)C(=O)N1C1C(C(C)=O)=C(O)C(=O)N1c1ccc(F)cc1` | `CHEMBL4129140` | 0.4187 | 0.1392 | 4.6140 | `O=C1C=CCN1C1C=CC(=O)N1c1ccccc1` |

## Borderline additions

These are the next candidates if you want to push toward a fuller set.

| Rank | Product SMILES | Parent | Orthogonal score | Max sim to ZnF ref | Pred pIC50 |
| --- | --- | --- | ---: | ---: | ---: |
| 16 | `CC(=O)C1=C(O)C(=O)N([C@@H]2CCNC[C@@H]2C)C1c1ccc(Cl)cc1` | `CHEMBL4129140` | 0.4388 | 0.2308 | 4.7192 |
| 17 | `O=C(O)CCc1ccccc1-c1ccccc1N1CCC(O)CC1` | `CHEMBL5410606` | 0.5254 | 0.2319 | 4.7378 |
| 18 | `CC(=O)C1=C(O)C(=O)N(c2ccccc2)C1c1ccc(Cl)cc1` | `CHEMBL4129140` | 0.5152 | 0.2429 | 4.8639 |
| 19 | `CC(=O)C1=C(O)C(=O)N(CCC(=O)O)C1c1ccc(Cl)cc1` | `CHEMBL4129140` | 0.5543 | 0.2571 | 4.9085 |

Only 19 molecules in the current relaxed survivor pool were available for this non-ZnF-oriented ranking. There is no 20th candidate without relaxing the filter further.
