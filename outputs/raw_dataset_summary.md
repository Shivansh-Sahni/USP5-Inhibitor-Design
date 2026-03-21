# Raw Dataset Summary

This summary was extracted directly from the uploaded [`First.csv`](/Users/shivanshsahni/Documents/New%20project/data/raw/First.csv) using the standard library only, before RDKit-based cleaning.

- Total rows: 26
- Columns present: `id`, `pIC50`, `ic50`, `smiles`
- Missing SMILES rows: 0
- Unique source IDs: 24
- Unique raw SMILES strings: 24
- Rows with measured `ic50 > 0`: 16
- Rows with active/no-numeric-IC50 label `ic50 == 0`: 4
- Rows with inactive label `ic50 == -1`: 6
- `pIC50` range: 2.5 to 6.1
- Mean `pIC50`: 4.25
- Median `pIC50`: 4.605
- Mean measured-row `pIC50`: 4.782
- Mean assigned-label-row `pIC50`: 3.4

One duplicate raw SMILES was already visible before canonicalization:

- `O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1`
  - `CHEMBL5278336`, `pIC50 = 6.1`, `ic50 = 800`
  - `CHEMBL5278336`, `pIC50 = 4.34`, `ic50 = 46000`
  - `CHEMBL5278336`, `pIC50 = 5.15`, `ic50 = 7000`

This confirms that duplicate handling must be explicit and that repeated molecules may have conflicting potency reports.
