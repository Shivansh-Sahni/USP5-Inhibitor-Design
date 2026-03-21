from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


DESCRIPTOR_FUNCTIONS = {
    "molecular_weight": Descriptors.MolWt,
    "logp": Descriptors.MolLogP,
    "tpsa": rdMolDescriptors.CalcTPSA,
    "hbd": Lipinski.NumHDonors,
    "hba": Lipinski.NumHAcceptors,
    "rotatable_bonds": Lipinski.NumRotatableBonds,
    "ring_count": Lipinski.RingCount,
    "heavy_atom_count": Lipinski.HeavyAtomCount,
    "fraction_csp3": Lipinski.FractionCSP3,
}


def _morgan_bitvector(mol: Chem.Mol, radius: int, n_bits: int) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    bitvector = generator.GetFingerprint(mol)
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(bitvector, arr)
    return arr


def build_descriptor_table(modeling_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in modeling_df.iterrows():
        feature_row = {"canonical_smiles": row["canonical_smiles"]}
        for name, func in DESCRIPTOR_FUNCTIONS.items():
            feature_row[name] = float(func(row["mol"]))
        rows.append(feature_row)
    return pd.DataFrame(rows)


def build_fingerprint_table(modeling_df: pd.DataFrame, radius: int = 2, n_bits: int = 1024) -> pd.DataFrame:
    rows = []
    for _, row in modeling_df.iterrows():
        fp = _morgan_bitvector(row["mol"], radius=radius, n_bits=n_bits)
        fp_row = {"canonical_smiles": row["canonical_smiles"]}
        fp_row.update({f"fp_{index}": int(value) for index, value in enumerate(fp)})
        rows.append(fp_row)
    return pd.DataFrame(rows)


def build_similarity_outputs(modeling_df: pd.DataFrame, radius: int = 2, n_bits: int = 1024) -> tuple[pd.DataFrame, pd.DataFrame]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    bitvectors = []
    labels = modeling_df["canonical_smiles"].tolist()
    ids = modeling_df["representative_id"].tolist()

    for mol in modeling_df["mol"]:
        bitvectors.append(generator.GetFingerprint(mol))

    matrix = np.eye(len(bitvectors))
    for i, j in combinations(range(len(bitvectors)), 2):
        score = DataStructs.TanimotoSimilarity(bitvectors[i], bitvectors[j])
        matrix[i, j] = score
        matrix[j, i] = score

    similarity_matrix = pd.DataFrame(matrix, index=labels, columns=labels)

    summaries = []
    for idx, row in enumerate(matrix):
        others = np.delete(row, idx)
        if len(others) == 0:
            max_similarity = np.nan
            mean_similarity = np.nan
            nearest_index = None
        else:
            nearest_index = np.argsort(row)[-2]
            max_similarity = float(row[nearest_index])
            mean_similarity = float(np.mean(others))

        summaries.append(
            {
                "canonical_smiles": labels[idx],
                "representative_id": ids[idx],
                "max_tanimoto_to_other": max_similarity,
                "mean_tanimoto_to_others": mean_similarity,
                "nearest_neighbor_smiles": None if nearest_index is None else labels[nearest_index],
                "nearest_neighbor_id": None if nearest_index is None else ids[nearest_index],
            }
        )

    return similarity_matrix, pd.DataFrame(summaries)


def build_scaffold_summary(modeling_df: pd.DataFrame) -> pd.DataFrame:
    scaffolds = []
    for _, row in modeling_df.iterrows():
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=row["mol"])
        scaffolds.append({"canonical_smiles": row["canonical_smiles"], "murcko_scaffold": scaffold})

    scaffold_df = pd.DataFrame(scaffolds)
    counts = scaffold_df["murcko_scaffold"].value_counts(dropna=False).rename_axis("murcko_scaffold").reset_index(name="count")
    return scaffold_df.merge(counts, on="murcko_scaffold", how="left")
