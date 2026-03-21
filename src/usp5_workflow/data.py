from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem


REQUIRED_COLUMNS = ["id", "pIC50", "ic50", "smiles"]


@dataclass
class DatasetBundle:
    raw: pd.DataFrame
    annotated: pd.DataFrame
    modeling: pd.DataFrame


def load_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df.copy()


def annotate_labels(df: pd.DataFrame) -> pd.DataFrame:
    annotated = df.copy()
    annotated["source_id"] = annotated["id"].astype(str)
    annotated["pIC50"] = pd.to_numeric(annotated["pIC50"], errors="coerce")
    annotated["ic50"] = pd.to_numeric(annotated["ic50"], errors="coerce")
    annotated["smiles"] = annotated["smiles"].astype(str).str.strip()

    annotated["is_measured"] = annotated["ic50"] > 0
    annotated["is_active_no_ic50"] = annotated["ic50"] == 0
    annotated["is_inactive"] = annotated["ic50"] == -1
    annotated["is_assigned_label"] = annotated["ic50"].isin([0, -1])

    return annotated


def validate_and_canonicalize_smiles(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df.copy()
    mols = []
    canonical_smiles = []
    valid_flags = []

    for smiles in cleaned["smiles"]:
        mol = Chem.MolFromSmiles(smiles)
        mols.append(mol)
        valid = mol is not None
        valid_flags.append(valid)
        canonical_smiles.append(Chem.MolToSmiles(mol) if valid else None)

    cleaned["is_valid_smiles"] = valid_flags
    cleaned["canonical_smiles"] = canonical_smiles
    cleaned["mol"] = mols
    return cleaned


def _aggregate_duplicate_group(group: pd.DataFrame) -> pd.Series:
    measured = group[group["is_measured"]]
    chosen = measured if not measured.empty else group
    modeling_target = chosen["pIC50"].median()

    if not measured.empty and group["is_assigned_label"].any():
        target_origin = "mixed_duplicate_measured_median"
    elif not measured.empty:
        target_origin = "measured_median"
    else:
        target_origin = "assigned_median"

    return pd.Series(
        {
            "canonical_smiles": group.name,
            "compound_ids": ";".join(group["source_id"].astype(str)),
            "representative_id": str(group["source_id"].iloc[0]),
            "row_count": int(len(group)),
            "duplicate_group_has_conflict": bool(group["pIC50"].nunique() > 1),
            "pIC50_values": ";".join(group["pIC50"].map(lambda x: f"{x:.4g}")),
            "ic50_values": ";".join(group["ic50"].map(lambda x: f"{x:.4g}")),
            "target_pIC50": float(modeling_target),
            "target_origin": target_origin,
            "n_measured_rows": int(group["is_measured"].sum()),
            "n_assigned_rows": int(group["is_assigned_label"].sum()),
            "has_measured_row": bool(group["is_measured"].any()),
            "has_active_no_ic50_row": bool(group["is_active_no_ic50"].any()),
            "has_inactive_row": bool(group["is_inactive"].any()),
            "mol": group["mol"].iloc[0],
        }
    )


def deduplicate_for_modeling(df: pd.DataFrame) -> pd.DataFrame:
    valid = df[df["is_valid_smiles"]].copy()
    grouped_rows = []
    for canonical_smiles, group in valid.groupby("canonical_smiles", sort=False, dropna=False):
        group = group.copy()
        group.name = canonical_smiles
        grouped_rows.append(_aggregate_duplicate_group(group))
    modeling = pd.DataFrame(grouped_rows)
    return modeling


def prepare_dataset(path: Path) -> DatasetBundle:
    raw = load_dataset(path)
    annotated = annotate_labels(raw)
    annotated = validate_and_canonicalize_smiles(annotated)

    invalid_pIC50 = annotated["pIC50"].isna().sum()
    invalid_ic50 = annotated["ic50"].isna().sum()
    if invalid_pIC50 or invalid_ic50:
        raise ValueError(
            f"Numeric parsing failed for {invalid_pIC50} pIC50 values and {invalid_ic50} ic50 values."
        )

    modeling = deduplicate_for_modeling(annotated)
    return DatasetBundle(raw=raw, annotated=annotated, modeling=modeling)
