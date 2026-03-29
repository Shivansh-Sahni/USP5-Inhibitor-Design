from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Contrib.SA_Score import sascorer


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
LIBRARY_PATH = OUTPUT_DIR / "enumeration_similarity_window_0p55_0p95.csv"
MODELING_PATH = OUTPUT_DIR / "modeling_dataset.csv"

RDLogger.DisableLog("rdApp.*")
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def load_parent_activity() -> dict[str, float]:
    rows = list(csv.DictReader(MODELING_PATH.open()))
    out = {}
    for row in rows:
        if row["has_measured_row"] == "True" or row["has_active_no_ic50_row"] == "True":
            out[row["representative_id"]] = float(row["target_pIC50"])
    return out


def build_alert_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
    return FilterCatalog(params)


def closeness(value: float, target: float, tolerance: float) -> float:
    score = 1.0 - abs(value - target) / tolerance
    return max(0.0, min(1.0, score))


def prepare_dataframe() -> pd.DataFrame:
    parent_activity = load_parent_activity()
    alert_catalog = build_alert_catalog()
    rows = []
    for row in csv.DictReader(LIBRARY_PATH.open()):
        mol = Chem.MolFromSmiles(row["product_smiles"])
        if mol is None:
            continue
        fp = FP_GEN.GetFingerprint(mol)
        rows.append(
            {
                **row,
                "mol": mol,
                "fp": fp,
                "mw": Descriptors.MolWt(mol),
                "logp": Descriptors.MolLogP(mol),
                "tpsa": rdMolDescriptors.CalcTPSA(mol),
                "hbd": Lipinski.NumHDonors(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "rotb": Lipinski.NumRotatableBonds(mol),
                "rings": Lipinski.RingCount(mol),
                "fsp3": Lipinski.FractionCSP3(mol),
                "formal_charge": Chem.GetFormalCharge(mol),
                "qed": QED.qed(mol),
                "sa_score": sascorer.calculateScore(mol),
                "alert_count": len(alert_catalog.GetMatches(mol)),
                "nearest_parent_pIC50": parent_activity.get(row["nearest_original_id"], np.nan),
                "sim": float(row["max_tanimoto_to_any_original"]),
            }
        )
    return pd.DataFrame(rows)


def run_pipeline(df: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    stages: list[tuple[str, int]] = [("start_similarity_window_library", int(df["product_smiles"].nunique()))]

    stage1 = df[df["alert_count"] == 0].copy()
    stages.append(("stage1_alert_free_pains_brenk_nih", int(stage1["product_smiles"].nunique())))

    stage2 = stage1[
        stage1["mw"].between(250, 550)
        & stage1["logp"].between(1, 5)
        & (stage1["hbd"] <= 5)
        & (stage1["hba"] <= 10)
    ].copy()
    stages.append(("stage2_lead_like_core_properties", int(stage2["product_smiles"].nunique())))

    stage3 = stage2[
        stage2["tpsa"].between(40, 140)
        & (stage2["rotb"] <= 10)
        & (stage2["rings"] <= 6)
    ].copy()
    stages.append(("stage3_veber_surface_area_flexibility", int(stage3["product_smiles"].nunique())))

    # Local ADMET proxies only: no dedicated hERG/Ames models available in the repo.
    stage4 = stage3[
        (stage3["qed"] >= 0.45)
        & (stage3["formal_charge"].abs() <= 1)
        & stage3["sim"].between(0.60, 0.90)
    ].copy()
    stages.append(("stage4_admet_proxy_qed_charge_similarity_focus", int(stage4["product_smiles"].nunique())))

    return stage4, stages


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    potency_min = df["nearest_parent_pIC50"].min()
    potency_max = df["nearest_parent_pIC50"].max()
    potency_span = max(0.001, potency_max - potency_min)

    ranked = df.copy()
    ranked["sim_score"] = ranked["sim"].apply(lambda x: closeness(x, 0.75, 0.20))
    ranked["potency_score"] = (ranked["nearest_parent_pIC50"] - potency_min) / potency_span
    ranked["qed_score"] = ranked["qed"]
    ranked["sa_norm"] = ranked["sa_score"].apply(lambda x: closeness(x, 3.0, 3.5))
    ranked["logp_score"] = ranked["logp"].apply(lambda x: closeness(x, 2.8, 2.2))
    ranked["tpsa_score"] = ranked["tpsa"].apply(lambda x: closeness(x, 85.0, 55.0))
    ranked["mw_score"] = ranked["mw"].apply(lambda x: closeness(x, 420.0, 170.0))

    ranked["composite_score"] = (
        0.25 * ranked["sim_score"]
        + 0.25 * ranked["potency_score"]
        + 0.20 * ranked["qed_score"]
        + 0.15 * ranked["sa_norm"]
        + 0.05 * ranked["logp_score"]
        + 0.05 * ranked["tpsa_score"]
        + 0.05 * ranked["mw_score"]
    )
    ranked = ranked.sort_values(
        ["composite_score", "qed", "sim", "nearest_parent_pIC50"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    return ranked


def diversity_pick(df: pd.DataFrame, n: int = 10, sim_ceiling: float = 0.75) -> pd.DataFrame:
    selected_rows = []
    selected_fps = []
    selected_parents = set()

    for row in df.itertuples(index=False):
        fp = row.fp
        max_sim = max((DataStructs.TanimotoSimilarity(fp, prev) for prev in selected_fps), default=0.0)
        parent_bonus = row.nearest_original_id not in selected_parents
        if len(selected_rows) < 3:
            selected_rows.append(row)
            selected_fps.append(fp)
            selected_parents.add(row.nearest_original_id)
        elif max_sim < sim_ceiling or parent_bonus:
            selected_rows.append(row)
            selected_fps.append(fp)
            selected_parents.add(row.nearest_original_id)
        if len(selected_rows) >= n:
            break

    if len(selected_rows) < n:
        existing = {row.product_smiles for row in selected_rows}
        for row in df.itertuples(index=False):
            if row.product_smiles in existing:
                continue
            selected_rows.append(row)
            if len(selected_rows) >= n:
                break

    return pd.DataFrame(selected_rows)


def save_outputs(final_pool: pd.DataFrame, stages: list[tuple[str, int]], leads: pd.DataFrame) -> None:
    pipeline_path = OUTPUT_DIR / "lead_selection_pipeline_counts.csv"
    final_pool_path = OUTPUT_DIR / "lead_selection_filtered_pool.csv"
    leads_path = OUTPUT_DIR / "final_leads.csv"
    summary_path = OUTPUT_DIR / "lead_selection_summary.md"

    pd.DataFrame(stages, columns=["stage", "remaining_unique_compounds"]).to_csv(pipeline_path, index=False)

    drop_cols = ["mol", "fp"]
    final_pool.drop(columns=drop_cols, errors="ignore").to_csv(final_pool_path, index=False)
    leads.drop(columns=drop_cols, errors="ignore").to_csv(leads_path, index=False)

    lines = [
        "# Lead Selection Summary",
        "",
        "Starting library: similarity-window enumeration (`0.55-0.95` vs original positives).",
        "",
        "## Stage counts",
        "",
    ]
    for stage, count in stages:
        lines.append(f"- `{stage}`: {count}")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Stage 1 removes PAINS, BRENK, and NIH alerts.",
            "- Stage 2 applies core lead-like property limits.",
            "- Stage 3 applies Veber/surface-area/flexibility limits.",
            "- Stage 4 is an ADMET proxy stage using QED, charge sanity, and a tighter similarity focus.",
            "- Dedicated hERG, Ames, or Kd models were not available locally, so these were not run as true predictive models.",
            "",
            "## Final leads",
            "",
        ]
    )
    for row in leads.itertuples(index=False):
        lines.append(
            f"- `{row.product_smiles}` | score {row.composite_score:.4f} | nearest parent `{row.nearest_original_id}` "
            f"(pIC50 {row.nearest_parent_pIC50:.2f}) | sim {row.sim:.3f} | MW {row.mw:.1f} | cLogP {row.logp:.2f} "
            f"| TPSA {row.tpsa:.1f} | QED {row.qed:.3f} | SA {row.sa_score:.2f} | transform `{row.transform}`"
        )
    summary_path.write_text("\n".join(lines))


def main() -> None:
    df = prepare_dataframe()
    final_pool, stages = run_pipeline(df)
    ranked = rank_candidates(final_pool)
    leads = diversity_pick(ranked, n=10, sim_ceiling=0.75)
    save_outputs(ranked, stages, leads)

    print("stage_counts")
    for stage, count in stages:
        print(stage, count)
    print("\nfinal_leads")
    cols = ["product_smiles", "nearest_original_id", "sim", "mw", "logp", "tpsa", "qed", "sa_score", "composite_score", "transform"]
    print(leads[cols].to_string(index=False))


if __name__ == "__main__":
    main()
