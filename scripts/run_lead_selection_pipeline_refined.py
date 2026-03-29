from __future__ import annotations

import csv
from pathlib import Path

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


def load_parent_meta() -> dict[str, dict[str, object]]:
    meta: dict[str, dict[str, object]] = {}
    for row in csv.DictReader(MODELING_PATH.open()):
        if row["has_measured_row"] == "True" or row["has_active_no_ic50_row"] == "True":
            meta[row["representative_id"]] = {
                "parent_pIC50": float(row["target_pIC50"]),
                "target_origin": row["target_origin"],
                "has_measured_row": row["has_measured_row"] == "True",
                "duplicate_conflict": row["duplicate_group_has_conflict"] == "True",
            }
    return meta


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
    parent_meta = load_parent_meta()
    alert_catalog = build_alert_catalog()
    rows = []
    for row in csv.DictReader(LIBRARY_PATH.open()):
        mol = Chem.MolFromSmiles(row["product_smiles"])
        if mol is None:
            continue
        meta = parent_meta[row["nearest_original_id"]]
        rows.append(
            {
                **row,
                **meta,
                "mol": mol,
                "fp": FP_GEN.GetFingerprint(mol),
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
                "sim": float(row["max_tanimoto_to_any_original"]),
            }
        )
    return pd.DataFrame(rows)


def run_pipeline(df: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    stages: list[tuple[str, int]] = [("start_similarity_window_library", int(df["product_smiles"].nunique()))]

    stage1 = df[df["alert_count"] == 0].copy()
    stages.append(("stage1_alert_free_pains_brenk_nih", int(stage1["product_smiles"].nunique())))

    stage2 = stage1[
        (stage1["has_measured_row"])
        & (stage1["parent_pIC50"] >= 4.5)
    ].copy()
    stages.append(("stage2_measured_parent_and_parent_pIC50_ge_4p5", int(stage2["product_smiles"].nunique())))

    stage3 = stage2[
        (stage2["mw"] <= 500)
        & stage2["logp"].between(1.0, 4.5)
        & stage2["tpsa"].between(50.0, 130.0)
        & (stage2["hbd"] <= 4)
        & (stage2["hba"] <= 9)
    ].copy()
    stages.append(("stage3_lead_like_property_window", int(stage3["product_smiles"].nunique())))

    stage4 = stage3[
        (stage3["rotb"] <= 8)
        & (stage3["formal_charge"].abs() <= 1)
        & (stage3["qed"] >= 0.50)
        & (stage3["sa_score"] <= 4.0)
    ].copy()
    stages.append(("stage4_developability_admet_proxies", int(stage4["product_smiles"].nunique())))

    stage5 = stage4[stage4["sim"].between(0.65, 0.88)].copy()
    stages.append(("stage5_similarity_sweet_spot", int(stage5["product_smiles"].nunique())))

    return stage5, stages


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    ranked = df.copy()
    ranked["potency_score"] = ranked["parent_pIC50"].apply(lambda x: closeness(x, 5.1, 0.8))
    ranked["sim_score"] = ranked["sim"].apply(lambda x: closeness(x, 0.75, 0.13))
    ranked["qed_score"] = ranked["qed"]
    ranked["sa_norm"] = ranked["sa_score"].apply(lambda x: closeness(x, 2.7, 1.5))
    ranked["mw_score"] = ranked["mw"].apply(lambda x: closeness(x, 420.0, 120.0))
    ranked["logp_score"] = ranked["logp"].apply(lambda x: closeness(x, 2.8, 1.7))
    ranked["tpsa_score"] = ranked["tpsa"].apply(lambda x: closeness(x, 95.0, 35.0))
    ranked["rot_score"] = ranked["rotb"].apply(lambda x: closeness(x, 6.0, 4.0))
    ranked["conflict_penalty"] = ranked["duplicate_conflict"].map({True: -0.05, False: 0.0})

    ranked["composite_score"] = (
        0.22 * ranked["potency_score"]
        + 0.22 * ranked["sim_score"]
        + 0.18 * ranked["qed_score"]
        + 0.12 * ranked["sa_norm"]
        + 0.10 * ranked["mw_score"]
        + 0.06 * ranked["logp_score"]
        + 0.05 * ranked["tpsa_score"]
        + 0.05 * ranked["rot_score"]
        + ranked["conflict_penalty"]
    )

    ranked = ranked.sort_values(
        ["composite_score", "qed", "sim", "parent_pIC50"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    return ranked


def diversity_pick(df: pd.DataFrame, n: int = 10, sim_ceiling: float = 0.72, per_parent_cap: int = 3) -> pd.DataFrame:
    selected = []
    selected_fps = []
    parent_counts: dict[str, int] = {}

    for row in df.itertuples(index=False):
        if parent_counts.get(row.nearest_original_id, 0) >= per_parent_cap:
            continue
        max_sim_to_selected = max(
            (DataStructs.TanimotoSimilarity(row.fp, prev_fp) for prev_fp in selected_fps),
            default=0.0,
        )
        need_new_parent = row.nearest_original_id not in parent_counts
        if len(selected) < 4 or max_sim_to_selected < sim_ceiling or need_new_parent:
            selected.append(row)
            selected_fps.append(row.fp)
            parent_counts[row.nearest_original_id] = parent_counts.get(row.nearest_original_id, 0) + 1
        if len(selected) >= n:
            break

    if len(selected) < n:
        existing = {row.product_smiles for row in selected}
        for row in df.itertuples(index=False):
            if row.product_smiles in existing:
                continue
            if parent_counts.get(row.nearest_original_id, 0) >= per_parent_cap:
                continue
            selected.append(row)
            parent_counts[row.nearest_original_id] = parent_counts.get(row.nearest_original_id, 0) + 1
            if len(selected) >= n:
                break

    return pd.DataFrame(selected)


def save_outputs(final_pool: pd.DataFrame, stages: list[tuple[str, int]], leads: pd.DataFrame) -> None:
    counts_path = OUTPUT_DIR / "lead_selection_refined_counts.csv"
    pool_path = OUTPUT_DIR / "lead_selection_refined_pool.csv"
    leads_path = OUTPUT_DIR / "final_leads_refined.csv"
    summary_path = OUTPUT_DIR / "lead_selection_refined_summary.md"

    pd.DataFrame(stages, columns=["stage", "remaining_unique_compounds"]).to_csv(counts_path, index=False)
    final_pool.drop(columns=["mol", "fp"], errors="ignore").to_csv(pool_path, index=False)
    leads.drop(columns=["mol", "fp"], errors="ignore").to_csv(leads_path, index=False)

    lines = [
        "# Refined Lead Selection Summary",
        "",
        "This refined funnel is stricter than the first pass and is designed to behave more like a real early lead-selection workflow.",
        "",
        "## Stage counts",
        "",
    ]
    for stage, count in stages:
        lines.append(f"- `{stage}`: {count}")
    lines.extend(
        [
            "",
            "## Stage logic",
            "",
            "- Stage 1: remove PAINS, BRENK, and NIH alerts.",
            "- Stage 2: require the nearest original parent to be a measured positive with pIC50 >= 4.5.",
            "- Stage 3: tighter lead-like property window.",
            "- Stage 4: developability/ADMET proxy filters using flexibility, charge, QED, and SA score.",
            "- Stage 5: keep compounds in a focused similarity sweet spot (0.65 to 0.88).",
            "",
            "## Final leads",
            "",
        ]
    )
    for row in leads.itertuples(index=False):
        lines.append(
            f"- `{row.product_smiles}` | score {row.composite_score:.4f} | nearest parent `{row.nearest_original_id}` "
            f"(pIC50 {row.parent_pIC50:.2f}) | sim {row.sim:.3f} | MW {row.mw:.1f} | cLogP {row.logp:.2f} | "
            f"TPSA {row.tpsa:.1f} | QED {row.qed:.3f} | SA {row.sa_score:.2f} | transform `{row.transform}`"
        )
    summary_path.write_text("\n".join(lines))


def main() -> None:
    df = prepare_dataframe()
    final_pool, stages = run_pipeline(df)
    ranked = rank_candidates(final_pool)
    leads = diversity_pick(ranked, n=10, sim_ceiling=0.72, per_parent_cap=3)
    save_outputs(ranked, stages, leads)

    print("stage_counts")
    for stage, count in stages:
        print(stage, count)
    print("\nfinal_leads_refined")
    cols = ["product_smiles", "nearest_original_id", "sim", "mw", "logp", "tpsa", "qed", "sa_score", "composite_score", "transform"]
    print(leads[cols].to_string(index=False))


if __name__ == "__main__":
    main()
