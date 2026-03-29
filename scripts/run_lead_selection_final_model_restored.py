from __future__ import annotations

from pathlib import Path

import joblib
import pandas as pd

import lead_selection_multistructure_common as mod


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"

COUNTS_PATH = OUTPUT_DIR / "lead_selection_final_model_restored_counts.csv"
STRICT_POOL_PATH = OUTPUT_DIR / "lead_selection_final_model_restored_strict_pool.csv"
RELAXED_POOL_PATH = OUTPUT_DIR / "lead_selection_final_model_restored_relaxed_pool.csv"
PRIMARY_LEADS_PATH = OUTPUT_DIR / "final_leads_final_model_restored.csv"
BACKUP_LEADS_PATH = OUTPUT_DIR / "backup_leads_final_model_restored.csv"
SUMMARY_PATH = OUTPUT_DIR / "lead_selection_final_model_restored_summary.md"

POTENCY_CUTOFF = 4.60


def add_final_model_predictions(df: pd.DataFrame, training_features: pd.DataFrame) -> pd.DataFrame:
    model = joblib.load(OUTPUT_DIR / "final_model" / "final_model.joblib")
    X = df[list(mod.DESCRIPTOR_FUNCS.keys())].copy()
    df = df.copy()
    df["pred_pIC50"] = model.predict(X)
    df["pred_ic50_uM"] = 10 ** (6 - df["pred_pIC50"])

    desc_train = training_features[list(mod.DESCRIPTOR_FUNCS.keys())]
    mean = desc_train.mean()
    std = desc_train.std(ddof=0).replace(0, 1.0)
    z = (X - mean) / std
    df["descriptor_ad_distance"] = ((z**2).mean(axis=1)) ** 0.5
    return df


def score_frames(*frames: pd.DataFrame) -> None:
    for frame in frames:
        if frame.empty:
            continue
        frame["admet_score"] = (
            mod.normalize_lower(frame["AMES"])
            + mod.normalize_lower(frame["hERG"])
            + mod.normalize_lower(frame["ClinTox"])
            + mod.normalize_lower(frame["DILI"])
            + mod.normalize_higher(frame["HIA_Hou"])
            + mod.normalize_higher(frame["Bioavailability_Ma"])
            + mod.normalize_higher(frame["Solubility_AqSolDB"])
        ) / 7.0
        frame["property_score"] = (
            mod.closeness(frame["logp"], 2.7, 1.8)
            + mod.closeness(frame["tpsa"], 95.0, 45.0)
            + mod.closeness(frame["labute_asa"], 175.0, 45.0)
            + mod.closeness(frame["rot"], 6.0, 4.0)
            + mod.normalize_higher(frame["qed"])
        ) / 5.0
        frame["binding_component"] = (
            0.45 * mod.normalize_higher(frame["best_binding_score"])
            + 0.25 * mod.normalize_higher(frame["mean_top2_binding_score"])
            + 0.20 * mod.normalize_higher(frame["best_shape_tanimoto"])
            + 0.10 * mod.normalize_higher(frame["best_pharmacophore_score"])
        )
        frame["restored_composite_score"] = (
            0.35 * mod.normalize_higher(frame["pred_pIC50"])
            + 0.10 * mod.normalize_lower(frame["descriptor_ad_distance"])
            + 0.20 * frame["admet_score"]
            + 0.10 * frame["property_score"]
            + 0.25 * frame["binding_component"]
        )


def main() -> None:
    candidates = mod.prepare_candidates()
    modeling, training_features = mod.prepare_training_set()
    screened = add_final_model_predictions(candidates, training_features)

    stage_counts: list[tuple[str, int]] = [("start_broad_enumeration_plus_original_positives", int(screened["product_smiles"].nunique()))]

    stage1 = screened[
        (screened["pred_pIC50"] >= POTENCY_CUTOFF)
        & (screened["descriptor_ad_distance"] <= 4.5)
    ].copy()
    stage_counts.append((f"stage1_original_final_model_predicted_pIC50_ge_{POTENCY_CUTOFF:.2f}", int(stage1["product_smiles"].nunique())))

    stage2 = stage1[stage1["pains_brenk_alerts"] == 0].copy()
    stage_counts.append(("stage2_pains_brenk_free", int(stage2["product_smiles"].nunique())))

    stage3 = stage2[
        (stage2["lipinski_violations"] <= 1)
        & stage2["tpsa"].between(40.0, 140.0)
        & stage2["labute_asa"].between(130.0, 235.0)
        & stage2["mw"].between(250.0, 550.0)
        & (stage2["rot"] <= 10)
        & (stage2["formal_charge"].abs() <= 1)
    ].copy()
    stage_counts.append(("stage3_lipinski_psa_surface_area_flexibility", int(stage3["product_smiles"].nunique())))

    stage4 = mod.run_admet(stage3)
    stage4 = stage4[
        (stage4["AMES"] <= 0.20)
        & (stage4["hERG"] <= 0.40)
        & (stage4["ClinTox"] <= 0.50)
        & (stage4["HIA_Hou"] >= 0.55)
        & (stage4["Bioavailability_Ma"] >= 0.45)
        & (stage4["Caco2_Wang"] >= -5.80)
        & (stage4["Solubility_AqSolDB"] >= -4.20)
    ].copy()
    stage_counts.append(("stage4_admet_ai_multigate", int(stage4["product_smiles"].nunique())))

    stage5 = mod.run_multistructure_scoring(stage4)
    strict = stage5[
        (stage5["best_binding_score"] >= 0.68)
        & (stage5["mean_top2_binding_score"] >= 0.56)
        & (stage5["best_shape_tanimoto"] >= 0.58)
        & (stage5["best_pharmacophore_score"] >= 0.50)
        & (stage5["min_clash_count"] <= 2)
    ].copy()
    relaxed = stage5[
        (stage5["best_binding_score"] >= 0.56)
        & (stage5["best_shape_tanimoto"] >= 0.48)
        & (stage5["best_pharmacophore_score"] >= 0.40)
        & (stage5["min_clash_count"] <= 3)
    ].copy()
    stage_counts.append(("stage5_multistructure_template_docking_strict", int(strict["product_smiles"].nunique())))
    stage_counts.append(("stage5b_multistructure_template_docking_relaxed", int(relaxed["product_smiles"].nunique())))

    score_frames(stage5, strict, relaxed)
    strict = strict.sort_values(["restored_composite_score", "best_binding_score", "pred_pIC50"], ascending=[False, False, False]).reset_index(drop=True)
    relaxed = relaxed.sort_values(["restored_composite_score", "best_binding_score", "pred_pIC50"], ascending=[False, False, False]).reset_index(drop=True)

    primary = mod.select_portfolio(strict, mod.PRIMARY_LIMIT, dominant_scaffold_cap=4)
    dominant_parent = primary["primary_parent_id"].mode().iat[0] if not primary.empty else None
    exploratory = stage5.copy()
    if dominant_parent is not None:
        exploratory = exploratory[exploratory["primary_parent_id"] != dominant_parent].copy()
    exploratory = exploratory.sort_values(
        ["restored_composite_score", "best_binding_score", "pred_pIC50"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    backup_source = exploratory if not exploratory.empty else relaxed
    backup = mod.select_portfolio(
        backup_source,
        mod.BACKUP_LIMIT,
        dominant_scaffold_cap=None,
        exclude_smiles=set(primary["product_smiles"]),
    )

    pd.DataFrame(stage_counts, columns=["stage", "remaining_unique_compounds"]).to_csv(COUNTS_PATH, index=False)
    strict.drop(columns=["mol", "fp"], errors="ignore").to_csv(STRICT_POOL_PATH, index=False)
    relaxed.drop(columns=["mol", "fp"], errors="ignore").to_csv(RELAXED_POOL_PATH, index=False)
    primary.drop(columns=["mol", "fp"], errors="ignore").to_csv(PRIMARY_LEADS_PATH, index=False)
    backup.drop(columns=["mol", "fp"], errors="ignore").to_csv(BACKUP_LEADS_PATH, index=False)

    lines = [
        "# Restored Final-Model Lead Selection",
        "",
        "This run restores the original saved ExtraTrees regression model as the only potency model. The downstream workflow is still the stronger version: broad library input, PAINS/BRENK cleanup, lead-like property filters, ADMET-AI, and multistructure USP5 3D binding plausibility.",
        "",
        "## Stage counts",
        "",
    ]
    for stage, count in stage_counts:
        lines.append(f"- `{stage}`: {count}")
    lines.extend(
        [
            "",
            "## Primary leads",
            "",
        ]
    )
    for row in primary.itertuples(index=False):
        lines.append(
            f"- `{row.product_smiles}` | score {row.restored_composite_score:.4f} | predicted pIC50 {row.pred_pIC50:.3f} "
            f"(pred IC50 {row.pred_ic50_uM:.2f} uM) | AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | "
            f"best binding {row.best_binding_score:.3f} | scaffold `{row.scaffold}` | parent `{row.primary_parent_id}`"
        )
    lines.extend(["", "## Orthogonal backups", ""])
    if backup.empty:
        lines.append("- No non-dominant backup leads survived after excluding the dominant parent series.")
    else:
        for row in backup.itertuples(index=False):
            lines.append(
                f"- `{row.product_smiles}` | score {row.restored_composite_score:.4f} | predicted pIC50 {row.pred_pIC50:.3f} | "
                f"AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | best binding {row.best_binding_score:.3f} | "
                f"scaffold `{row.scaffold}` | parent `{row.primary_parent_id}`"
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The original final model was preserved intact. The main scientific upgrades happen after potency scoring: better developability triage, ADMET-AI, and multistructure 3D evidence. This means the final list reflects the trusted model while still avoiding overreliance on a single chemotype in the final program recommendation.",
        ]
    )
    SUMMARY_PATH.write_text("\n".join(lines))

    print("stage_counts")
    for stage, count in stage_counts:
        print(stage, count)
    print("\nprimary_leads")
    if primary.empty:
        print("none")
    else:
        print(primary[["product_smiles", "restored_composite_score", "pred_pIC50", "best_binding_score", "primary_parent_id"]].to_string(index=False))
    print("\nbackup_leads")
    if backup.empty:
        print("none")
    else:
        print(backup[["product_smiles", "restored_composite_score", "pred_pIC50", "best_binding_score", "primary_parent_id"]].to_string(index=False))


if __name__ == "__main__":
    main()
