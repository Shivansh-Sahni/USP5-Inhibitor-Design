from __future__ import annotations

from pathlib import Path

import joblib
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

import lead_selection_multistructure_common as mod


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"

COUNTS_PATH = OUTPUT_DIR / "lead_selection_counts.csv"
STRICT_POOL_PATH = OUTPUT_DIR / "lead_selection_strict_pool.csv"
RELAXED_POOL_PATH = OUTPUT_DIR / "lead_selection_relaxed_pool.csv"
PRIMARY_LEADS_PATH = OUTPUT_DIR / "final_leads.csv"
BACKUP_LEADS_PATH = OUTPUT_DIR / "backup_leads.csv"
SUMMARY_PATH = OUTPUT_DIR / "lead_selection_summary.md"
DEPRIORITIZED_ZNF_PATH = OUTPUT_DIR / "deprioritized_znf_like_pool.csv"

POTENCY_CUTOFF = 4.60
NOVELTY_SIMILARITY_CUTOFF = 0.85
ZNF_SIMILARITY_CUTOFF = 0.30
BLOCKED_PRIMARY_PARENTS = {"CHEMBL5278336"}
PRIMARY_PARENT_CAP = 2
BACKUP_PARENT_CAP = 2
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


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


def add_novelty_columns(df: pd.DataFrame) -> pd.DataFrame:
    modeling = pd.read_csv(OUTPUT_DIR / "modeling_dataset.csv")
    known_smiles = modeling["canonical_smiles"].tolist()
    known_fps = [FP_GEN.GetFingerprint(Chem.MolFromSmiles(smiles)) for smiles in known_smiles]

    znf_smiles = [
        item["template_smiles"]
        for item in mod.STRUCTURES
        if item["name"].startswith("7MS")
    ]
    znf_smiles.extend(
        modeling[modeling["representative_id"].isin(BLOCKED_PRIMARY_PARENTS)]["canonical_smiles"].tolist()
    )
    znf_fps = [FP_GEN.GetFingerprint(Chem.MolFromSmiles(smiles)) for smiles in znf_smiles]

    df = df.copy()
    max_known = []
    max_znf = []
    for row in df.itertuples(index=False):
        fp = row.fp
        max_known.append(max(DataStructs.TanimotoSimilarity(fp, ref) for ref in known_fps))
        max_znf.append(max(DataStructs.TanimotoSimilarity(fp, ref) for ref in znf_fps))
    df["max_similarity_to_existing"] = max_known
    df["max_similarity_to_znf_reference"] = max_znf
    df["novelty_score"] = (
        0.55 * mod.normalize_lower(df["max_similarity_to_existing"])
        + 0.45 * mod.normalize_lower(df["max_similarity_to_znf_reference"])
    )
    return df


def score_novelty_frames(*frames: pd.DataFrame) -> None:
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
        frame["orthogonal_composite_score"] = (
            0.35 * mod.normalize_higher(frame["pred_pIC50"])
            + 0.15 * mod.normalize_lower(frame["descriptor_ad_distance"])
            + 0.20 * frame["admet_score"]
            + 0.10 * frame["property_score"]
            + 0.20 * frame["novelty_score"]
        )


def main() -> None:
    candidates = mod.prepare_candidates()
    _, training_features = mod.prepare_training_set()
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

    stage5 = add_novelty_columns(stage4)
    deprioritized_znf = stage5[
        stage5["primary_parent_id"].isin(BLOCKED_PRIMARY_PARENTS)
        | stage5["is_original_positive"]
        | (stage5["max_similarity_to_znf_reference"] >= ZNF_SIMILARITY_CUTOFF)
    ].copy()

    strict = stage5[
        (~stage5["primary_parent_id"].isin(BLOCKED_PRIMARY_PARENTS))
        & (~stage5["is_original_positive"])
        & (stage5["max_similarity_to_znf_reference"] < ZNF_SIMILARITY_CUTOFF)
        & (stage5["max_similarity_to_existing"] <= NOVELTY_SIMILARITY_CUTOFF)
    ].copy()
    relaxed = stage5[
        (~stage5["primary_parent_id"].isin(BLOCKED_PRIMARY_PARENTS))
        & (stage5["max_similarity_to_znf_reference"] < 0.35)
    ].copy()
    stage_counts.append(("stage5_non_znf_and_not_original", int(strict["product_smiles"].nunique())))
    stage_counts.append(("stage5b_non_znf_relaxed_pool", int(relaxed["product_smiles"].nunique())))

    score_novelty_frames(stage5, strict, relaxed, deprioritized_znf)
    strict = strict.sort_values(
        ["orthogonal_composite_score", "novelty_score", "pred_pIC50"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    relaxed = relaxed.sort_values(
        ["orthogonal_composite_score", "novelty_score", "pred_pIC50"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    deprioritized_znf = deprioritized_znf.sort_values(
        ["pred_pIC50", "max_similarity_to_znf_reference"],
        ascending=[False, False],
    ).reset_index(drop=True)

    primary = mod.select_portfolio(
        strict,
        mod.PRIMARY_LIMIT,
        dominant_scaffold_cap=2,
        parent_cap=PRIMARY_PARENT_CAP,
    )
    backup_source = relaxed.copy()
    if not primary.empty:
        backup_source = backup_source[~backup_source["product_smiles"].isin(set(primary["product_smiles"]))].copy()
    backup = mod.select_portfolio(
        backup_source,
        mod.BACKUP_LIMIT,
        dominant_scaffold_cap=None,
        exclude_smiles=set(primary["product_smiles"]),
        parent_cap=BACKUP_PARENT_CAP,
    )

    pd.DataFrame(stage_counts, columns=["stage", "remaining_unique_compounds"]).to_csv(COUNTS_PATH, index=False)
    strict.drop(columns=["mol", "fp"], errors="ignore").to_csv(STRICT_POOL_PATH, index=False)
    relaxed.drop(columns=["mol", "fp"], errors="ignore").to_csv(RELAXED_POOL_PATH, index=False)
    primary.drop(columns=["mol", "fp"], errors="ignore").to_csv(PRIMARY_LEADS_PATH, index=False)
    backup.drop(columns=["mol", "fp"], errors="ignore").to_csv(BACKUP_LEADS_PATH, index=False)
    deprioritized_znf.drop(columns=["mol", "fp"], errors="ignore").to_csv(DEPRIORITIZED_ZNF_PATH, index=False)

    lines = [
        "# Restored Final-Model Lead Selection",
        "",
        "This run preserves the original saved ExtraTrees regression model but changes the final prioritization goal. The primary leads are now intentionally filtered away from the known ZnF-UBD-like chemistry and re-ranked for novelty relative to the existing molecules, while still preserving potency, ADMET, and lead-like property constraints.",
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
            f"- `{row.product_smiles}` | orthogonal score {row.orthogonal_composite_score:.4f} | predicted pIC50 {row.pred_pIC50:.3f} "
            f"(pred IC50 {row.pred_ic50_uM:.2f} uM) | AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | "
            f"max existing similarity {row.max_similarity_to_existing:.3f} | max ZnF-reference similarity {row.max_similarity_to_znf_reference:.3f} | "
            f"scaffold `{row.scaffold}` | parent `{row.primary_parent_id}`"
        )
    lines.extend(["", "## Orthogonal backups", ""])
    if backup.empty:
        lines.append("- No orthogonal backup leads survived after the non-ZnF / novelty re-ranking.")
    else:
        for row in backup.itertuples(index=False):
            lines.append(
                f"- `{row.product_smiles}` | orthogonal score {row.orthogonal_composite_score:.4f} | predicted pIC50 {row.pred_pIC50:.3f} | "
                f"AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | max existing similarity {row.max_similarity_to_existing:.3f} | "
                f"max ZnF-reference similarity {row.max_similarity_to_znf_reference:.3f} | scaffold `{row.scaffold}` | parent `{row.primary_parent_id}`"
            )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"The original final model was preserved intact. The key change is post-model: known ZnF-like chemistry is explicitly deprioritized using a reference-similarity filter (`< {ZNF_SIMILARITY_CUTOFF:.2f}` to ZnF templates), original known molecules are excluded from primary leads, and novelty versus the existing dataset is rewarded. This shifts the primary list toward orthogonal chemotypes with more distinct scaffolds and less dependence on the previously favored CHEMBL5278336-like series.",
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
        print(
            primary[
                [
                    "product_smiles",
                    "orthogonal_composite_score",
                    "pred_pIC50",
                    "max_similarity_to_existing",
                    "max_similarity_to_znf_reference",
                    "primary_parent_id",
                ]
            ].to_string(index=False)
        )
    print("\nbackup_leads")
    if backup.empty:
        print("none")
    else:
        print(
            backup[
                [
                    "product_smiles",
                    "orthogonal_composite_score",
                    "pred_pIC50",
                    "max_similarity_to_existing",
                    "max_similarity_to_znf_reference",
                    "primary_parent_id",
                ]
            ].to_string(index=False)
        )


if __name__ == "__main__":
    main()
