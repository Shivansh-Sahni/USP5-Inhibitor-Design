from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from usp5_workflow.data import prepare_dataset
from usp5_workflow.features import (
    build_descriptor_table,
    build_fingerprint_table,
    build_scaffold_summary,
    build_similarity_outputs,
)
from usp5_workflow.modeling import build_ranked_compounds, run_loocv


def _write_summary(
    output_dir: Path,
    annotated: pd.DataFrame,
    modeling: pd.DataFrame,
    similarity_summary: pd.DataFrame,
    scaffold_summary: pd.DataFrame,
    metrics: pd.DataFrame,
) -> None:
    total_rows = len(annotated)
    valid_rows = int(annotated["is_valid_smiles"].sum())
    unique_molecules = len(modeling)
    duplicate_groups = int((modeling["row_count"] > 1).sum())
    measured_rows = int(annotated["is_measured"].sum())
    assigned_rows = int(annotated["is_assigned_label"].sum())
    duplicate_conflicts = int(modeling["duplicate_group_has_conflict"].sum())
    top_scaffolds = scaffold_summary["murcko_scaffold"].nunique(dropna=False)

    tanimoto_values = similarity_summary["max_tanimoto_to_other"].dropna()
    best_model = metrics.iloc[0]

    lines = [
        "# USP5 Baseline Analysis Summary",
        "",
        "## Dataset",
        f"- Total input rows: {total_rows}",
        f"- Valid SMILES rows: {valid_rows}",
        f"- Unique molecules after canonical deduplication: {unique_molecules}",
        f"- Measured rows (`ic50 > 0`): {measured_rows}",
        f"- Assigned-label rows (`ic50 in {{0, -1}}`): {assigned_rows}",
        f"- Duplicate molecule groups: {duplicate_groups}",
        f"- Duplicate groups with conflicting pIC50 values: {duplicate_conflicts}",
        "",
        "## Chemical diversity",
        f"- Unique Bemis-Murcko scaffolds: {top_scaffolds}",
        f"- Mean nearest-neighbor Tanimoto: {tanimoto_values.mean():.3f}" if not tanimoto_values.empty else "- Mean nearest-neighbor Tanimoto: not available",
        f"- Max nearest-neighbor Tanimoto: {tanimoto_values.max():.3f}" if not tanimoto_values.empty else "- Max nearest-neighbor Tanimoto: not available",
        "",
        "## Modeling",
        f"- Best LOOCV model: {best_model['model']}",
        f"- Best feature set: {best_model['feature_set']}",
        f"- MAE: {best_model['mae']:.3f}",
        f"- RMSE: {best_model['rmse']:.3f}",
        f"- R2: {best_model['r2']:.3f}",
        "",
        "## Interpretation",
        "- This is a very small dataset, so LOOCV estimates can still be unstable.",
        "- The best setup was selected by comparing many feature and model configurations on the same dataset, so the top score is optimistic.",
        "- Assigned labels are mixed with measured potency values, so predictive performance should be treated as provisional.",
        "- Duplicate molecules with inconsistent activities suggest assay or curation variability and should be reviewed before trusting fine-grained rankings.",
        "- The next rational improvement is a measured-only sensitivity analysis and more experimentally measured compounds, not a more complex model.",
    ]

    (output_dir / "analysis_summary.md").write_text("\n".join(lines))


def run_workflow(
    input_path: Path,
    output_dir: Path,
    fingerprint_sizes: list[int],
    random_seed: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset = prepare_dataset(input_path)
    descriptor_df = build_descriptor_table(dataset.modeling)
    fingerprint_tables = {
        size: build_fingerprint_table(dataset.modeling, radius=2, n_bits=size) for size in fingerprint_sizes
    }
    similarity_matrix, similarity_summary = build_similarity_outputs(dataset.modeling, radius=2, n_bits=1024)
    scaffold_summary = build_scaffold_summary(dataset.modeling)
    feature_sets: dict[str, pd.DataFrame] = {
        "descriptors": descriptor_df,
        "similarity_summary": similarity_summary[
            ["canonical_smiles", "max_tanimoto_to_other", "mean_tanimoto_to_others"]
        ],
        "descriptors_plus_similarity": descriptor_df.merge(
            similarity_summary[["canonical_smiles", "max_tanimoto_to_other", "mean_tanimoto_to_others"]],
            on="canonical_smiles",
            how="left",
        ),
    }
    for size, fp_table in fingerprint_tables.items():
        feature_sets[f"morgan_{size}"] = fp_table
        feature_sets[f"descriptors_plus_morgan_{size}"] = descriptor_df.merge(
            fp_table,
            on="canonical_smiles",
            how="left",
        )
        feature_sets[f"descriptors_plus_similarity_plus_morgan_{size}"] = feature_sets[
            "descriptors_plus_similarity"
        ].merge(fp_table, on="canonical_smiles", how="left")

    model_result = run_loocv(dataset.modeling, feature_sets, random_seed=random_seed)
    ranked = build_ranked_compounds(dataset.modeling, model_result)

    dataset.annotated.drop(columns=["mol"]).to_csv(output_dir / "annotated_rows.csv", index=False)
    dataset.modeling.drop(columns=["mol"]).to_csv(output_dir / "modeling_dataset.csv", index=False)
    descriptor_df.to_csv(output_dir / "descriptor_features.csv", index=False)
    for feature_set_name, feature_df in feature_sets.items():
        safe_name = feature_set_name.replace("/", "_")
        feature_df.to_csv(output_dir / f"feature_set_{safe_name}.csv", index=False)
    for size, table in fingerprint_tables.items():
        table.to_csv(output_dir / f"fingerprint_features_{size}.csv", index=False)
    similarity_matrix.to_csv(output_dir / "tanimoto_similarity_matrix.csv")
    similarity_summary.to_csv(output_dir / "similarity_summary.csv", index=False)
    scaffold_summary.to_csv(output_dir / "scaffold_summary.csv", index=False)
    model_result.metrics.to_csv(output_dir / "model_metrics.csv", index=False)
    model_result.predictions.to_csv(output_dir / "loocv_predictions.csv", index=False)
    model_result.full_fit_predictions.to_csv(output_dir / "full_fit_predictions.csv", index=False)
    ranked.to_csv(output_dir / "ranked_compounds.csv", index=False)

    _write_summary(
        output_dir=output_dir,
        annotated=dataset.annotated,
        modeling=dataset.modeling,
        similarity_summary=similarity_summary,
        scaffold_summary=scaffold_summary,
        metrics=model_result.metrics,
    )
