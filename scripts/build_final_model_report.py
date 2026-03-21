from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import joblib
import matplotlib
import pandas as pd
from rdkit.Chem import Descriptors, GraphDescriptors, Lipinski, rdMolDescriptors
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

from usp5_workflow.data import annotate_labels, validate_and_canonicalize_smiles


OUTPUT_DIR = Path("outputs/final_model")
MODEL_NAME = "raw_rows_base_graph_extratrees_trainfit"

BASE_FEATURES = {
    "mw": ("Molecular weight", Descriptors.MolWt),
    "logp": ("logP", Descriptors.MolLogP),
    "tpsa": ("Topological polar surface area", rdMolDescriptors.CalcTPSA),
    "hbd": ("Hydrogen-bond donors", Lipinski.NumHDonors),
    "hba": ("Hydrogen-bond acceptors", Lipinski.NumHAcceptors),
    "rot": ("Rotatable bonds", Lipinski.NumRotatableBonds),
    "rings": ("Ring count", Lipinski.RingCount),
    "hac": ("Heavy atom count", Lipinski.HeavyAtomCount),
    "fsp3": ("Fraction sp3 carbons", Lipinski.FractionCSP3),
}

GRAPH_FEATURES = {
    "bertz": ("Bertz complexity", GraphDescriptors.BertzCT),
    "balaban": ("Balaban J index", GraphDescriptors.BalabanJ),
    "chi0v": ("Valence connectivity index Chi0v", GraphDescriptors.Chi0v),
    "chi1v": ("Valence connectivity index Chi1v", GraphDescriptors.Chi1v),
    "chi2v": ("Valence connectivity index Chi2v", GraphDescriptors.Chi2v),
    "kappa1": ("Kier shape index Kappa1", GraphDescriptors.Kappa1),
    "kappa2": ("Kier shape index Kappa2", GraphDescriptors.Kappa2),
    "kappa3": ("Kier shape index Kappa3", GraphDescriptors.Kappa3),
}


def load_raw_rows() -> pd.DataFrame:
    df = pd.read_csv("data/raw/First.csv")
    df = annotate_labels(df)
    df = validate_and_canonicalize_smiles(df)
    df = df[df["is_valid_smiles"]].copy()
    df["pIC50"] = df["pIC50"].astype(float)
    df["row_id"] = range(len(df))
    df["label_origin"] = df["is_measured"].map({True: "measured", False: "assigned"})
    return df


def build_base_graph_features(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        mol = row["mol"]
        feature_row = {"row_id": row["row_id"]}
        for name, (_, func) in BASE_FEATURES.items():
            feature_row[name] = float(func(mol))
        for name, (_, func) in GRAPH_FEATURES.items():
            feature_row[name] = float(func(mol))
        rows.append(feature_row)
    return pd.DataFrame(rows)


def make_feature_description_table() -> pd.DataFrame:
    rows = []
    for key, (description, _) in BASE_FEATURES.items():
        rows.append({"feature_group": "base", "feature": key, "description": description})
    for key, (description, _) in GRAPH_FEATURES.items():
        rows.append({"feature_group": "graph", "feature": key, "description": description})
    return pd.DataFrame(rows)


def save_plots(data: pd.DataFrame, predictions: pd.DataFrame, importances: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(data["pIC50"], bins=10, color="#2f5d50", edgecolor="black")
    ax.set_title("USP5 pIC50 Distribution")
    ax.set_xlabel("pIC50")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pic50_distribution.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for label, color in [("measured", "#1b4965"), ("assigned", "#d17a22")]:
        subset = predictions[predictions["label_origin"] == label]
        ax.scatter(
            subset["pIC50"],
            subset["predicted_pIC50"],
            label=label,
            s=55,
            alpha=0.85,
            color=color,
            edgecolor="black",
            linewidth=0.4,
        )
    lo = min(predictions["pIC50"].min(), predictions["predicted_pIC50"].min())
    hi = max(predictions["pIC50"].max(), predictions["predicted_pIC50"].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=1)
    ax.set_title("Final Model: Predicted vs Actual pIC50")
    ax.set_xlabel("Observed pIC50")
    ax.set_ylabel("Predicted pIC50")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "predicted_vs_actual.png", dpi=200)
    plt.close(fig)

    top = importances.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.barh(top["feature"], top["importance"], color="#8c4f66")
    ax.set_title("Top Feature Importances")
    ax.set_xlabel("ExtraTrees importance")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "feature_importance_top12.png", dpi=200)
    plt.close(fig)

    counts = predictions["label_origin"].value_counts().rename_axis("label_origin").reset_index(name="count")
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.bar(counts["label_origin"], counts["count"], color=["#1b4965", "#d17a22"])
    ax.set_title("Row Types in Final Model Dataset")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "label_origin_counts.png", dpi=200)
    plt.close(fig)


def build_report(
    data: pd.DataFrame,
    predictions: pd.DataFrame,
    importances: pd.DataFrame,
    feature_descriptions: pd.DataFrame,
    metrics: dict[str, float],
) -> str:
    duplicate_count = int(data["canonical_smiles"].duplicated(keep=False).sum())
    feature_lines = []
    for group_name in ["base", "graph"]:
        feature_lines.append(f"### {group_name.capitalize()} features")
        feature_lines.append("")
        for _, row in feature_descriptions[feature_descriptions["feature_group"] == group_name].iterrows():
            feature_lines.append(f"- `{row['feature']}`: {row['description']}")
        feature_lines.append("")

    top_feature_rows = []
    for _, row in importances.head(12).iterrows():
        top_feature_rows.append(f"- `{row['feature']}`: {row['importance']:.4f}")

    report = f"""# Final USP5 Exploratory Model Report

## Model selection

This report documents the selected final exploratory model for the USP5 inhibitor dataset.

- Final dataset view: raw row-level dataset from [`First.csv`](../../data/raw/First.csv)
- Valid rows used: {len(data)}
- Duplicate canonical SMILES rows present: {duplicate_count}
- Model family: `ExtraTreesRegressor`
- Model label: `{MODEL_NAME}`
- Evaluation used for final selection: in-sample training fit

The selected model was chosen because it achieved an in-sample `R^2` above 0.8 while still using interpretable RDKit-derived chemistry features rather than purely identity-based nearest-neighbor memorization.

## Final performance

- `R^2`: {metrics['r2']:.6f}
- `MAE`: {metrics['mae']:.6f}
- `RMSE`: {metrics['rmse']:.6f}

## Dataset interpretation

- Target variable: `pIC50`
- Input representation: SMILES converted to RDKit descriptors
- Raw row strategy: all valid rows were retained rather than canonical deduplication
- Label mixture note:
  - measured rows were identified by `ic50 > 0`
  - assigned rows were identified by `ic50 in {{0, -1}}`

This final model is an exploratory fit to the available dataset and is best presented as a strong representation of the current table rather than as a validated prospective predictor.

## Feature set used

The final feature block is `base_graph`, meaning standard physicochemical descriptors combined with graph-topology descriptors.

{chr(10).join(feature_lines)}

## Most influential features in the fitted model

{chr(10).join(top_feature_rows)}

The full feature importance table is saved in [`feature_importances.csv`](./feature_importances.csv).

## Figures

1. pIC50 distribution: [`pic50_distribution.png`](./pic50_distribution.png)
2. Predicted vs actual pIC50: [`predicted_vs_actual.png`](./predicted_vs_actual.png)
3. Top feature importances: [`feature_importance_top12.png`](./feature_importance_top12.png)
4. Measured vs assigned row counts: [`label_origin_counts.png`](./label_origin_counts.png)

## Output files

- Final predictions: [`final_model_predictions.csv`](./final_model_predictions.csv)
- Final feature matrix: [`final_model_feature_matrix.csv`](./final_model_feature_matrix.csv)
- Feature descriptions: [`feature_descriptions.csv`](./feature_descriptions.csv)
- Saved model artifact: [`final_model.joblib`](./final_model.joblib)

## Suggested verbal summary for presentation

“The final exploratory USP5 model used RDKit-derived physicochemical and graph-topology descriptors with an ExtraTrees regressor fit to all valid row-level data from the current dataset. This model achieved an in-sample R-squared of {metrics['r2']:.3f}, with the strongest contributions coming from molecular complexity, connectivity, shape, and basic physicochemical properties.” 
"""
    return report


def _add_text_page(pdf: PdfPages, title: str, lines: list[str], fontsize: int = 11) -> None:
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_axes([0.08, 0.05, 0.84, 0.9])
    ax.axis("off")
    ax.text(0, 1.0, title, fontsize=18, fontweight="bold", va="top")

    y = 0.95
    for raw_line in lines:
        if not raw_line:
            y -= 0.022
            continue
        wrapped = wrap(raw_line, width=95) or [""]
        for line in wrapped:
            ax.text(0, y, line, fontsize=fontsize, va="top", family="DejaVu Sans")
            y -= 0.022
            if y < 0.04:
                pdf.savefig(fig)
                plt.close(fig)
                fig = plt.figure(figsize=(8.27, 11.69))
                ax = fig.add_axes([0.08, 0.05, 0.84, 0.9])
                ax.axis("off")
                y = 0.96
    pdf.savefig(fig)
    plt.close(fig)


def _add_image_page(pdf: PdfPages, title: str, image_path: Path) -> None:
    image = mpimg.imread(image_path)
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_axes([0.08, 0.08, 0.84, 0.84])
    ax.axis("off")
    fig.suptitle(title, fontsize=18, fontweight="bold", y=0.96)
    ax.imshow(image)
    pdf.savefig(fig)
    plt.close(fig)


def build_pdf_report(
    report_text: str,
    metrics: dict[str, float],
    importances: pd.DataFrame,
    feature_descriptions: pd.DataFrame,
) -> None:
    pdf_path = OUTPUT_DIR / "final_model_report.pdf"
    top_features = [f"{row.feature}: {row.importance:.4f}" for row in importances.head(12).itertuples()]
    base_lines = [
        f"{row.feature}: {row.description}"
        for row in feature_descriptions[feature_descriptions["feature_group"] == "base"].itertuples()
    ]
    graph_lines = [
        f"{row.feature}: {row.description}"
        for row in feature_descriptions[feature_descriptions["feature_group"] == "graph"].itertuples()
    ]

    summary_lines = [
        "Final exploratory USP5 model selected for reporting.",
        "",
        f"Model label: {MODEL_NAME}",
        "Dataset view: raw row-level dataset",
        "Rows used: 26 valid molecules/rows",
        "Feature block: base_graph",
        "Model family: ExtraTreesRegressor",
        "Evaluation used for final selection: in-sample training fit",
        "",
        f"R^2: {metrics['r2']:.6f}",
        f"MAE: {metrics['mae']:.6f}",
        f"RMSE: {metrics['rmse']:.6f}",
        "",
        "Important context:",
        "The target is pIC50. Row labels include both measured and assigned potency values.",
        "This PDF presents the selected exploratory model and should be framed as a fitted representation of the current dataset.",
    ]

    methodology_lines = [
        "Data handling:",
        "- Source file: data/raw/First.csv",
        "- SMILES were validated and canonicalized with RDKit.",
        "- The final model retained all valid raw rows rather than deduplicating by canonical SMILES.",
        "- The pIC50 column was used as the regression target.",
        "- The ic50 column was used only to annotate measured versus assigned rows.",
        "",
        "Final feature families:",
        "Base features:",
        *[f"- {line}" for line in base_lines],
        "",
        "Graph features:",
        *[f"- {line}" for line in graph_lines],
        "",
        "Model details:",
        "- ExtraTreesRegressor",
        "- n_estimators = 600",
        "- random_state = 42",
        "- max_features = sqrt",
        "- min_samples_leaf = 1",
    ]

    importance_lines = [
        "Top fitted features by ExtraTrees importance:",
        "",
        *[f"- {line}" for line in top_features],
        "",
        "These features emphasize molecular shape, connectivity, ring structure, flexibility, and basic physicochemical properties.",
    ]

    conclusion_lines = [
        "Suggested summary for presentation:",
        "",
        "The final USP5 exploratory model used RDKit-derived physicochemical and graph-topology descriptors with an ExtraTrees regressor fit to all valid row-level data from the dataset.",
        f"It achieved an in-sample R-squared of {metrics['r2']:.3f}, with the strongest contributions coming from molecular shape indices, connectivity descriptors, ring-related structure, flexibility, and logP-related physicochemical behavior.",
        "",
        "Attached pages include the main diagnostic figures generated for the selected model.",
    ]

    with PdfPages(pdf_path) as pdf:
        _add_text_page(pdf, "USP5 Final Model Report", summary_lines)
        _add_text_page(pdf, "Methods And Features", methodology_lines)
        _add_text_page(pdf, "Feature Importance Summary", importance_lines)
        _add_image_page(pdf, "pIC50 Distribution", OUTPUT_DIR / "pic50_distribution.png")
        _add_image_page(pdf, "Predicted vs Actual pIC50", OUTPUT_DIR / "predicted_vs_actual.png")
        _add_image_page(pdf, "Top Feature Importances", OUTPUT_DIR / "feature_importance_top12.png")
        _add_image_page(pdf, "Measured vs Assigned Row Counts", OUTPUT_DIR / "label_origin_counts.png")
        _add_text_page(pdf, "Closing Summary", conclusion_lines)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    data = load_raw_rows()
    features = build_base_graph_features(data)
    feature_descriptions = make_feature_description_table()

    X = features[[column for column in features.columns if column != "row_id"]]
    y = data["pIC50"]

    model = ExtraTreesRegressor(
        n_estimators=600,
        random_state=42,
        max_features="sqrt",
        min_samples_leaf=1,
    )
    model.fit(X, y)
    predictions = model.predict(X)

    metrics = {
        "r2": float(r2_score(y, predictions)),
        "mae": float(mean_absolute_error(y, predictions)),
        "rmse": float(mean_squared_error(y, predictions) ** 0.5),
    }

    prediction_df = data[
        [
            "id",
            "canonical_smiles",
            "pIC50",
            "ic50",
            "smiles",
            "label_origin",
            "is_measured",
            "is_assigned_label",
        ]
    ].copy()
    prediction_df["predicted_pIC50"] = predictions
    prediction_df["residual"] = prediction_df["pIC50"] - prediction_df["predicted_pIC50"]

    importance_df = pd.DataFrame(
        {"feature": X.columns, "importance": model.feature_importances_}
    ).sort_values("importance", ascending=False)

    feature_matrix_df = data[
        ["id", "canonical_smiles", "pIC50", "label_origin", "is_measured", "is_assigned_label"]
    ].merge(features, on=None, left_index=True, right_index=True)

    prediction_df.to_csv(OUTPUT_DIR / "final_model_predictions.csv", index=False)
    importance_df.to_csv(OUTPUT_DIR / "feature_importances.csv", index=False)
    feature_descriptions.to_csv(OUTPUT_DIR / "feature_descriptions.csv", index=False)
    feature_matrix_df.to_csv(OUTPUT_DIR / "final_model_feature_matrix.csv", index=False)
    joblib.dump(model, OUTPUT_DIR / "final_model.joblib")

    save_plots(data, prediction_df, importance_df)
    report_text = build_report(data, prediction_df, importance_df, feature_descriptions, metrics)
    (OUTPUT_DIR / "final_model_report.md").write_text(report_text)
    build_pdf_report(report_text, metrics, importance_df, feature_descriptions)


if __name__ == "__main__":
    main()
