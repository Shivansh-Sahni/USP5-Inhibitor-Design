from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from admet_ai import ADMETModel
from rdkit import Chem, DataStructs, RDConfig, RDLogger
from rdkit.Chem import (
    AllChem,
    ChemicalFeatures,
    Descriptors,
    GraphDescriptors,
    Lipinski,
    MolSurf,
    QED,
    rdMolAlign,
    rdMolDescriptors,
    rdShapeHelpers,
)
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from usp5_workflow.data import prepare_dataset


RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"

BROAD_LIBRARY_PATH = OUTPUT_DIR / "enumeration_library_10_methods.csv"
MODELING_PATH = OUTPUT_DIR / "modeling_dataset.csv"
RAW_DATA_PATH = ROOT / "data" / "raw" / "First.csv"
FINAL_MODEL_PATH = OUTPUT_DIR / "final_model" / "final_model.joblib"

STRUCTURES = [
    {
        "name": "7MS7_ZQ1",
        "pdb_path": ROOT / "data" / "structures_7ms7.pdb",
        "chain": "A",
        "resname": "ZQ1",
        "resseq": 302,
        "template_smiles": "O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1",
    },
    {
        "name": "7MS6_ZPV",
        "pdb_path": ROOT / "data" / "structures_7ms6.pdb",
        "chain": "A",
        "resname": "ZPV",
        "resseq": 301,
        "template_smiles": "O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccccc3)CC2)cc1F",
    },
    {
        "name": "7MS5_ZOG",
        "pdb_path": ROOT / "data" / "structures_7ms5.pdb",
        "chain": "A",
        "resname": "ZOG",
        "resseq": 302,
        "template_smiles": "O=C(O)CCC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(F)c(F)c3)CC2)cc1",
    },
    {
        "name": "6DXT_HHY",
        "pdb_path": ROOT / "data" / "structures_6dxt.pdb",
        "chain": "A",
        "resname": "HHY",
        "resseq": 701,
        "template_smiles": "O=C(O)CCc1nnc(-c2ccccc2)o1",
    },
]

COUNTS_PATH = OUTPUT_DIR / "lead_selection_robust_counts.csv"
STRICT_POOL_PATH = OUTPUT_DIR / "lead_selection_robust_strict_pool.csv"
RELAXED_POOL_PATH = OUTPUT_DIR / "lead_selection_robust_relaxed_pool.csv"
PRIMARY_LEADS_PATH = OUTPUT_DIR / "final_leads_robust.csv"
BACKUP_LEADS_PATH = OUTPUT_DIR / "backup_leads_robust.csv"
MODEL_METRICS_PATH = OUTPUT_DIR / "potency_consensus_model_metrics.csv"
MODEL_PREDICTIONS_PATH = OUTPUT_DIR / "potency_consensus_model_training_predictions.csv"
SUMMARY_PATH = OUTPUT_DIR / "lead_selection_robust_summary.md"

FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(
    str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
)
FEATURE_WEIGHTS = {
    "NegIonizable": 2.0,
    "Acceptor": 1.5,
    "Donor": 1.5,
    "Aromatic": 1.0,
    "Hydrophobe": 0.5,
}
FEATURE_MATCH_CUTOFF = {
    "NegIonizable": 1.8,
    "Acceptor": 1.8,
    "Donor": 1.8,
    "Aromatic": 2.0,
    "Hydrophobe": 2.2,
}
FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

DESCRIPTOR_FUNCS = {
    "mw": Descriptors.MolWt,
    "logp": Descriptors.MolLogP,
    "tpsa": rdMolDescriptors.CalcTPSA,
    "hbd": Lipinski.NumHDonors,
    "hba": Lipinski.NumHAcceptors,
    "rot": Lipinski.NumRotatableBonds,
    "rings": Lipinski.RingCount,
    "hac": Lipinski.HeavyAtomCount,
    "fsp3": Lipinski.FractionCSP3,
    "bertz": GraphDescriptors.BertzCT,
    "balaban": GraphDescriptors.BalabanJ,
    "chi0v": GraphDescriptors.Chi0v,
    "chi1v": GraphDescriptors.Chi1v,
    "chi2v": GraphDescriptors.Chi2v,
    "kappa1": GraphDescriptors.Kappa1,
    "kappa2": GraphDescriptors.Kappa2,
    "kappa3": GraphDescriptors.Kappa3,
}

PRIMARY_LIMIT = 8
BACKUP_LIMIT = 6


@dataclass
class ReferenceTemplate:
    name: str
    ligand: Chem.Mol
    pocket_coords: np.ndarray
    features: list[tuple[str, np.ndarray, float]]


def build_alert_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def build_feature_row(mol: Chem.Mol) -> dict[str, float]:
    return {name: float(func(mol)) for name, func in DESCRIPTOR_FUNCS.items()}


def lipinski_violations(row: pd.Series) -> int:
    return int(row["mw"] > 500) + int(row["logp"] > 5) + int(row["hbd"] > 5) + int(row["hba"] > 10)


def normalize_higher(series: pd.Series) -> pd.Series:
    lo = float(series.min())
    hi = float(series.max())
    if math.isclose(lo, hi):
        return pd.Series(np.ones(len(series)), index=series.index)
    return (series - lo) / (hi - lo)


def normalize_lower(series: pd.Series) -> pd.Series:
    return 1.0 - normalize_higher(series)


def closeness(series: pd.Series, target: float, tolerance: float) -> pd.Series:
    return (1.0 - (series - target).abs() / tolerance).clip(0.0, 1.0)


def murcko_scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)


def load_parent_meta() -> pd.DataFrame:
    df = pd.read_csv(MODELING_PATH)
    positives = df[(df["has_measured_row"]) | (df["has_active_no_ic50_row"])].copy()
    positives["product_smiles"] = positives["canonical_smiles"]
    positives["method"] = "original_positive"
    positives["parent_id"] = positives["representative_id"]
    positives["parent_smiles"] = positives["canonical_smiles"]
    positives["transform"] = "original_positive"
    return positives


def aggregate_library() -> pd.DataFrame:
    broad = pd.read_csv(BROAD_LIBRARY_PATH)
    positives = load_parent_meta()[
        [
            "product_smiles",
            "method",
            "parent_id",
            "parent_smiles",
            "transform",
            "representative_id",
            "target_pIC50",
            "target_origin",
        ]
    ].copy()
    positives["representative_id"] = positives["representative_id"].astype(str)
    broad["representative_id"] = broad["parent_id"].astype(str)
    broad["target_pIC50"] = np.nan
    broad["target_origin"] = ""

    combined = pd.concat(
        [
            broad[
                ["product_smiles", "method", "parent_id", "parent_smiles", "transform", "representative_id", "target_pIC50", "target_origin"]
            ],
            positives,
        ],
        ignore_index=True,
    )

    grouped_rows = []
    for product_smiles, group in combined.groupby("product_smiles", sort=False):
        grouped_rows.append(
            {
                "product_smiles": product_smiles,
                "parent_ids": ";".join(sorted(set(group["parent_id"].astype(str)))),
                "primary_parent_id": str(group["parent_id"].iloc[0]),
                "methods": ";".join(sorted(set(group["method"].astype(str)))),
                "transforms": ";".join(sorted(set(group["transform"].astype(str)))),
                "n_unique_parents": int(group["parent_id"].astype(str).nunique()),
                "n_unique_methods": int(group["method"].astype(str).nunique()),
                "is_original_positive": bool((group["method"] == "original_positive").any()),
            }
        )
    return pd.DataFrame(grouped_rows)


def prepare_candidates() -> pd.DataFrame:
    alert_catalog = build_alert_catalog()
    library = aggregate_library()
    rows = []
    for row in library.itertuples(index=False):
        mol = Chem.MolFromSmiles(row.product_smiles)
        if mol is None:
            continue
        feature_row = build_feature_row(mol)
        feature_row.update(
            {
                "product_smiles": row.product_smiles,
                "mol": mol,
                "fp": FP_GEN.GetFingerprint(mol),
                "parent_ids": row.parent_ids,
                "primary_parent_id": row.primary_parent_id,
                "methods": row.methods,
                "transforms": row.transforms,
                "n_unique_parents": row.n_unique_parents,
                "n_unique_methods": row.n_unique_methods,
                "is_original_positive": row.is_original_positive,
                "labute_asa": float(MolSurf.LabuteASA(mol)),
                "formal_charge": float(Chem.GetFormalCharge(mol)),
                "qed": float(QED.qed(mol)),
                "pains_brenk_alerts": len(alert_catalog.GetMatches(mol)),
                "scaffold": murcko_scaffold(row.product_smiles),
            }
        )
        rows.append(feature_row)
    df = pd.DataFrame(rows)
    df["lipinski_violations"] = df.apply(lipinski_violations, axis=1)
    return df


def prepare_training_set() -> tuple[pd.DataFrame, pd.DataFrame]:
    bundle = prepare_dataset(RAW_DATA_PATH)
    modeling = bundle.modeling.copy().reset_index(drop=True)
    features = []
    for row in modeling.itertuples(index=False):
        feature_row = build_feature_row(row.mol)
        feature_row["canonical_smiles"] = row.canonical_smiles
        features.append(feature_row)
    feature_df = pd.DataFrame(features)
    return modeling, feature_df


def run_admet(df: pd.DataFrame) -> pd.DataFrame:
    admet = ADMETModel(num_workers=0)
    predictions = admet.predict(df["product_smiles"].tolist())
    keep = [
        "AMES",
        "hERG",
        "ClinTox",
        "DILI",
        "HIA_Hou",
        "Bioavailability_Ma",
        "Pgp_Broccatelli",
        "Caco2_Wang",
        "Solubility_AqSolDB",
        "Lipophilicity_AstraZeneca",
    ]
    return df.join(predictions[keep], on="product_smiles")


def build_consensus_models() -> dict[str, object]:
    final_extra_trees = joblib.load(FINAL_MODEL_PATH)
    return {
        "final_extratrees_trainfit": final_extra_trees,
        "random_forest_loocv_best": RandomForestRegressor(
            n_estimators=200,
            random_state=42,
            min_samples_leaf=1,
            max_features=0.3,
            max_depth=None,
        ),
        "knn_k5_distance": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", KNeighborsRegressor(n_neighbors=5, weights="distance")),
            ]
        ),
        "elastic_net_conservative": Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", ElasticNet(alpha=0.5, l1_ratio=0.2, random_state=42, max_iter=20000)),
            ]
        ),
    }


def evaluate_models(modeling: pd.DataFrame, feature_df: pd.DataFrame, models: dict[str, object]) -> tuple[pd.DataFrame, pd.DataFrame]:
    X = feature_df[[c for c in feature_df.columns if c != "canonical_smiles"]]
    y = modeling["target_pIC50"].to_numpy()
    loo = LeaveOneOut()

    metric_rows = []
    prediction_rows = []
    for model_name, model in models.items():
        preds = np.zeros(len(y), dtype=float)
        for train_idx, test_idx in loo.split(X):
            estimator = model
            estimator.fit(X.iloc[train_idx], y[train_idx])
            preds[test_idx] = estimator.predict(X.iloc[test_idx])
        rmse = float(np.sqrt(mean_squared_error(y, preds)))
        metric_rows.append(
            {
                "model": model_name,
                "mae": float(mean_absolute_error(y, preds)),
                "rmse": rmse,
                "r2": float(r2_score(y, preds)),
            }
        )
        for idx, pred in enumerate(preds):
            prediction_rows.append(
                {
                    "model": model_name,
                    "representative_id": modeling.loc[idx, "representative_id"],
                    "canonical_smiles": modeling.loc[idx, "canonical_smiles"],
                    "observed_pIC50": y[idx],
                    "loocv_prediction": float(pred),
                    "residual": float(y[idx] - pred),
                }
            )
    return pd.DataFrame(metric_rows).sort_values(["mae", "rmse"]).reset_index(drop=True), pd.DataFrame(prediction_rows)


def apply_models(df: pd.DataFrame, training_features: pd.DataFrame, models: dict[str, object]) -> pd.DataFrame:
    train_X = training_features[[c for c in training_features.columns if c not in {"canonical_smiles", "target_pIC50"}]]
    screen_X = df[list(DESCRIPTOR_FUNCS.keys())].copy()

    for model_name, model in models.items():
        estimator = model
        estimator.fit(train_X, training_features["target_pIC50"])
        df[f"pred_{model_name}"] = estimator.predict(screen_X)

    pred_cols = [c for c in df.columns if c.startswith("pred_")]
    df["pred_pIC50_consensus"] = df[pred_cols].mean(axis=1)
    df["pred_pIC50_std"] = df[pred_cols].std(axis=1, ddof=0)
    df["pred_ic50_uM_consensus"] = 10 ** (6 - df["pred_pIC50_consensus"])

    train_desc = training_features[list(DESCRIPTOR_FUNCS.keys())].copy()
    mean = train_desc.mean()
    std = train_desc.std(ddof=0).replace(0, 1.0)
    z = (screen_X - mean) / std
    df["descriptor_ad_distance"] = np.sqrt((z**2).mean(axis=1))
    return df


def extract_reference_ligand(pdb_path: Path, chain: str, resname: str, resseq: int, template_smiles: str) -> Chem.Mol:
    lines = pdb_path.read_text().splitlines()
    ligand_lines = [
        line
        for line in lines
        if line.startswith("HETATM")
        and line[17:20].strip() == resname
        and line[21].strip() == chain
        and int(line[22:26]) == resseq
    ]
    atom_serials = {int(line[6:11]) for line in ligand_lines}
    conect_lines = []
    for line in lines:
        if line.startswith("CONECT"):
            items = [int(tok) for tok in line.split()[1:]]
            if items and items[0] in atom_serials:
                conect_lines.append(line)
    pdb_block = "\n".join(ligand_lines + conect_lines) + "\nEND\n"
    pdb_mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if pdb_mol is None:
        raise RuntimeError(f"Could not parse ligand {resname} from {pdb_path.name}.")
    template = Chem.MolFromSmiles(template_smiles)
    return AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)


def extract_pocket_coordinates(pdb_path: Path, ligand_coords: np.ndarray, chain: str, cutoff: float = 6.0) -> np.ndarray:
    coords: list[list[float]] = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21].strip() != chain:
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element.upper() == "H":
            continue
        xyz = np.array(
            [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]
        )
        if np.min(np.linalg.norm(ligand_coords - xyz, axis=1)) <= cutoff:
            coords.append(xyz.tolist())
    return np.array(coords, dtype=float)


def get_feature_points(mol: Chem.Mol, conf_id: int) -> list[tuple[str, np.ndarray, float]]:
    feats = []
    for feat in FEATURE_FACTORY.GetFeaturesForMol(mol):
        family = feat.GetFamily()
        if family not in FEATURE_WEIGHTS:
            continue
        pos = feat.GetPos(confId=conf_id)
        feats.append(
            (
                family,
                np.array([pos.x, pos.y, pos.z], dtype=float),
                FEATURE_WEIGHTS[family],
            )
        )
    return feats


def feature_overlap_score(
    reference_features: list[tuple[str, np.ndarray, float]],
    probe_features: list[tuple[str, np.ndarray, float]],
) -> float:
    total_weight = sum(weight for _, _, weight in reference_features)
    if total_weight == 0:
        return 0.0
    matched_weight = 0.0
    used_probe: set[int] = set()
    for family, ref_pos, weight in reference_features:
        cutoff = FEATURE_MATCH_CUTOFF[family]
        best_idx = None
        best_dist = None
        for idx, (probe_family, probe_pos, _) in enumerate(probe_features):
            if idx in used_probe or probe_family != family:
                continue
            dist = float(np.linalg.norm(ref_pos - probe_pos))
            if dist <= cutoff and (best_dist is None or dist < best_dist):
                best_idx = idx
                best_dist = dist
        if best_idx is not None:
            used_probe.add(best_idx)
            matched_weight += weight
    return matched_weight / total_weight


def clash_count(mol: Chem.Mol, conf_id: int, pocket_coords: np.ndarray, threshold: float = 1.6) -> int:
    conf = mol.GetConformer(conf_id)
    clashes = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz = np.array([pos.x, pos.y, pos.z], dtype=float)
        if np.min(np.linalg.norm(pocket_coords - xyz, axis=1)) < threshold:
            clashes += 1
    return clashes


def load_reference_templates() -> list[ReferenceTemplate]:
    refs = []
    for item in STRUCTURES:
        ligand = extract_reference_ligand(
            item["pdb_path"],
            item["chain"],
            item["resname"],
            item["resseq"],
            item["template_smiles"],
        )
        conf = ligand.GetConformer()
        ligand_coords = np.array(
            [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(ligand.GetNumAtoms())],
            dtype=float,
        )
        pocket_coords = extract_pocket_coordinates(item["pdb_path"], ligand_coords, item["chain"])
        refs.append(
            ReferenceTemplate(
                name=item["name"],
                ligand=ligand,
                pocket_coords=pocket_coords,
                features=get_feature_points(ligand, 0),
            )
        )
    return refs


def score_against_reference(mol: Chem.Mol, ref: ReferenceTemplate, num_confs: int = 12) -> dict[str, float]:
    work = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.pruneRmsThresh = 0.35
    params.randomSeed = 0xC0D3
    params.numThreads = 0
    conf_ids = list(AllChem.EmbedMultipleConfs(work, numConfs=num_confs, params=params))
    if not conf_ids:
        return {"shape_tanimoto": np.nan, "pharmacophore_score": np.nan, "clash_count": np.nan, "binding_score": np.nan}
    try:
        AllChem.MMFFOptimizeMoleculeConfs(work, numThreads=0)
    except Exception:
        pass

    ref_crippen = rdMolDescriptors._CalcCrippenContribs(ref.ligand)
    probe_crippen = rdMolDescriptors._CalcCrippenContribs(work)
    best = None
    for conf_id in conf_ids:
        try:
            o3a = rdMolAlign.GetCrippenO3A(work, ref.ligand, probe_crippen, ref_crippen, conf_id, 0)
            o3a.Align()
        except Exception:
            continue
        shape = 1.0 - rdShapeHelpers.ShapeTanimotoDist(work, ref.ligand, confId1=conf_id, confId2=0)
        pharm = feature_overlap_score(ref.features, get_feature_points(work, conf_id))
        clashes = clash_count(work, conf_id, ref.pocket_coords)
        binding = 0.60 * shape + 0.40 * pharm - 0.10 * clashes
        result = {
            "shape_tanimoto": float(shape),
            "pharmacophore_score": float(pharm),
            "clash_count": int(clashes),
            "binding_score": float(binding),
        }
        if best is None or result["binding_score"] > best["binding_score"]:
            best = result
    if best is None:
        return {"shape_tanimoto": np.nan, "pharmacophore_score": np.nan, "clash_count": np.nan, "binding_score": np.nan}
    return best


def run_multistructure_scoring(df: pd.DataFrame) -> pd.DataFrame:
    references = load_reference_templates()
    rows = []
    for row in df.itertuples(index=False):
        per_ref = {}
        binding_scores = []
        shape_scores = []
        pharm_scores = []
        clash_scores = []
        for ref in references:
            result = score_against_reference(row.mol, ref)
            per_ref.update(
                {
                    f"{ref.name}_shape": result["shape_tanimoto"],
                    f"{ref.name}_pharm": result["pharmacophore_score"],
                    f"{ref.name}_clash": result["clash_count"],
                    f"{ref.name}_binding": result["binding_score"],
                }
            )
            if not np.isnan(result["binding_score"]):
                binding_scores.append(result["binding_score"])
                shape_scores.append(result["shape_tanimoto"])
                pharm_scores.append(result["pharmacophore_score"])
                clash_scores.append(result["clash_count"])
        binding_scores_sorted = sorted(binding_scores, reverse=True)
        mean_top2 = float(np.mean(binding_scores_sorted[:2])) if binding_scores_sorted else np.nan
        rows.append(
            {
                "product_smiles": row.product_smiles,
                "best_binding_score": max(binding_scores) if binding_scores else np.nan,
                "mean_top2_binding_score": mean_top2,
                "best_shape_tanimoto": max(shape_scores) if shape_scores else np.nan,
                "best_pharmacophore_score": max(pharm_scores) if pharm_scores else np.nan,
                "min_clash_count": min(clash_scores) if clash_scores else np.nan,
                **per_ref,
            }
        )
    return df.merge(pd.DataFrame(rows), on="product_smiles", how="left")


def select_portfolio(
    df: pd.DataFrame,
    n: int,
    dominant_scaffold_cap: int | None = None,
    exclude_smiles: set[str] | None = None,
    parent_cap: int | None = None,
) -> pd.DataFrame:
    exclude_smiles = exclude_smiles or set()
    selected = []
    selected_fps = []
    scaffold_counts: dict[str, int] = {}
    parent_counts: dict[str, int] = {}
    if dominant_scaffold_cap is not None and not df.empty:
        dominant_scaffold = df["scaffold"].mode().iat[0]
    else:
        dominant_scaffold = None

    for row in df.itertuples(index=False):
        if row.product_smiles in exclude_smiles:
            continue
        if dominant_scaffold is not None and row.scaffold == dominant_scaffold and scaffold_counts.get(row.scaffold, 0) >= dominant_scaffold_cap:
            continue
        if parent_cap is not None and parent_counts.get(row.primary_parent_id, 0) >= parent_cap:
            continue
        max_sim = max((DataStructs.TanimotoSimilarity(row.fp, fp) for fp in selected_fps), default=0.0)
        if len(selected) < 3 or max_sim < 0.78 or row.scaffold not in scaffold_counts:
            selected.append(row)
            selected_fps.append(row.fp)
            scaffold_counts[row.scaffold] = scaffold_counts.get(row.scaffold, 0) + 1
            parent_counts[row.primary_parent_id] = parent_counts.get(row.primary_parent_id, 0) + 1
        if len(selected) >= n:
            break

    if len(selected) < n:
        current = {row.product_smiles for row in selected}
        for row in df.itertuples(index=False):
            if row.product_smiles in exclude_smiles or row.product_smiles in current:
                continue
            if dominant_scaffold is not None and row.scaffold == dominant_scaffold and scaffold_counts.get(row.scaffold, 0) >= dominant_scaffold_cap:
                continue
            if parent_cap is not None and parent_counts.get(row.primary_parent_id, 0) >= parent_cap:
                continue
            selected.append(row)
            scaffold_counts[row.scaffold] = scaffold_counts.get(row.scaffold, 0) + 1
            parent_counts[row.primary_parent_id] = parent_counts.get(row.primary_parent_id, 0) + 1
            if len(selected) >= n:
                break
    return pd.DataFrame(selected)


def main() -> None:
    candidates = prepare_candidates()
    modeling, training_features = prepare_training_set()
    training_with_target = training_features.merge(
        modeling[["canonical_smiles", "target_pIC50"]],
        on="canonical_smiles",
        how="left",
    )

    models = build_consensus_models()
    metric_df, loocv_predictions = evaluate_models(modeling, training_features, models)
    metric_df.to_csv(MODEL_METRICS_PATH, index=False)
    loocv_predictions.to_csv(MODEL_PREDICTIONS_PATH, index=False)

    screened = apply_models(candidates, training_with_target, models)
    screened["potency_support_score"] = (
        normalize_higher(screened["pred_pIC50_consensus"])
        + normalize_lower(screened["pred_pIC50_std"])
        + normalize_lower(screened["descriptor_ad_distance"])
    ) / 3.0

    stage_counts = [("start_broad_enumeration_plus_original_positives", int(screened["product_smiles"].nunique()))]

    stage1 = screened[
        (screened["pred_pIC50_consensus"] >= 4.25)
        & (screened["pred_pIC50_std"] <= 0.75)
        & (screened["descriptor_ad_distance"] <= 4.5)
    ].copy()
    stage_counts.append(("stage1_consensus_potency_and_uncertainty", int(stage1["product_smiles"].nunique())))

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

    stage4 = run_admet(stage3)
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

    stage5 = run_multistructure_scoring(stage4)
    strict = stage5[
        (stage5["best_binding_score"] >= 0.68)
        & (stage5["mean_top2_binding_score"] >= 0.56)
        & (stage5["best_shape_tanimoto"] >= 0.58)
        & (stage5["best_pharmacophore_score"] >= 0.50)
        & (stage5["min_clash_count"] <= 2)
    ].copy()
    relaxed = stage5[
        (stage5["best_binding_score"] >= 0.58)
        & (stage5["best_shape_tanimoto"] >= 0.50)
        & (stage5["best_pharmacophore_score"] >= 0.40)
        & (stage5["min_clash_count"] <= 3)
    ].copy()
    stage_counts.append(("stage5_multistructure_template_docking_strict", int(strict["product_smiles"].nunique())))
    stage_counts.append(("stage5b_multistructure_template_docking_relaxed", int(relaxed["product_smiles"].nunique())))

    for frame in [stage5, strict, relaxed]:
        frame["admet_score"] = (
            normalize_lower(frame["AMES"])
            + normalize_lower(frame["hERG"])
            + normalize_lower(frame["ClinTox"])
            + normalize_lower(frame["DILI"])
            + normalize_higher(frame["HIA_Hou"])
            + normalize_higher(frame["Bioavailability_Ma"])
            + normalize_higher(frame["Solubility_AqSolDB"])
        ) / 7.0
        frame["property_score"] = (
            closeness(frame["logp"], 2.7, 1.8)
            + closeness(frame["tpsa"], 95.0, 45.0)
            + closeness(frame["labute_asa"], 175.0, 45.0)
            + closeness(frame["rot"], 6.0, 4.0)
            + normalize_higher(frame["qed"])
        ) / 5.0
        frame["binding_component"] = (
            0.45 * normalize_higher(frame["best_binding_score"])
            + 0.25 * normalize_higher(frame["mean_top2_binding_score"])
            + 0.20 * normalize_higher(frame["best_shape_tanimoto"])
            + 0.10 * normalize_higher(frame["best_pharmacophore_score"])
        )
        frame["robust_composite_score"] = (
            0.25 * normalize_higher(frame["pred_pIC50_consensus"])
            + 0.15 * normalize_lower(frame["pred_pIC50_std"])
            + 0.10 * normalize_lower(frame["descriptor_ad_distance"])
            + 0.20 * frame["admet_score"]
            + 0.10 * frame["property_score"]
            + 0.20 * frame["binding_component"]
        )

    strict = strict.sort_values(["robust_composite_score", "best_binding_score", "pred_pIC50_consensus"], ascending=[False, False, False]).reset_index(drop=True)
    relaxed = relaxed.sort_values(["robust_composite_score", "best_binding_score", "pred_pIC50_consensus"], ascending=[False, False, False]).reset_index(drop=True)

    primary = select_portfolio(strict, PRIMARY_LIMIT, dominant_scaffold_cap=4)
    dominant_parent = primary["primary_parent_id"].mode().iat[0] if not primary.empty else None
    exploratory = stage5.copy()
    if dominant_parent is not None:
        exploratory = exploratory[exploratory["primary_parent_id"] != dominant_parent].copy()
    exploratory = exploratory.sort_values(
        ["robust_composite_score", "best_binding_score", "pred_pIC50_consensus"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    backup_source = exploratory if not exploratory.empty else relaxed
    backup = select_portfolio(backup_source, BACKUP_LIMIT, dominant_scaffold_cap=None, exclude_smiles=set(primary["product_smiles"]))

    pd.DataFrame(stage_counts, columns=["stage", "remaining_unique_compounds"]).to_csv(COUNTS_PATH, index=False)
    strict.drop(columns=["mol", "fp"], errors="ignore").to_csv(STRICT_POOL_PATH, index=False)
    relaxed.drop(columns=["mol", "fp"], errors="ignore").to_csv(RELAXED_POOL_PATH, index=False)
    primary.drop(columns=["mol", "fp"], errors="ignore").to_csv(PRIMARY_LEADS_PATH, index=False)
    backup.drop(columns=["mol", "fp"], errors="ignore").to_csv(BACKUP_LEADS_PATH, index=False)

    lines = [
        "# Robust Multi-Evidence USP5 Lead Selection",
        "",
        "This report intentionally shifts away from over-reliance on any single signal. The workflow uses a broader starting library, a consensus potency score with uncertainty, PAINS/BRENK cleanup, lead-like properties, ADMET-AI, and multi-structure 3D binding plausibility against multiple USP5 ZnF-UBD co-crystals.",
        "",
        "## Stage counts",
        "",
    ]
    for stage, count in stage_counts:
        lines.append(f"- `{stage}`: {count}")

    lines.extend(
        [
            "",
            "## Potency modeling note",
            "",
            "The potency model was treated conservatively. Leave-one-compound-out performance on the deduplicated dataset is modest, so consensus potency was used as one screen among several rather than as decisive proof of activity.",
            "",
            "## Structural references",
            "",
            "- `6DXT` small-molecule oxadiazole acid binder.",
            "- `7MS5` difluorophenyl-piperidine sulfonamide keto-acid binder.",
            "- `7MS6` fluorobenzoyl-glycine sulfonamide binder.",
            "- `7MS7` chlorophenyl-piperidine sulfonamide glycine binder.",
            "",
            "## Primary leads",
            "",
        ]
    )
    for row in primary.itertuples(index=False):
        lines.append(
            f"- `{row.product_smiles}` | score {row.robust_composite_score:.4f} | consensus pIC50 {row.pred_pIC50_consensus:.3f} | "
            f"potency std {row.pred_pIC50_std:.3f} | AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | "
            f"best binding {row.best_binding_score:.3f} | best shape {row.best_shape_tanimoto:.3f} | "
            f"best pharm {row.best_pharmacophore_score:.3f} | scaffold `{row.scaffold}` | methods `{row.methods}`"
        )
    lines.extend(["", "## Diverse backup leads", ""])
    if backup.empty:
        lines.append("- No orthogonal backup leads survived after excluding the dominant parent series.")
    else:
        for row in backup.itertuples(index=False):
            lines.append(
                f"- `{row.product_smiles}` | score {row.robust_composite_score:.4f} | consensus pIC50 {row.pred_pIC50_consensus:.3f} | "
                f"AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | best binding {row.best_binding_score:.3f} | "
                f"scaffold `{row.scaffold}` | primary parent `{row.primary_parent_id}`"
            )
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- The potency models remain low-data models and should not be overinterpreted.",
            "- The 3D stage is a multistructure template-docking surrogate, not a production docking campaign with force-field or MD refinement.",
            "- If one chemotype dominates strict survivors, that should be treated as current evidence concentration, not proof that other scaffolds are inactive.",
            "",
            "## Literature anchors",
            "",
            "- Mann et al. 2019, J Med Chem. Discovery of Small Molecule Antagonists of the USP5 Zinc Finger Ubiquitin-Binding Domain. PMID `31663737`.",
            "- Wang et al. 2024, Comput Biol Med. Structure-based virtual screening of novel USP5 inhibitors targeting the zinc finger ubiquitin-binding domain. PMID `38603896`.",
            "- RCSB structures used: `6DXT`, `7MS5`, `7MS6`, `7MS7`.",
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
        print(primary[["product_smiles", "robust_composite_score", "pred_pIC50_consensus", "pred_pIC50_std", "best_binding_score", "scaffold"]].to_string(index=False))
    print("\nbackup_leads")
    if backup.empty:
        print("none")
    else:
        print(backup[["product_smiles", "robust_composite_score", "pred_pIC50_consensus", "best_binding_score", "scaffold"]].to_string(index=False))


if __name__ == "__main__":
    main()
