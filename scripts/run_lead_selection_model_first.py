from __future__ import annotations

import csv
import math
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from admet_ai import ADMETModel
from rdkit import Chem, DataStructs, RDConfig
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
from rdkit.Chem import rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
LIBRARY_PATH = OUTPUT_DIR / "enumeration_similarity_window_0p55_0p95.csv"
MODEL_PATH = OUTPUT_DIR / "final_model" / "final_model.joblib"
STRUCTURE_PATH = ROOT / "data" / "structures_7ms7.pdb"

COUNTS_PATH = OUTPUT_DIR / "lead_selection_model_first_counts.csv"
POOL_PATH = OUTPUT_DIR / "lead_selection_model_first_pool.csv"
LEADS_PATH = OUTPUT_DIR / "final_leads_model_first.csv"
SUMMARY_PATH = OUTPUT_DIR / "lead_selection_model_first_summary.md"

POTENCY_CUTOFF = 4.60
TPSA_RANGE = (45.0, 140.0)
LABUTE_ASA_RANGE = (140.0, 235.0)
ROTATABLE_MAX = 10
SHAPE_MIN = 0.55
PHARM_MIN = 0.45
CLASH_MAX = 2

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

MODEL_FEATURES = {
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


def build_alert_catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def lipinski_violations(mw: float, logp: float, hbd: float, hba: float) -> int:
    return int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)


def normalize_higher(series: pd.Series) -> pd.Series:
    lo = float(series.min())
    hi = float(series.max())
    if math.isclose(lo, hi):
        return pd.Series(np.ones(len(series)), index=series.index)
    return (series - lo) / (hi - lo)


def normalize_lower(series: pd.Series) -> pd.Series:
    return 1.0 - normalize_higher(series)


def closeness(series: pd.Series, target: float, tolerance: float) -> pd.Series:
    values = 1.0 - (series - target).abs() / tolerance
    return values.clip(lower=0.0, upper=1.0)


def load_library() -> pd.DataFrame:
    alert_catalog = build_alert_catalog()
    model = joblib.load(MODEL_PATH)

    rows = []
    for row in csv.DictReader(LIBRARY_PATH.open()):
        mol = Chem.MolFromSmiles(row["product_smiles"])
        if mol is None:
            continue

        record = {
            **row,
            "mol": mol,
            "fp": FP_GEN.GetFingerprint(mol),
            "mw": float(Descriptors.MolWt(mol)),
            "logp": float(Descriptors.MolLogP(mol)),
            "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
            "labute_asa": float(MolSurf.LabuteASA(mol)),
            "hbd": float(Lipinski.NumHDonors(mol)),
            "hba": float(Lipinski.NumHAcceptors(mol)),
            "rot": float(Lipinski.NumRotatableBonds(mol)),
            "rings": float(Lipinski.RingCount(mol)),
            "formal_charge": float(Chem.GetFormalCharge(mol)),
            "qed": float(QED.qed(mol)),
            "pains_brenk_alerts": len(alert_catalog.GetMatches(mol)),
        }
        for name, func in MODEL_FEATURES.items():
            if name in record:
                continue
            record[name] = float(func(mol))
        rows.append(record)

    df = pd.DataFrame(rows)
    df["pred_pIC50"] = model.predict(df[model.feature_names_in_])
    df["pred_ic50_uM"] = 10 ** (6 - df["pred_pIC50"])
    df["lipinski_violations"] = df.apply(
        lambda r: lipinski_violations(r["mw"], r["logp"], r["hbd"], r["hba"]), axis=1
    )
    return df


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


def extract_reference_ligand(
    pdb_path: Path,
    template_smiles: str,
    chain: str = "A",
    resname: str = "ZQ1",
    resseq: int = 302,
) -> Chem.Mol:
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
        if not line.startswith("CONECT"):
            continue
        ints = [int(token) for token in line.split()[1:]]
        if ints and ints[0] in atom_serials:
            conect_lines.append(line)
    pdb_block = "\n".join(ligand_lines + conect_lines) + "\nEND\n"
    pdb_mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if pdb_mol is None:
        raise RuntimeError("Could not parse reference ligand from 7MS7 PDB.")
    template = Chem.MolFromSmiles(template_smiles)
    return AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)


def extract_pocket_coordinates(
    pdb_path: Path,
    ligand_coords: np.ndarray,
    chain: str = "A",
    cutoff: float = 6.0,
) -> np.ndarray:
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
        best_idx = None
        best_dist = None
        cutoff = FEATURE_MATCH_CUTOFF[family]
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
    count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz = np.array([pos.x, pos.y, pos.z], dtype=float)
        if np.min(np.linalg.norm(pocket_coords - xyz, axis=1)) < threshold:
            count += 1
    return count


def score_binding_pose(
    mol: Chem.Mol,
    reference: Chem.Mol,
    reference_features: list[tuple[str, np.ndarray, float]],
    pocket_coords: np.ndarray,
    num_confs: int = 12,
) -> dict[str, float]:
    work_mol = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.pruneRmsThresh = 0.35
    params.randomSeed = 0xC0D3
    params.numThreads = 0
    conf_ids = list(AllChem.EmbedMultipleConfs(work_mol, numConfs=num_confs, params=params))
    if not conf_ids:
        return {
            "shape_tanimoto": np.nan,
            "pharmacophore_score": np.nan,
            "clash_count": np.nan,
            "binding_score": np.nan,
            "best_conf_id": np.nan,
        }

    try:
        AllChem.MMFFOptimizeMoleculeConfs(work_mol, numThreads=0)
    except Exception:
        pass

    best: dict[str, float] | None = None
    ref_crippen = rdMolDescriptors._CalcCrippenContribs(reference)
    probe_crippen = rdMolDescriptors._CalcCrippenContribs(work_mol)
    for conf_id in conf_ids:
        try:
            o3a = rdMolAlign.GetCrippenO3A(work_mol, reference, probe_crippen, ref_crippen, conf_id, 0)
            o3a.Align()
        except Exception:
            continue

        shape = 1.0 - rdShapeHelpers.ShapeTanimotoDist(work_mol, reference, confId1=conf_id, confId2=0)
        probe_features = get_feature_points(work_mol, conf_id)
        pharm = feature_overlap_score(reference_features, probe_features)
        clashes = clash_count(work_mol, conf_id, pocket_coords)
        binding_score = 0.65 * shape + 0.35 * pharm - 0.10 * clashes
        result = {
            "shape_tanimoto": float(shape),
            "pharmacophore_score": float(pharm),
            "clash_count": int(clashes),
            "binding_score": float(binding_score),
            "best_conf_id": int(conf_id),
        }
        if best is None or result["binding_score"] > best["binding_score"]:
            best = result

    if best is None:
        return {
            "shape_tanimoto": np.nan,
            "pharmacophore_score": np.nan,
            "clash_count": np.nan,
            "binding_score": np.nan,
            "best_conf_id": np.nan,
        }
    return best


def run_pipeline(df: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    stages: list[tuple[str, int]] = [("start_focused_library", int(df["product_smiles"].nunique()))]

    stage1 = df[df["pred_pIC50"] >= POTENCY_CUTOFF].copy()
    stages.append((f"stage1_predicted_pIC50_ge_{POTENCY_CUTOFF:.2f}", int(stage1["product_smiles"].nunique())))

    stage2 = stage1[stage1["pains_brenk_alerts"] == 0].copy()
    stages.append(("stage2_pains_brenk_free", int(stage2["product_smiles"].nunique())))

    stage3 = stage2[
        (stage2["lipinski_violations"] <= 1)
        & stage2["tpsa"].between(*TPSA_RANGE)
        & stage2["labute_asa"].between(*LABUTE_ASA_RANGE)
        & (stage2["rot"] <= ROTATABLE_MAX)
        & (stage2["formal_charge"].abs() <= 1)
    ].copy()
    stages.append(("stage3_lipinski_psa_surface_area_flexibility", int(stage3["product_smiles"].nunique())))

    stage4 = run_admet(stage3)
    stage4 = stage4[
        (stage4["AMES"] <= 0.20)
        & (stage4["hERG"] <= 0.40)
        & (stage4["ClinTox"] <= 0.50)
        & (stage4["HIA_Hou"] >= 0.60)
        & (stage4["Bioavailability_Ma"] >= 0.50)
    ].copy()
    stages.append(("stage4_admet_ai_core_gates", int(stage4["product_smiles"].nunique())))

    reference_smiles = "O=C(O)CNC(=O)c1ccc(S(=O)(=O)N2CCC(c3ccc(Cl)cc3)CC2)cn1"
    reference = extract_reference_ligand(STRUCTURE_PATH, reference_smiles)
    ref_conf = reference.GetConformer()
    ligand_coords = np.array(
        [
            [ref_conf.GetAtomPosition(i).x, ref_conf.GetAtomPosition(i).y, ref_conf.GetAtomPosition(i).z]
            for i in range(reference.GetNumAtoms())
        ],
        dtype=float,
    )
    pocket_coords = extract_pocket_coordinates(STRUCTURE_PATH, ligand_coords)
    reference_features = get_feature_points(reference, 0)

    binding_rows = []
    for row in stage4.itertuples(index=False):
        binding_rows.append(
            {
                "product_smiles": row.product_smiles,
                **score_binding_pose(row.mol, reference, reference_features, pocket_coords),
            }
        )
    binding_df = pd.DataFrame(binding_rows)
    stage5 = stage4.merge(binding_df, on="product_smiles", how="left")
    stage5 = stage5[
        (stage5["shape_tanimoto"] >= SHAPE_MIN)
        & (stage5["pharmacophore_score"] >= PHARM_MIN)
        & (stage5["clash_count"] <= CLASH_MAX)
    ].copy()
    stages.append(("stage5_template_docking_and_3d_pharmacophore", int(stage5["product_smiles"].nunique())))

    return stage5, stages


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    ranked = df.copy()
    ranked["potency_score"] = normalize_higher(ranked["pred_pIC50"])
    ranked["admet_score"] = (
        normalize_lower(ranked["AMES"])
        + normalize_lower(ranked["hERG"])
        + normalize_lower(ranked["ClinTox"])
        + normalize_lower(ranked["DILI"])
        + normalize_higher(ranked["HIA_Hou"])
        + normalize_higher(ranked["Bioavailability_Ma"])
        + normalize_lower(ranked["Pgp_Broccatelli"])
        + normalize_higher(ranked["Solubility_AqSolDB"])
    ) / 8.0
    ranked["property_score"] = (
        closeness(ranked["logp"], 2.7, 1.8)
        + closeness(ranked["tpsa"], 95.0, 40.0)
        + closeness(ranked["labute_asa"], 175.0, 45.0)
        + closeness(ranked["rot"], 6.0, 4.0)
        + normalize_higher(ranked["qed"])
    ) / 5.0
    ranked["binding_component"] = (
        0.45 * normalize_higher(ranked["binding_score"])
        + 0.35 * normalize_higher(ranked["shape_tanimoto"])
        + 0.20 * normalize_higher(ranked["pharmacophore_score"])
    )
    ranked["composite_score"] = (
        0.35 * ranked["potency_score"]
        + 0.25 * ranked["admet_score"]
        + 0.15 * ranked["property_score"]
        + 0.25 * ranked["binding_component"]
    )
    ranked = ranked.sort_values(
        ["composite_score", "pred_pIC50", "binding_score", "admet_score"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    return ranked


def diversity_pick(df: pd.DataFrame, n: int = 10, sim_ceiling: float = 0.78, per_parent_cap: int | None = None) -> pd.DataFrame:
    selected = []
    selected_fps = []
    parent_counts: dict[str, int] = {}

    for row in df.itertuples(index=False):
        if per_parent_cap is not None and parent_counts.get(row.parent_id, 0) >= per_parent_cap:
            continue
        max_sim = max(
            (DataStructs.TanimotoSimilarity(row.fp, fp) for fp in selected_fps),
            default=0.0,
        )
        if len(selected) < 3 or max_sim < sim_ceiling:
            selected.append(row)
            selected_fps.append(row.fp)
            parent_counts[row.parent_id] = parent_counts.get(row.parent_id, 0) + 1
        if len(selected) >= n:
            break

    if len(selected) < n:
        selected_smiles = {row.product_smiles for row in selected}
        for row in df.itertuples(index=False):
            if row.product_smiles in selected_smiles:
                continue
            if per_parent_cap is not None and parent_counts.get(row.parent_id, 0) >= per_parent_cap:
                continue
            selected.append(row)
            parent_counts[row.parent_id] = parent_counts.get(row.parent_id, 0) + 1
            if len(selected) >= n:
                break

    return pd.DataFrame(selected)


def save_outputs(pool: pd.DataFrame, stages: list[tuple[str, int]], leads: pd.DataFrame) -> None:
    pd.DataFrame(stages, columns=["stage", "remaining_unique_compounds"]).to_csv(COUNTS_PATH, index=False)
    pool.drop(columns=["mol", "fp"], errors="ignore").to_csv(POOL_PATH, index=False)
    leads.drop(columns=["mol", "fp"], errors="ignore").to_csv(LEADS_PATH, index=False)

    lines = [
        "# Model-First Lead Selection Summary",
        "",
        "This funnel starts with the trained ExtraTrees potency model, then moves through PAINS/BRENK cleanup, lead-like property filters, ADMET-AI, and a structure-aware 3D binding plausibility screen against the USP5 `7MS7` co-crystal ligand pose.",
        "",
        "## Stage counts",
        "",
    ]
    for stage, count in stages:
        lines.append(f"- `{stage}`: {count}")
    lines.extend(
        [
            "",
            "## Filter logic",
            "",
            f"- Stage 1: keep compounds with predicted `pIC50 >= {POTENCY_CUTOFF:.2f}` from the final ExtraTrees regression model.",
            "- Stage 2: remove compounds with `PAINS` or `BRENK` alerts. `NIH` alerts were not used.",
            f"- Stage 3: require `Lipinski violations <= 1`, `TPSA {TPSA_RANGE[0]:.0f}-{TPSA_RANGE[1]:.0f}`, `Labute ASA {LABUTE_ASA_RANGE[0]:.0f}-{LABUTE_ASA_RANGE[1]:.0f}`, `rotatable bonds <= {ROTATABLE_MAX}`, and `formal charge within -1 to +1`.",
            "- Stage 4: run `admet-ai` and keep compounds with low `AMES`, low `hERG`, low `ClinTox`, plus acceptable `HIA` and `oral bioavailability`.",
            "- Stage 5: align conformers into the bound USP5 `7MS7` ligand pose, require good 3D shape overlap, feature overlap, and low receptor clash count.",
            "",
            "## Final leads",
            "",
        ]
    )
    for row in leads.itertuples(index=False):
        lines.append(
            f"- `{row.product_smiles}` | score {row.composite_score:.4f} | predicted pIC50 {row.pred_pIC50:.3f} "
            f"(pred IC50 {row.pred_ic50_uM:.2f} uM) | AMES {row.AMES:.3f} | hERG {row.hERG:.3f} | "
            f"HIA {row.HIA_Hou:.3f} | Bioavailability {row.Bioavailability_Ma:.3f} | shape {row.shape_tanimoto:.3f} | "
            f"pharmacophore {row.pharmacophore_score:.3f} | clashes {int(row.clash_count)} | transform `{row.transform}`"
        )
    SUMMARY_PATH.write_text("\n".join(lines))


def main() -> None:
    df = load_library()
    final_pool, stages = run_pipeline(df)
    ranked = rank_candidates(final_pool)
    leads = diversity_pick(ranked, n=10, sim_ceiling=0.78, per_parent_cap=None)
    save_outputs(ranked, stages, leads)

    print("stage_counts")
    for stage, count in stages:
        print(stage, count)
    print("\nfinal_leads_model_first")
    cols = [
        "product_smiles",
        "pred_pIC50",
        "pred_ic50_uM",
        "AMES",
        "hERG",
        "HIA_Hou",
        "Bioavailability_Ma",
        "shape_tanimoto",
        "pharmacophore_score",
        "clash_count",
        "binding_score",
        "composite_score",
        "transform",
    ]
    print(leads[cols].to_string(index=False))


if __name__ == "__main__":
    main()
