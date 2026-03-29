from __future__ import annotations

import csv
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
MODELING_DATASET = OUTPUT_DIR / "modeling_dataset.csv"
RDLogger.DisableLog("rdApp.*")

SIM_MIN = 0.55
SIM_MAX = 0.95
MAX_DEPTH = 3
MAX_NEW_PER_PARENT = 350
MAX_QUEUE_PER_PARENT = 220


@dataclass(frozen=True)
class Parent:
    representative_id: str
    smiles: str
    mol: Chem.Mol


FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def canonicalize(mol: Chem.Mol | None) -> str | None:
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        return None
    smiles = Chem.MolToSmiles(mol, canonical=True)
    if "." in smiles:
        return None
    return smiles


def reaction_products(mol: Chem.Mol, smarts: str) -> list[str]:
    rxn = AllChem.ReactionFromSmarts(smarts)
    out: list[str] = []
    for products in rxn.RunReactants((mol,)):
        cs = canonicalize(products[0])
        if cs:
            out.append(cs)
    return out


def reaction_products_2(mol_a: Chem.Mol, mol_b: Chem.Mol, smarts: str) -> list[str]:
    rxn = AllChem.ReactionFromSmarts(smarts)
    out: list[str] = []
    for products in rxn.RunReactants((mol_a, mol_b)):
        cs = canonicalize(products[0])
        if cs:
            out.append(cs)
    return out


def load_positive_parents() -> list[Parent]:
    rows = list(csv.DictReader(MODELING_DATASET.open()))
    parents = []
    for row in rows:
        if row["has_measured_row"] != "True" and row["has_active_no_ic50_row"] != "True":
            continue
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            continue
        parents.append(Parent(row["representative_id"], row["canonical_smiles"], mol))
    return parents


def max_similarity_to_originals(smiles: str, original_fps: list[tuple[str, str, object]]) -> tuple[float, str, str]:
    mol = Chem.MolFromSmiles(smiles)
    fp = FP_GEN.GetFingerprint(mol)
    best_sim = -1.0
    best_id = ""
    best_smiles = ""
    for oid, osmi, ofp in original_fps:
        sim = DataStructs.TanimotoSimilarity(fp, ofp)
        if sim > best_sim:
            best_sim = sim
            best_id = oid
            best_smiles = osmi
    return best_sim, best_id, best_smiles


def local_transforms(mol: Chem.Mol) -> dict[str, list[str]]:
    products: dict[str, list[str]] = defaultdict(list)

    unary_transforms = {
        "mmp_aryl_F_to_Cl": "[c:1][F:2]>>[c:1]Cl",
        "mmp_aryl_Cl_to_F": "[c:1][Cl:2]>>[c:1]F",
        "mmp_aryl_Br_to_Cl": "[c:1][Br:2]>>[c:1]Cl",
        "mmp_aryl_Cl_to_Br": "[c:1][Cl:2]>>[c:1]Br",
        "mmp_aryl_OCH3_to_OH": "[c:1]OC>>[c:1]O",
        "mmp_aryl_OH_to_OCH3": "[c:1][OH:2]>>[c:1]OC",
        "heteroatom_CH_to_N": "[cH:1]>>[n:1]",
        "bioisostere_acid_to_amide": "[C:1](=O)[OH:2]>>[C:1](=O)N",
        "bioisostere_acid_to_thioacid": "[C:1](=O)[OH:2]>>[C:1](=S)[OH:2]",
        "bioisostere_amide_to_thioamide": "[C:1](=O)[N:2]>>[C:1](=S)[N:2]",
        "linker_carbonylCH2N_to_carbamate": "[C:1](=O)[CH2:2][N:3]>>[C:1](=O)O[CH2:2][N:3]",
        "linker_benzylNH_to_oxybenzylNH": "[c:1][CH2:2][NH:3][C:4](=O)>>[c:1]O[CH2:2][NH:3][C:4](=O)",
        "ring_pyrrolidine_to_piperidine": "[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1>>[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1",
        "ring_piperidine_to_pyrrolidine": "[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1>>[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1",
        "ring_piperidine_to_morpholine": "[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1>>[N:1]1[CH2:2][CH2:3]O[CH2:5][CH2:6]1",
    }
    for label, smarts in unary_transforms.items():
        products[label].extend(reaction_products(mol, smarts))

    amines = [
        ("N", "rxn_acid_to_amide_ammonia"),
        ("NC", "rxn_acid_to_amide_methylamine"),
        ("NCC", "rxn_acid_to_amide_ethylamine"),
        ("NCCO", "rxn_acid_to_amide_ethanolamine"),
        ("N1CCOCC1", "rxn_acid_to_amide_morpholine"),
        ("N1CCCCC1", "rxn_acid_to_amide_piperidine"),
        ("N1CCCC1", "rxn_acid_to_amide_pyrrolidine"),
    ]
    acid_to_amide = "[C:1](=O)[OH:2].[N:3]>>[C:1](=O)[N:3]"
    for amine_smiles, label in amines:
        products[label].extend(reaction_products_2(mol, Chem.MolFromSmiles(amine_smiles), acid_to_amide))

    alkyl_halides = [
        ("CBr", "rxn_phenol_o_methyl"),
        ("CCBr", "rxn_phenol_o_ethyl"),
        ("CC(C)Br", "rxn_phenol_o_isopropyl"),
        ("CCCOBr", "rxn_phenol_o_hydroxypropyl"),
    ]
    phenol_alkylation = "[c:1][OH:2].[C:3][Br,Cl,I:4]>>[c:1]O[C:3]"
    for halide_smiles, label in alkyl_halides:
        products[label].extend(reaction_products_2(mol, Chem.MolFromSmiles(halide_smiles), phenol_alkylation))

    return products


def enumerate_similarity_window(parents: list[Parent]) -> pd.DataFrame:
    original_fps = [(parent.representative_id, parent.smiles, FP_GEN.GetFingerprint(parent.mol)) for parent in parents]
    records: list[dict[str, object]] = []
    global_seen: set[str] = set(parent.smiles for parent in parents)

    for parent in parents:
        parent_records = 0
        queue: deque[tuple[str, int]] = deque([(parent.smiles, 0)])
        local_seen: set[str] = {parent.smiles}
        while queue and parent_records < MAX_NEW_PER_PARENT:
            current_smiles, depth = queue.popleft()
            if depth >= MAX_DEPTH:
                continue
            mol = Chem.MolFromSmiles(current_smiles)
            if mol is None:
                continue
            transform_map = local_transforms(mol)
            for label, products in transform_map.items():
                for product in products:
                    if product in local_seen:
                        continue
                    local_seen.add(product)
                    max_sim, nearest_id, nearest_smiles = max_similarity_to_originals(product, original_fps)
                    if not (SIM_MIN <= max_sim <= SIM_MAX):
                        continue
                    if product in global_seen:
                        continue
                    global_seen.add(product)
                    parent_records += 1
                    records.append(
                        {
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "generation_depth": depth + 1,
                            "transform": label,
                            "max_tanimoto_to_any_original": round(max_sim, 4),
                            "nearest_original_id": nearest_id,
                            "nearest_original_smiles": nearest_smiles,
                        }
                    )
                    if depth + 1 < MAX_DEPTH and len(queue) < MAX_QUEUE_PER_PARENT:
                        queue.append((product, depth + 1))
                    if parent_records >= MAX_NEW_PER_PARENT:
                        break
                if parent_records >= MAX_NEW_PER_PARENT:
                    break
    return pd.DataFrame(records)


def save_outputs(df: pd.DataFrame) -> None:
    library_path = OUTPUT_DIR / "enumeration_similarity_window_0p55_0p95.csv"
    counts_path = OUTPUT_DIR / "enumeration_similarity_window_counts.csv"
    summary_path = OUTPUT_DIR / "enumeration_similarity_window_summary.md"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = df.sort_values(["parent_id", "generation_depth", "transform", "product_smiles"]).reset_index(drop=True)
    df.to_csv(library_path, index=False)

    counts = df.groupby("transform")["product_smiles"].nunique().rename("unique_products").sort_values(ascending=False).reset_index()
    counts.to_csv(counts_path, index=False)

    lines = [
        "# Similarity-Window Enumeration Summary",
        "",
        "Target window: Morgan Tanimoto 0.55 to 0.95 versus the 13 original positive parents.",
        "",
        f"- Unique products: {df['product_smiles'].nunique()}",
        f"- Mean similarity to nearest original: {df['max_tanimoto_to_any_original'].mean():.4f}",
        f"- Median similarity to nearest original: {df['max_tanimoto_to_any_original'].median():.4f}",
        f"- Min / max similarity: {df['max_tanimoto_to_any_original'].min():.4f} / {df['max_tanimoto_to_any_original'].max():.4f}",
        "",
        "## By parent",
        "",
    ]
    by_parent = df.groupby("parent_id")["product_smiles"].nunique().sort_values(ascending=False)
    for parent_id, count in by_parent.items():
        lines.append(f"- `{parent_id}`: {count}")
    lines.extend(["", "## By transform", ""])
    for row in counts.itertuples(index=False):
        lines.append(f"- `{row.transform}`: {row.unique_products}")
    summary_path.write_text("\n".join(lines))


def main() -> None:
    parents = load_positive_parents()
    df = enumerate_similarity_window(parents)
    save_outputs(df)
    print(f"unique_products={df['product_smiles'].nunique()}")
    print(f"mean_similarity={df['max_tanimoto_to_any_original'].mean():.4f}")
    print(f"median_similarity={df['max_tanimoto_to_any_original'].median():.4f}")
    print(df.groupby('transform')['product_smiles'].nunique().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
