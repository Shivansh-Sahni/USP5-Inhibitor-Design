from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, BRICS


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs"
MODELING_DATASET = OUTPUT_DIR / "modeling_dataset.csv"
RDLogger.DisableLog("rdApp.*")


@dataclass(frozen=True)
class Parent:
    representative_id: str
    smiles: str
    mol: Chem.Mol


def load_positive_parents() -> list[Parent]:
    rows = list(csv.DictReader(MODELING_DATASET.open()))
    parents = []
    for row in rows:
        if row["has_measured_row"] != "True" and row["has_active_no_ic50_row"] != "True":
            continue
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            continue
        parents.append(
            Parent(
                representative_id=row["representative_id"],
                smiles=row["canonical_smiles"],
                mol=mol,
            )
        )
    return parents


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


def classify_brics_fragments(parents: list[Parent]) -> tuple[list[tuple[str, Chem.Mol]], list[tuple[str, Chem.Mol]]]:
    core_frags: dict[str, Chem.Mol] = {}
    side_frags: dict[str, Chem.Mol] = {}
    for parent in parents:
        for frag_smiles in BRICS.BRICSDecompose(parent.mol):
            frag = Chem.MolFromSmiles(frag_smiles)
            if frag is None:
                continue
            dummy_count = sum(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms())
            ring_count = frag.GetRingInfo().NumRings()
            heavy_atoms = frag.GetNumHeavyAtoms()
            if dummy_count >= 2 and ring_count >= 1 and heavy_atoms >= 5:
                core_frags.setdefault(frag_smiles, frag)
            elif dummy_count == 1 and heavy_atoms >= 2:
                side_frags.setdefault(frag_smiles, frag)
    return list(core_frags.items()), list(side_frags.items())


def limited_brics_build(fragments: list[Chem.Mol], limit: int) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    try:
        generator = BRICS.BRICSBuild(fragments, uniquify=True)
        for mol in generator:
            cs = canonicalize(mol)
            if cs and cs not in seen:
                seen.add(cs)
                out.append(cs)
                if len(out) >= limit:
                    break
    except Exception:
        return out
    return out


def enumerate_r_group(parents: list[Parent], side_frags: list[tuple[str, Chem.Mol]]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    global_side_pool = [frag for _, frag in side_frags[:6]]
    for parent in parents:
        parent_frags = [Chem.MolFromSmiles(s) for s in BRICS.BRICSDecompose(parent.mol)]
        parent_frags = [frag for frag in parent_frags if frag is not None]
        cores = []
        for frag in parent_frags:
            dummy_count = sum(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms())
            if dummy_count >= 2 and frag.GetRingInfo().NumRings() >= 1:
                cores.append(frag)
        for core in cores[:2]:
            products = limited_brics_build([core, *global_side_pool], limit=25)
            for product in products:
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "r_group_enumeration",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": "fixed_core_with_swapped_side_fragments",
                        }
                    )
    return records


def enumerate_reaction_based(parents: list[Parent]) -> list[dict[str, str]]:
    amines = [
        ("N", "ammonia"),
        ("NC", "methylamine"),
        ("NCC", "ethylamine"),
        ("NCCO", "ethanolamine"),
        ("N1CCCCC1", "piperidine"),
    ]
    alkyl_halides = [
        ("CBr", "methyl_bromide"),
        ("CCBr", "ethyl_bromide"),
        ("CC(C)Br", "isopropyl_bromide"),
    ]
    acid_to_amide = "[C:1](=O)[OH:2].[N:3]>>[C:1](=O)[N:3]"
    phenol_alkylation = "[c:1][OH:2].[C:3][Br,Cl,I:4]>>[c:1]O[C:3]"
    records: list[dict[str, str]] = []
    for parent in parents:
        for amine_smiles, label in amines:
            for product in reaction_products_2(parent.mol, Chem.MolFromSmiles(amine_smiles), acid_to_amide):
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "reaction_based_enumeration",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": f"acid_to_amide_with_{label}",
                        }
                    )
        for halide_smiles, label in alkyl_halides:
            for product in reaction_products_2(parent.mol, Chem.MolFromSmiles(halide_smiles), phenol_alkylation):
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "reaction_based_enumeration",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": f"phenol_o_alkylation_with_{label}",
                        }
                    )
    return records


def enumerate_reagent_pool_combinatorics(core_frags: list[tuple[str, Chem.Mol]], side_frags: list[tuple[str, Chem.Mol]]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    build_pool = [frag for _, frag in core_frags[:3]] + [frag for _, frag in side_frags[:6]]
    for product in limited_brics_build(build_pool, limit=120):
        records.append(
            {
                "method": "reagent_pool_combinatorics",
                "parent_id": "multiple",
                "parent_smiles": "multiple",
                "product_smiles": product,
                "transform": "global_brics_pool_build",
            }
        )
    return records


def enumerate_matched_pairs(parents: list[Parent]) -> list[dict[str, str]]:
    transforms = {
        "aryl_F_to_Cl": "[c:1][F:2]>>[c:1]Cl",
        "aryl_Cl_to_F": "[c:1][Cl:2]>>[c:1]F",
        "aryl_Br_to_Cl": "[c:1][Br:2]>>[c:1]Cl",
        "aryl_F_to_CN": "[c:1][F:2]>>[c:1]C#N",
        "aryl_OCH3_to_OH": "[c:1]OC>>[c:1]O",
        "aryl_OH_to_OCH3": "[c:1][OH:2]>>[c:1]OC",
    }
    records: list[dict[str, str]] = []
    for parent in parents:
        for label, smarts in transforms.items():
            for product in reaction_products(parent.mol, smarts):
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "matched_molecular_pair_expansion",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": label,
                        }
                    )
    return records


def enumerate_bioisosteres(parents: list[Parent]) -> list[dict[str, str]]:
    transforms = {
        "acid_to_amide": "[C:1](=O)[OH:2]>>[C:1](=O)N",
        "acid_to_thioacid": "[C:1](=O)[OH:2]>>[C:1](=S)[OH:2]",
        "amide_to_thioamide": "[C:1](=O)[N:2]>>[C:1](=S)[N:2]",
    }
    records: list[dict[str, str]] = []
    for parent in parents:
        for label, smarts in transforms.items():
            for product in reaction_products(parent.mol, smarts):
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "bioisosteric_swaps",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": label,
                        }
                    )
    return records


def enumerate_linker_scans(parents: list[Parent]) -> list[dict[str, str]]:
    transforms = {
        "benzylamine_to_oxybenzylamine": "[c:1][CH2:2][NH:3][C:4](=O)>>[c:1]O[CH2:2][NH:3][C:4](=O)",
        "carbonyl_ch2_amine_to_carbamate_like": "[C:1](=O)[CH2:2][N:3]>>[C:1](=O)O[CH2:2][N:3]",
        "biaryl_ethyl_to_biaryl_oxyethyl": "[c:1][CH2:2][CH2:3][c:4]>>[c:1]O[CH2:2][CH2:3][c:4]",
    }
    records: list[dict[str, str]] = []
    for parent in parents:
        for label, smarts in transforms.items():
            for product in reaction_products(parent.mol, smarts):
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "linker_scans",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": label,
                        }
                    )
    return records


def enumerate_ring_scans(parents: list[Parent]) -> list[dict[str, str]]:
    transforms = {
        "pyrrolidine_to_piperidine": "[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1>>[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1",
        "piperidine_to_pyrrolidine": "[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1>>[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1",
        "piperidine_to_morpholine": "[N:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1>>[N:1]1[CH2:2][CH2:3]O[CH2:5][CH2:6]1",
    }
    records: list[dict[str, str]] = []
    for parent in parents:
        for label, smarts in transforms.items():
            for product in reaction_products(parent.mol, smarts):
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "ring_size_and_ring_system_scans",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": label,
                        }
                    )
    return records


def enumerate_heteroatom_walks(parents: list[Parent]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    rxn = "[cH:1]>>[n:1]"
    for parent in parents:
        seen_for_parent: set[str] = set()
        for product in reaction_products(parent.mol, rxn):
            if product != parent.smiles and product not in seen_for_parent:
                seen_for_parent.add(product)
                records.append(
                    {
                        "method": "heteroatom_walks",
                        "parent_id": parent.representative_id,
                        "parent_smiles": parent.smiles,
                        "product_smiles": product,
                        "transform": "aromatic_CH_to_N",
                    }
                )
            if len(seen_for_parent) >= 40:
                break
    return records


def enumerate_scaffold_hops(parents: list[Parent], core_frags: list[tuple[str, Chem.Mol]]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    parent_core_map: dict[str, list[Chem.Mol]] = {}
    parent_side_map: dict[str, list[Chem.Mol]] = {}
    for parent in parents:
        parent_core_map[parent.representative_id] = []
        parent_side_map[parent.representative_id] = []
        for frag_smiles in BRICS.BRICSDecompose(parent.mol):
            frag = Chem.MolFromSmiles(frag_smiles)
            if frag is None:
                continue
            dummy_count = sum(atom.GetAtomicNum() == 0 for atom in frag.GetAtoms())
            ring_count = frag.GetRingInfo().NumRings()
            if dummy_count >= 2 and ring_count >= 1:
                parent_core_map[parent.representative_id].append(frag)
            elif dummy_count == 1:
                parent_side_map[parent.representative_id].append(frag)

    for parent in parents:
        side_frags = parent_side_map[parent.representative_id][:4]
        alt_cores = []
        for core_smiles, core in core_frags:
            native_smiles = {Chem.MolToSmiles(c) for c in parent_core_map[parent.representative_id]}
            if core_smiles not in native_smiles:
                alt_cores.append(core)
        for alt_core in alt_cores[:2]:
            products = limited_brics_build([alt_core, *side_frags], limit=20)
            for product in products:
                if product != parent.smiles:
                    records.append(
                        {
                            "method": "scaffold_hopping",
                            "parent_id": parent.representative_id,
                            "parent_smiles": parent.smiles,
                            "product_smiles": product,
                            "transform": "alternate_core_with_parent_sidechains",
                        }
                    )
    return records


def enumerate_fragment_growing(core_frags: list[tuple[str, Chem.Mol]], side_frags: list[tuple[str, Chem.Mol]]) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    limited_cores = [frag for _, frag in core_frags[:3]]
    limited_sides = [frag for _, frag in side_frags[:6]]
    for core in limited_cores:
        core_smiles = Chem.MolToSmiles(core)
        for side in limited_sides:
            products = limited_brics_build([core, side], limit=4)
            for product in products:
                records.append(
                    {
                        "method": "fragment_growing",
                        "parent_id": core_smiles,
                        "parent_smiles": core_smiles,
                        "product_smiles": product,
                        "transform": "core_plus_one_side_fragment",
                    }
                )
    return records


def deduplicate_records(records: list[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[tuple[str, str]] = set()
    deduped: list[dict[str, str]] = []
    for record in records:
        key = (record["method"], record["product_smiles"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(record)
    return deduped


def save_outputs(records: list[dict[str, str]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    library_path = OUTPUT_DIR / "enumeration_library_10_methods.csv"
    counts_path = OUTPUT_DIR / "enumeration_method_counts.csv"
    summary_path = OUTPUT_DIR / "enumeration_method_counts.md"

    library_df = pd.DataFrame(records).sort_values(["method", "product_smiles"]).reset_index(drop=True)
    library_df.to_csv(library_path, index=False)

    counts = (
        library_df.groupby("method")["product_smiles"]
        .nunique()
        .rename("unique_products")
        .sort_values(ascending=False)
        .reset_index()
    )
    counts.to_csv(counts_path, index=False)

    total_unique = library_df["product_smiles"].nunique()
    lines = [
        "# Enumeration Method Counts",
        "",
        f"- Total unique products across all methods: {total_unique}",
        "",
        "## Per-method counts",
        "",
    ]
    for row in counts.itertuples():
        lines.append(f"- `{row.method}`: {row.unique_products}")
    summary_path.write_text("\n".join(lines))


def main() -> None:
    parents = load_positive_parents()
    core_frags, side_frags = classify_brics_fragments(parents)

    all_records: list[dict[str, str]] = []
    all_records.extend(enumerate_r_group(parents, side_frags))
    all_records.extend(enumerate_reaction_based(parents))
    all_records.extend(enumerate_reagent_pool_combinatorics(core_frags, side_frags))
    all_records.extend(enumerate_matched_pairs(parents))
    all_records.extend(enumerate_bioisosteres(parents))
    all_records.extend(enumerate_linker_scans(parents))
    all_records.extend(enumerate_ring_scans(parents))
    all_records.extend(enumerate_heteroatom_walks(parents))
    all_records.extend(enumerate_scaffold_hops(parents, core_frags))
    all_records.extend(enumerate_fragment_growing(core_frags, side_frags))

    all_records = deduplicate_records(all_records)
    save_outputs(all_records)

    counts = defaultdict(int)
    unique_smiles = set()
    for record in all_records:
        counts[record["method"]] += 1
        unique_smiles.add(record["product_smiles"])

    print(f"total_unique_products={len(unique_smiles)}")
    for method in sorted(counts):
        print(f"{method}={counts[method]}")


if __name__ == "__main__":
    main()
