from __future__ import annotations

from pathlib import Path
from textwrap import wrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


OUTPUT_DIR = Path("outputs")
PDF_PATH = OUTPUT_DIR / "enumeration_techniques_python_guide.pdf"
MD_PATH = OUTPUT_DIR / "enumeration_techniques_python_guide.md"


TECHNIQUES: list[dict[str, object]] = [
    {
        "name": "R-Group Enumeration",
        "what": "Keep a core scaffold fixed and systematically substitute allowed groups at one or more attachment points.",
        "python": [
            "Represent the scaffold with labeled attachment atoms such as [*:1], [*:2].",
            "Store allowed substituent pools as SMILES strings or reaction-ready fragments.",
            "Loop over all allowed combinations and attach them to the scaffold.",
            "Canonicalize, deduplicate, and record the parent scaffold plus substituent IDs.",
        ],
        "why": "This is the most direct way to explore local SAR around a known active core.",
        "snippet": [
            "scaffold = 'core-[*:1]-[*:2]'",
            "r1_pool = ['F', 'Cl', 'CN']",
            "r2_pool = ['piperidine', 'morpholine', 'pyrrolidine']",
            "for r1 in r1_pool:",
            "    for r2 in r2_pool:",
            "        analog = attach_groups(scaffold, {1: r1, 2: r2})",
            "        save(analog)",
        ],
    },
    {
        "name": "Reaction-Based Enumeration",
        "what": "Generate compounds only through plausible synthetic reactions such as amide coupling, reductive amination, or Suzuki coupling.",
        "python": [
            "Define reactions as transforms or SMARTS-like templates.",
            "Load reagent lists for each reaction partner.",
            "Apply each reaction to compatible reagent pairs.",
            "Filter invalid products and annotate the reaction used.",
        ],
        "why": "This keeps the library grounded in chemistry that can actually be made.",
        "snippet": [
            "amide_rxn = define_reaction('acid + amine -> amide')",
            "for acid in acid_pool:",
            "    for amine in amine_pool:",
            "        product = run_reaction(amide_rxn, acid, amine)",
            "        if product is not None:",
            "            save(product, route='amide_coupling')",
        ],
    },
    {
        "name": "Reagent-Pool Combinatorics",
        "what": "Build large libraries by crossing compatible sets of reagents across one or more steps.",
        "python": [
            "Create reagent tables grouped by role such as acids, amines, aryl halides, boronic acids.",
            "Define which pools can be crossed in each step.",
            "Enumerate all compatible combinations.",
            "Track reagent lineage so every product can be traced back to building blocks.",
        ],
        "why": "This is the fastest route from a few parent chemotypes to thousands of virtual analogs.",
        "snippet": [
            "for acid in acid_pool:",
            "    for amine in amine_pool:",
            "        intermediate = couple(acid, amine)",
            "        for aryl_halide in aryl_halide_pool:",
            "            final = diversify(intermediate, aryl_halide)",
            "            save(final)",
        ],
    },
    {
        "name": "Matched Molecular Pair Expansion",
        "what": "Apply small medicinal chemistry edits that are known to change potency, polarity, or selectivity in interpretable ways.",
        "python": [
            "Create a library of local transforms such as F to CN, phenyl to pyridyl, methyl to cyclopropyl.",
            "Search each parent molecule for transformable sites.",
            "Apply one change at a time to generate local analogs.",
            "Record the exact transform used for later SAR analysis.",
        ],
        "why": "This explores nearby chemical space without drifting too far from active parents.",
        "snippet": [
            "transforms = [('F', 'CN'), ('phenyl', 'pyridyl')]",
            "for mol in parents:",
            "    for old, new in transforms:",
            "        for analog in apply_local_transform(mol, old, new):",
            "            save(analog, transform=f'{old}->{new}')",
        ],
    },
    {
        "name": "Bioisosteric Swaps",
        "what": "Replace a functional group with a chemically different group that can play a similar role in binding or properties.",
        "python": [
            "Build dictionaries of common bioisostere replacements.",
            "Match eligible motifs in each molecule.",
            "Substitute the motif while preserving attachment geometry where possible.",
            "Recompute descriptors so the property effect is visible.",
        ],
        "why": "Bioisosteres are a standard way to improve potency, stability, permeability, or safety.",
        "snippet": [
            "bioisosteres = {'carboxylic_acid': ['tetrazole', 'acylsulfonamide']}",
            "for mol in parents:",
            "    for motif, replacements in bioisosteres.items():",
            "        for repl in replacements:",
            "            analogs = swap_motif(mol, motif, repl)",
            "            save_all(analogs)",
        ],
    },
    {
        "name": "Linker Scans",
        "what": "Change the connector between two motifs by varying length, flexibility, saturation, or heteroatom content.",
        "python": [
            "Define the two fixed endpoint motifs.",
            "Create a library of linker fragments such as CH2, O, NH, CH2CH2, piperazine.",
            "Attach each linker between the same endpoints.",
            "Store linker identity as a separate annotation field.",
        ],
        "why": "Linkers strongly affect geometry, entropy, polarity, and binding presentation.",
        "snippet": [
            "linkers = ['CH2', 'O', 'NH', 'CH2CH2', 'piperazine']",
            "for linker in linkers:",
            "    analog = connect(left_motif, linker, right_motif)",
            "    save(analog, linker=linker)",
        ],
    },
    {
        "name": "Ring-Size And Ring-System Scans",
        "what": "Replace one ring with a related ring system such as pyrrolidine to piperidine or phenyl to pyridyl.",
        "python": [
            "Define a ring replacement table for each chemotype.",
            "Match the ring position in the parent scaffold.",
            "Swap in alternate ring systems while preserving vectors.",
            "Check resulting valence and geometry sanity before saving.",
        ],
        "why": "Ring changes often shift potency, selectivity, and developability in a controlled way.",
        "snippet": [
            "ring_swaps = ['pyrrolidine', 'piperidine', 'morpholine']",
            "for ring in ring_swaps:",
            "    analog = replace_ring(parent, target_site='amine_ring', new_ring=ring)",
            "    save(analog)",
        ],
    },
    {
        "name": "Heteroatom Walks",
        "what": "Move or swap heteroatoms within a scaffold to tune electronics, polarity, hydrogen bonding, and binding interactions.",
        "python": [
            "Identify atom positions where carbon, nitrogen, oxygen, or sulfur swaps are reasonable.",
            "Generate all allowed single-site heteroatom variants.",
            "Discard unstable or chemically nonsensical products.",
            "Track which atom and position changed.",
        ],
        "why": "A heteroatom walk is a compact way to probe both binding and property shifts.",
        "snippet": [
            "for position in editable_positions(parent):",
            "    for atom in ['C', 'N', 'O', 'S']:",
            "        analog = mutate_atom(parent, position, atom)",
            "        if is_reasonable(analog):",
            "            save(analog)",
        ],
    },
    {
        "name": "Scaffold Hopping",
        "what": "Keep the key pharmacophore pattern but replace the central core with a different scaffold.",
        "python": [
            "Define the required pharmacophore points or anchor vectors.",
            "Create a set of alternative cores that present those vectors similarly.",
            "Attach the same side chains to each alternate core.",
            "Cluster and compare the resulting chemotypes by novelty and properties.",
        ],
        "why": "This is how a project escapes from one chemotype while preserving the core hypothesis.",
        "snippet": [
            "for new_core in alternate_cores:",
            "    hopped = transplant_sidechains(new_core, sidechains=parent_sidechains)",
            "    save(hopped, strategy='scaffold_hop')",
        ],
    },
    {
        "name": "Fragment Growing",
        "what": "Start from a smaller active motif and extend it into nearby space with additional fragments.",
        "python": [
            "Choose a minimal active fragment or anchor motif.",
            "Define growth vectors and allowed fragment additions.",
            "Add one fragment at a time, then optionally iterate to a second growth round.",
            "Score each product for size, polarity, and tractability as the fragment grows.",
        ],
        "why": "Fragment growing is useful when you know a minimal motif but need more interactions and potency.",
        "snippet": [
            "seed = 'minimal_active_fragment'",
            "for frag in fragment_pool:",
            "    analog = grow(seed, frag, vector='exit_1')",
            "    save(analog)",
        ],
    },
]


def _wrap_lines(lines: list[str], width: int = 92) -> list[str]:
    wrapped: list[str] = []
    for line in lines:
        if not line:
            wrapped.append("")
            continue
        wrapped.extend(wrap(line, width=width) or [""])
    return wrapped


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
        for line in wrap(raw_line, width=95) or [""]:
            ax.text(0, y, line, fontsize=fontsize, va="top", family="DejaVu Sans Mono" if line.startswith("    ") else "DejaVu Sans")
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


def build_markdown() -> str:
    lines = [
        "# Enumeration Techniques And How They Would Be Done With Python",
        "",
        "This guide explains the 10 main enumeration techniques selected for the USP5 virtual library strategy.",
        "",
    ]
    for idx, technique in enumerate(TECHNIQUES, start=1):
        lines.append(f"## {idx}. {technique['name']}")
        lines.append("")
        lines.append(f"**What it is:** {technique['what']}")
        lines.append("")
        lines.append(f"**Why use it:** {technique['why']}")
        lines.append("")
        lines.append("**How it would be done with Python:**")
        for step in technique["python"]:
            lines.append(f"- {step}")
        lines.append("")
        lines.append("**Minimal pseudocode:**")
        lines.append("```python")
        lines.extend(technique["snippet"])
        lines.append("```")
        lines.append("")
    return "\n".join(lines)


def build_pdf() -> None:
    with PdfPages(PDF_PATH) as pdf:
        intro_lines = _wrap_lines(
            [
                "This PDF describes the 10 main virtual-enumeration techniques selected for the USP5 project and how each one would be carried out with Python.",
                "",
                "The techniques are:",
                *[f"- {item['name']}" for item in TECHNIQUES],
                "",
                "These are implementation-oriented descriptions, not full production code. The idea is to show the role of each technique and the Python workflow behind it.",
            ]
        )
        _add_text_page(pdf, "Enumeration Techniques Overview", intro_lines)

        for idx, technique in enumerate(TECHNIQUES, start=1):
            lines = [
                f"What it is: {technique['what']}",
                "",
                f"Why use it: {technique['why']}",
                "",
                "How it would be done with Python:",
                *[f"- {step}" for step in technique["python"]],
                "",
                "Minimal pseudocode:",
                *[f"    {line}" for line in technique["snippet"]],
            ]
            _add_text_page(pdf, f"{idx}. {technique['name']}", lines)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    MD_PATH.write_text(build_markdown())
    build_pdf()
    print(f"Wrote {MD_PATH}")
    print(f"Wrote {PDF_PATH}")


if __name__ == "__main__":
    main()
