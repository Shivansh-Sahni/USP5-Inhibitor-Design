# Regression Features, Enumeration Methods, and Narrowing Pipeline

Generated on 2026-03-31.

## Project context

- Final potency model: saved `ExtraTreesRegressor`
- Final reported in-sample `R^2`: `0.893359`
- Broad enumerated library: `3264` unique compounds
- Final screening start set: `3274` compounds including original positives

## 1. Features used for regression

The final regression model used the `base_graph` feature block, meaning standard physicochemical descriptors plus graph-topology descriptors derived from RDKit.

- `mw` (Molecular weight): Tracks overall size and strongly affects permeability, potency trends, and developability.
- `logp` (logP): Captures lipophilicity, which often influences membrane penetration and nonspecific binding.
- `tpsa` (Topological polar surface area): Measures polarity and is useful for absorption and binding-environment balance.
- `hbd` (Hydrogen-bond donors): Counts donor groups that can drive binding or reduce permeability if too high.
- `hba` (Hydrogen-bond acceptors): Counts acceptor groups that affect recognition and physicochemical behavior.
- `rot` (Rotatable bonds): Approximates flexibility, which can affect entropy, oral exposure, and 3D fit.
- `rings` (Ring count): Tracks scaffold rigidity and shape complexity.
- `hac` (Heavy atom count): A simple size and complexity descriptor that often correlates with potency trends.
- `fsp3` (Fraction sp3 carbons): Measures saturation / 3D character and can help distinguish flat from more three-dimensional molecules.
- `bertz` (Bertz complexity): Summarizes structural complexity from the molecular graph.
- `balaban` (Balaban J index): A graph-connectivity descriptor that reflects topology.
- `chi0v` (Valence connectivity index Chi0v): Encodes atom connectivity at a local graph level.
- `chi1v` (Valence connectivity index Chi1v): Encodes one-bond connectivity patterns.
- `chi2v` (Valence connectivity index Chi2v): Encodes two-bond graph connectivity patterns.
- `kappa1` (Kier shape index Kappa1): Represents coarse molecular shape.
- `kappa2` (Kier shape index Kappa2): Represents shape and branching at a deeper level.
- `kappa3` (Kier shape index Kappa3): Represents higher-order shape complexity and was the top single feature in the final fit.

## 2. Enumeration techniques

The canonical enumeration workflow used 10 methods.

### R-group enumeration

- What it does: Keep a core fragment fixed and swap side fragments around it.
- How it was done in this project: Implemented with BRICS core fragments plus a pooled side-fragment set, rebuilding products while preserving the parent core.
- Unique products produced in the final broad library: `1205`

### Reaction-based enumeration

- What it does: Generate analogs only through specific chemistry transforms.
- How it was done in this project: Implemented with RDKit reaction SMARTS such as acid-to-amide conversions and phenol O-alkylation.
- Unique products produced in the final broad library: `31`

### Reagent-pool combinatorics

- What it does: Combine a shared pool of compatible core and side fragments to make many plausible analogs.
- How it was done in this project: Implemented with a limited BRICS fragment pool and BRICSBuild to create global combinations.
- Unique products produced in the final broad library: `1200`

### Matched molecular pair expansion

- What it does: Make small medicinal-chemistry edits such as halogen swaps or methoxy/hydroxyl interconversion.
- How it was done in this project: Implemented with targeted single-step SMARTS replacements like Cl to F, Br to Cl, and OMe to OH.
- Unique products produced in the final broad library: `25`

### Bioisosteric swaps

- What it does: Replace one functional group with a chemically related surrogate.
- How it was done in this project: Implemented with transforms such as acid to amide, acid to thioacid, and amide to thioamide.
- Unique products produced in the final broad library: `15`

### Linker scans

- What it does: Vary the connector between motifs to alter spacing, polarity, or flexibility.
- How it was done in this project: Implemented with predefined SMARTS transforms that insert oxygen or carbamate-like changes into existing linkers.
- Unique products produced in the final broad library: `5`

### Ring-size / ring-system scans

- What it does: Change ring size or ring composition to test shape and polarity changes.
- How it was done in this project: Implemented with SMARTS transforms such as pyrrolidine to piperidine and piperidine to morpholine.
- Unique products produced in the final broad library: `3`

### Heteroatom walks

- What it does: Move from aryl CH positions into heteroaryl nitrogens to tune polarity and recognition.
- How it was done in this project: Implemented as aromatic CH to N conversions, capped per parent to avoid runaway expansion.
- Unique products produced in the final broad library: `61`

### Scaffold hopping

- What it does: Preserve side chains but replace the central core with an alternative BRICS core.
- How it was done in this project: Implemented by pairing parent-derived side fragments with non-native core fragments from other positives.
- Unique products produced in the final broad library: `671`

### Fragment growing

- What it does: Start from a core fragment and add one side fragment outward.
- How it was done in this project: Implemented with limited BRICS builds of one core plus one side fragment to create smaller controlled expansions.
- Unique products produced in the final broad library: `430`

## 3. Narrowing-down techniques

The final narrowing funnel starts with potency prediction and then applies increasingly strict chemistry, ADMET, and 3D plausibility filters.

### Model-first potency screen

The saved ExtraTrees regression model scores the full library first, and only compounds with predicted pIC50 >= 4.60 and acceptable descriptor-space distance continue.

### PAINS and BRENK alert removal

Compounds with nuisance or problematic substructure alerts are removed before deeper triage.

### Lipinski-style developability filter

The pipeline keeps compounds with <= 1 Lipinski violation to avoid drifting too far into poor oral-drug-like space.

### Topological polar surface area filter

TPSA is constrained to 40-140 A^2 to balance permeability with needed polarity.

### Molecular surface area filter

Labute approximate surface area is constrained to 130-235 A^2 to avoid structures that are too small or too bulky for the intended pocket and exposure profile.

### Molecular weight filter

MW is constrained to 250-550 to remove very small fragments and oversized analogs.

### Flexibility and charge filter

Rotatable bonds are capped at 10 and absolute formal charge at 1 to keep molecules closer to plausible lead space.

### ADMET-AI multigate

The surviving set is screened with ADMET-AI for AMES, hERG, ClinTox, HIA, oral bioavailability, Caco-2, and solubility-related properties.

### Multistructure 3D binding plausibility

Compounds are compared against several USP5 ZnF-UBD co-crystal-inspired templates to assess binding-score, shape overlap, pharmacophore match, and steric clashes.

### Portfolio selection

The last step separates strict primary leads from orthogonal backups so the project does not overcommit to only one chemotype.

## 4. Final stage counts

- `start_broad_enumeration_plus_original_positives`: `3274`
- `stage1_original_final_model_predicted_pIC50_ge_4.60`: `1272`
- `stage2_pains_brenk_free`: `772`
- `stage3_lipinski_psa_surface_area_flexibility`: `59`
- `stage4_admet_ai_multigate`: `24`
- `stage5_multistructure_template_docking_strict`: `4`
- `stage5b_multistructure_template_docking_relaxed`: `5`

## 5. Verbal summary

The regression model uses interpretable physicochemical and graph-shape descriptors, the enumeration layer spans 10 complementary medicinal-chemistry-inspired expansion methods, and the narrowing layer combines potency prediction, property filtering, ADMET-AI, and multistructure 3D binding plausibility. Together, these steps mean most of the core computational pipeline has already been completed.
