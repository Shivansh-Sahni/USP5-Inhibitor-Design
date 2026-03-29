# Full Positive Expansion And Lead-Narrowing Strategy

## Scope

This document is a **strategy report**, not an implementation. It starts from **all positive molecules already present in the training database**, describes how to expand them into a **very large combinatorial chemistry database**, explains how to improve that database after expansion, and then catalogs the major techniques that can be used to **narrow from a huge library to lead candidates**.

The current project dataset has:

- `19` unique molecules in the deduplicated modeling table
- `13` positive molecules
- `10` positives with measured `ic50 > 0`
- `3` positives with `ic50 == 0` meaning active but no numeric IC50
- `6` inactive molecules with `ic50 == -1`

For this report, a molecule is treated as **positive** if it has either:

- a measured active row
- an active/no-numeric-IC50 row

## Positive Starting Set

### Measured positives

| ID | Target pIC50 | Note |
| --- | ---: | --- |
| `CHEMBL1923233` | `5.52` | potent measured singleton |
| `CHEMBL2012938` | `5.23` | measured active with mixed duplicate status |
| `CHEMBL5278336` | `5.15` | measured active with conflicting duplicates |
| `CHEMBL5410606` | `5.10` | tractable measured singleton |
| `CHEMBL4635160` | `5.00` | potent but structurally complex |
| `CHEMBL5276142` | `4.67` | compact analog-series member |
| `CHEMBL1410015` | `4.62` | compact analog-series member |
| `CHEMBL1493046` | `4.54` | compact analog-series member |
| `CHEMBL5267683` | `3.70` | weaker measured active |
| `CHEMBL4129140` | `3.70` | weaker measured active |

### Active without numeric IC50

| ID | Assigned pIC50 | Note |
| --- | ---: | --- |
| `25191039` | `5.0` | active label only; polyphenolic chemotype |
| `70691637` | `5.0` | active label only; very large conjugate |
| `145963683` | `5.0` | active label only; very large conjugate |

## Parent Families To Expand

The expansion should begin from **all 13 positives**, but not all parents should contribute equally. The right approach is to keep every positive in the parent pool while assigning a **generation weight** based on confidence, tractability, and medicinal chemistry usefulness.

### Family A: compact analog-ready heteroaryl amide series

Parents:

- `CHEMBL5276142`
- `CHEMBL1410015`
- `CHEMBL1493046`

Why it matters:

- This is the strongest local analog family in the dataset.
- It already shows real within-project variation.
- It has obvious R-group vectors and good property headroom.

Generation weight:

- **very high**

### Family B: sulfonamide-piperidine nicotinamide-like series

Parent:

- `CHEMBL5278336`

Why it matters:

- Potent measured activity
- multiple medicinal chemistry vectors
- good polarity balance

Generation weight:

- **high**

### Family C: benzylamide-fused heterocycle acid series

Parent:

- `CHEMBL5410606`

Why it matters:

- good potency
- acceptable physicochemical window
- several straightforward substitution vectors

Generation weight:

- **high**

### Family D: cyano-enamide / potentially reactive motif series

Parent:

- `CHEMBL1923233`

Why it matters:

- strongest measured potency
- clean measured duplicate agreement

Caution:

- may carry electrophilic or Michael-acceptor-like liability

Generation weight:

- **moderate**, but handled in a separate liability-aware branch

### Family E: polyphenolic and acylated polyphenolic actives

Parents:

- `CHEMBL2012938`
- `25191039`

Why it matters:

- both are positive
- they may represent a phenolic recognition motif worth testing

Caution:

- likely alert-prone
- lipophilic and/or heavy
- risk of redox or assay-interference behavior

Generation weight:

- **moderate to low**

### Family F: bulky or structurally complex actives

Parents:

- `CHEMBL4635160`
- `CHEMBL5267683`
- `CHEMBL4129140`

Why it matters:

- these provide orthogonal chemical space
- they may capture binding hypotheses missed by simpler scaffolds

Caution:

- higher complexity
- less favorable starting properties in some cases
- fewer local analog precedents

Generation weight:

- **moderate**

### Family G: very large conjugated actives

Parents:

- `70691637`
- `145963683`

Why it matters:

- they are part of the positive set and could encode useful pharmacophore combinations

Caution:

- highly non-lead-like
- likely poor as primary lead series
- should mainly inform motif extraction, fragment transfer, and pharmacophore mining

Generation weight:

- **low for direct analog enumeration**
- **high for motif extraction**

## Full Combinatorial Chemistry Expansion Pipeline

The goal is not to make one giant uncontrolled pile of molecules. The goal is to build a **large, information-rich, synthesis-aware library** that preserves the signal from the starting actives while exploring adjacent chemical space intelligently.

### Stage 1: parent standardization and confidence tagging

For each of the 13 positives:

- standardize SMILES
- remove salts and normalize charge states
- assign canonical tautomer state
- annotate parent family
- tag measurement confidence
- tag duplicate conflict status
- tag tractability level
- tag risk level such as reactive motif, polyphenol, oversized, or conjugate

Output:

- a curated parent table with confidence and liability metadata

### Stage 2: scaffold extraction and attachment-vector mapping

For each parent:

- define the core scaffold
- identify attachment vectors
- identify essential versus optional substituents
- isolate ring system, linker, and terminal-group vectors
- identify any likely pharmacophore anchors

Key vector classes:

- aryl substitution
- heteroaryl substitution
- cyclic amine replacement
- linker length and composition
- amide, sulfonamide, urea, carbamate, ester, and reversed-amide swaps
- acid and acid-isostere replacement
- ring contraction, ring expansion, and rigidification
- heteroatom walk
- halogen walk
- nitrile, methoxy, hydroxyl, amino, and small alkyl scans

Output:

- one vector map per parent family

### Stage 3: reaction-aware enumeration

Enumerate analogs using chemistry that a medicinal chemistry team could plausibly make. This is where the library becomes large.

Reaction-aware routes should include:

- amide coupling
- sulfonamide formation
- reductive amination
- SNAr
- Suzuki coupling
- Buchwald-Hartwig amination
- urea and carbamate formation
- heterocycle substitution
- ether formation
- alkylation of secondary amines
- esterification and hydrolysis pairs
- acid-isostere installation

The important rule is:

- only enumerate transformations that correspond to a plausible synthetic disconnection

This keeps the library large but still useful.

### Stage 4: medicinal chemistry transform expansion

After direct reaction-based enumeration, expand each family with known medicinal chemistry transforms:

- matched molecular pair transforms
- bioisosteric replacements
- ring opening and ring closure
- conformational restriction
- polarity tuning
- lipophilicity tuning
- acidity and basicity tuning
- aromatic to saturated replacement
- scaffold decoration
- scaffold simplification
- scaffold hopping
- vector inversion
- stereochemical expansion where justified

This is the stage that converts a decent library into a genuinely rich design space.

### Stage 5: motif extraction from large or unattractive positives

Do not throw away large or awkward positives. Instead:

- extract the minimum motif likely responsible for activity
- split conjugates into pharmacophore-bearing and delivery-like regions
- move useful motifs into smaller lead-like scaffolds
- transfer promising aryl, heteroaryl, linker, or cationic features into tractable cores

This is especially important for:

- `70691637`
- `145963683`
- `CHEMBL2012938`

These parents are more valuable as **motif donors** than as direct lead series.

### Stage 6: external analog and idea harvesting

To strengthen the database beyond purely internal enumeration, add conceptually aligned compounds from outside the current training set:

- nearest analog mining from public bioactivity databases
- patent analog mining
- make-on-demand vendor space analog mining
- fragment library overlay
- known USP or DUB chemotype mining
- ligand-based similarity expansion from the current positives
- pharmacophore-based analog harvesting

This stage improves both diversity and realism.

### Stage 7: AI-assisted proposal generation

Use AI as a proposal layer, not as the only generator.

Useful AI roles:

- propose substituent sets for each attachment vector
- suggest scaffold hops
- suggest bioisosteres
- suggest matched-pair-inspired edits
- identify likely dead-end chemotypes to avoid
- generate lead-like motif compressions from oversized positives
- propose routes and building-block-compatible analogs

AI should be constrained by:

- synthetic accessibility
- parent family rules
- property windows
- forbidden liabilities

### Stage 8: database consolidation

Once raw generation is complete:

- canonicalize
- deduplicate
- remove impossible structures
- remove unstable valence states
- collapse redundant stereoisomers if not relevant yet
- annotate parent lineage
- annotate transform type
- annotate synthetic route class
- annotate estimated tractability

At this point the library becomes a true working database rather than a raw enumeration dump.

## How To Reach Thousands And Thousands Of Molecules

The parent pool already supports a multi-thousand library if expansion is staged correctly.

### Reasonable scale-up model

A practical way to reach a large database is:

- `Family A`: `2,000` to `8,000` compounds
- `Family B`: `1,000` to `4,000` compounds
- `Family C`: `1,000` to `4,000` compounds
- `Family D`: `500` to `2,000` compounds
- `Family E`: `500` to `2,000` compounds
- `Family F`: `1,000` to `3,000` compounds
- `Family G`: `500` to `2,000` compounds direct plus many motif-transfer ideas

That already yields a raw design space on the order of:

- roughly `6,500` to `25,000+` discrete proposals

If the project wants a genuinely huge virtual space, additional combinatorial breadth can push this further:

- multiple reagent pools per vector
- multi-vector cross-combination
- scaffold hops per family
- motif transfer from bulky positives into compact cores
- make-on-demand building-block enumeration

That easily moves the project into:

- `25,000` to `100,000+` virtual compounds

For a small-data project, that is already more than enough. Bigger than that is only useful if the ranking and triage stack is strong.

## How To Improve The Database After Initial Expansion

Once the first large library exists, the next goal is not just more compounds. The goal is a **stronger** database.

### 1. Improve label quality

- separate measured positives from assigned actives in every downstream view
- down-weight uncertain positives
- flag duplicate conflicts
- add assay provenance
- distinguish orthogonal confirmation from single-assay signal

### 2. Improve parent diversity

- make sure every positive family contributes analogs
- prevent one easy chemistry series from dominating the database
- intentionally preserve orthogonal chemotypes

### 3. Improve tractability

- prefer routes with short synthesis paths
- track availability of reagents
- track make-on-demand feasibility
- remove chemistry that depends on unrealistic disconnections

### 4. Improve novelty

- cluster the library
- measure parent-distance and scaffold-distance
- cap oversampling of near-clones
- add scaffold hops deliberately

### 5. Improve property balance

- tune molecular weight distribution
- tune cLogP distribution
- tune TPSA distribution
- tune aromaticity and `Fsp3`
- tune rotatable bond counts

### 6. Improve liability awareness

- remove obvious structural alerts
- remove likely aggregators
- remove unstable or reactive motifs unless intentional
- explicitly tag covalent or electrophilic chemistry

### 7. Improve realism

- include tautomers and protonation-aware scoring metadata
- model stereochemistry when relevant
- include synthetic route annotations
- include building-block IDs where possible

### 8. Improve data richness

- attach descriptor blocks
- attach fingerprints
- attach scaffold IDs
- attach family IDs
- attach transform provenance
- attach route class
- attach predicted ADME/Tox fields
- attach docking or pharmacophore scores if later computed

### 9. Improve the active-learning loop

- use early test data to retrain prioritization
- enrich underexplored series
- penalize repeatedly failing transform classes
- expand around real SAR once experiments start coming back

### 10. Improve balance between hypothesis styles

The best large database should not be only one kind of chemistry. It should include:

- close analogs
- bioisosteres
- scaffold hops
- polarity-tuned variants
- rigidified analogs
- solubility-rescue analogs
- selectivity-improving analogs
- motif-transferred designs

## Techniques To Improve The Database Even Further

If the goal is the strongest possible pre-screening database, these are the main technique classes to consider.

### Chemistry-generation techniques

- classical R-group enumeration
- reaction-based enumeration
- reagent-pool combinatorics
- matched molecular pair expansion
- bioisosteric swap libraries
- linker scans
- ring-size scans
- heteroatom walks
- halogen walks
- conformational lock libraries
- scaffold-hopping libraries
- fragment growing
- fragment merging
- fragment linking
- pharmacophore grafting
- motif transfer from non-lead-like actives
- stereoisomer enumeration
- tautomer-aware enumeration where relevant

### Knowledge-driven augmentation techniques

- literature analog mining
- patent analog mining
- public bioactivity neighbor mining
- known-target-family chemotype mining
- competitor chemotype abstraction
- mechanism-informed hypothesis design

### Computation-driven augmentation techniques

- ligand similarity expansion
- pharmacophore expansion
- 3D shape expansion
- docking-seeded analog generation
- interaction-fingerprint-guided enumeration
- QSAR-guided generation
- uncertainty-guided generation
- generative model proposal layers
- retrosynthesis-constrained generation

### Portfolio-improvement techniques

- novelty quotas
- family-balance quotas
- route-feasibility quotas
- property-distribution balancing
- liability balancing
- diversity clustering
- active-learning refresh cycles

## Narrowing To Lead: Report Only, Do Not Implement

The database should be narrowed through a staged funnel. No single filter is enough. The best practice is to combine chemistry sense, developability, target relevance, and experimental risk.

## Stage-Wise Narrowing Funnel

### Stage 1: structural sanity and database hygiene

Remove:

- duplicates
- impossible valence states
- unstable or nonsensical protonation states
- clearly broken structures
- compounds outside intended chemistry space

### Stage 2: medicinal chemistry hard filters

Use property and chemistry guardrails such as:

- molecular weight limits
- cLogP or cLogD limits
- TPSA limits
- HBD and HBA limits
- rotatable bond limits
- aromatic ring count limits
- heavy atom count limits
- charge-state sanity
- ring strain sanity
- chemical stability sanity

### Stage 3: structural alert and nuisance filters

Flag or remove:

- PAINS motifs
- frequent hitters
- colloidal aggregators
- redox cyclers
- quinones
- catechols when problematic
- Michael acceptors when not intentional
- unstable esters or acylating motifs
- metal chelators when unwanted
- thiol-reactive motifs
- promiscuous electrophiles

### Stage 4: lead-likeness and developability filters

Prefer compounds with:

- balanced polarity
- manageable lipophilicity
- acceptable size
- reasonable solubility outlook
- manageable flexibility
- acceptable 3D character
- low unnecessary complexity

### Stage 5: synthetic feasibility filters

Keep compounds that are:

- synthesizable in a short route
- compatible with available reagents
- low to moderate cost
- likely to be stable during synthesis and storage
- easy to resupply and analog around

### Stage 6: family and diversity control

Do not let one close analog cloud dominate the shortlist.

Apply:

- scaffold clustering
- diversity picking
- per-family caps
- per-transform caps
- novelty quotas
- redundancy penalties

### Stage 7: target-relevance scoring

Rank by target-focused evidence such as:

- similarity to the strongest measured actives
- pharmacophore match
- docking score
- interaction fingerprint quality
- 3D shape overlay
- binding-site complementarity
- target-state compatibility
- water-network compatibility
- key residue engagement hypotheses

### Stage 8: model-based prioritization

Use predictive models carefully:

- QSAR potency prediction
- consensus scoring across multiple models
- uncertainty-aware ranking
- applicability-domain penalties
- conformal or calibrated confidence scoring
- ensemble rank aggregation

Because the dataset is small, model output should be used as one signal among many, not as the only decider.

### Stage 9: ADME prediction narrowing

Apply in silico screens for:

- solubility
- permeability
- microsomal stability
- hepatocyte stability
- plasma protein binding
- CYP inhibition and induction risk
- transporter liability
- hERG risk
- general ion-channel risk
- clearance risk
- metabolite liability

### Stage 10: toxicity and liability narrowing

Apply screens for:

- Ames mutagenicity risk
- genotoxicity
- hepatotoxicity risk
- phospholipidosis risk
- mitochondrial toxicity risk
- phototoxicity
- reactive metabolite risk
- off-target promiscuity
- covalent liability if not intended

### Stage 11: selectivity narrowing

Especially for USP5 work, narrow further using:

- related deubiquitinase selectivity predictions
- homolog panel strategy
- cysteine-reactivity controls
- counterscreens against assay-specific nuisance proteins
- broader off-target family awareness

### Stage 12: assay-risk narrowing

Prioritize compounds that are less likely to generate misleading hits:

- low fluorescence interference risk
- low quenching risk
- low aggregation risk
- low redox interference risk
- low detergent sensitivity
- orthogonal-assay compatibility

### Stage 13: experiment design narrowing

Before synthesis or purchase, prioritize compounds that maximize learning:

- series leaders
- edge-of-SAR compounds
- property-rescue analogs
- selectivity probes
- mechanism probes
- orthogonal chemotypes
- positive controls
- negative controls

### Stage 14: multi-parameter optimization ranking

Final lead selection should be based on a balanced score, not potency alone.

Common MPO dimensions:

- potency
- novelty
- tractability
- property balance
- selectivity
- ADME outlook
- tox outlook
- assay robustness
- patentability
- cost and cycle time

## All Major Narrowing Techniques To Consider

This is the broad practical menu of narrowing methods to consider after the large database exists.

### Chemistry and rules-based techniques

- Lipinski-style filters
- Veber-style filters
- lead-likeness windows
- hard alert exclusion
- soft alert review
- reactivity review
- acid-base balance review
- ring-system review

### Similarity and diversity techniques

- nearest-neighbor ranking
- cluster-center selection
- max-min diversity picking
- scaffold diversity picking
- family-balanced picking
- novelty-distance picking

### Ligand-based scoring techniques

- 2D similarity
- 3D similarity
- shape overlays
- pharmacophore matching
- matched-pair extrapolation
- local-SAR interpolation

### Structure-based techniques

- rigid docking
- flexible docking
- ensemble docking
- induced-fit docking
- interaction fingerprint scoring
- water-aware scoring
- MM-GBSA style rescoring
- free-energy methods for top subsets

### Machine-learning techniques

- QSAR ranking
- ensemble models
- classification plus regression stacking
- uncertainty-aware models
- active learning
- Bayesian optimization over chemical space
- generative reranking
- applicability-domain filtering

### Developability techniques

- solubility prediction
- permeability prediction
- metabolic stability prediction
- clearance prediction
- protein binding prediction
- transporter liability prediction
- formulation-risk screening

### Safety techniques

- hERG prediction
- Ames prediction
- reactive metabolite screening
- tox panel prediction
- off-target panel prediction
- chemotype blacklists

### Project-level decision techniques

- patent-space review
- reagent-availability review
- synthesis-cost review
- cycle-time review
- portfolio-balance review
- orthogonal-hypothesis coverage review

## Best Practical Recommendation

If this strategy were executed, the best workflow would be:

1. Start from **all 13 positives** as the parent universe.
2. Expand most aggressively around the tractable measured families, especially the compact analog-ready series and the better-balanced measured singletons.
3. Use oversized or awkward positives mainly as **motif donors** and pharmacophore sources rather than as primary lead series.
4. Build a first serious library in the range of **10,000 to 30,000** curated virtual compounds.
5. If needed, extend to **25,000 to 100,000+** only after the ranking stack, synthesis rules, and liability filters are mature enough to handle that scale.
6. Narrow using a **multi-stage funnel** that combines chemistry rules, diversity control, target relevance, predictive models, ADME/Tox risk, and synthesis practicality.

## Most Important Warning

The current dataset is still very small and mixes measured values with assigned labels. That means the expansion can absolutely be ambitious, but the ranking logic must stay honest:

- measured actives should drive confidence
- assigned actives should guide hypothesis generation, not dominate prioritization
- duplicate conflicts should reduce confidence
- huge libraries are only valuable if the narrowing funnel is disciplined

## One-Sentence Summary

Use **all 13 positive molecules** as the starting universe, expand them through a **reaction-aware, medicinal-chemistry-guided combinatorial pipeline** into a **curated multi-thousand to low-six-figure virtual database**, strengthen that database through diversity, tractability, liability, and data-quality improvements, and then narrow to lead through a **multi-parameter staged funnel** rather than any single score or filter.
