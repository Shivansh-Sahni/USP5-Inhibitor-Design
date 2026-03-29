# USP5 Enumeration Seed Recommendation

## Bottom line

If the goal is **combinatorial chemistry followed by narrowing**, the best current **base series for enumeration** is:

1. **Primary seed series:** `CHEMBL5276142`, `CHEMBL1410015`, `CHEMBL1493046`
2. **Secondary seed:** `CHEMBL5278336`
3. **Secondary seed:** `CHEMBL5410606`
4. **Conditional seed only:** `CHEMBL1923233`

If you want a **single first base molecule**, start from **`CHEMBL5276142`**.

## Why `CHEMBL5276142` should go first

This compound sits inside the only compact local analog set that already shows usable SAR-like variation in the current dataset:

- `CHEMBL5276142` | pIC50 `4.67` | MW `337.4` | cLogP `2.61` | TPSA `69.3`
- `CHEMBL1410015` | pIC50 `4.62` | MW `300.4` | cLogP `3.51` | TPSA `25.2`
- `CHEMBL1493046` | pIC50 `4.54` | MW `314.4` | cLogP `3.90` | TPSA `25.2`

Why this is the best enumeration base:

- It already has nearby analog support in the dataset, instead of being a potency singleton.
- The series has reasonable property headroom for analog growth.
- It exposes obvious medicinal chemistry vectors:
  - terminal aryl substituent
  - cyclic amine identity and ring size
  - amine polarity / hydroxylation state
  - linker substitution around the carbonyl-methylene region
- `CHEMBL5276142` is the best member of this local series and has the healthiest polarity balance.

## Secondary seeds worth focused enumeration

### `CHEMBL5278336`

- pIC50 `5.15`
- MW `437.9`
- cLogP `2.12`
- TPSA `116.7`

Why keep it:

- Potency is decent.
- Physicochemical profile is still usable.
- Several obvious vectors exist:
  - aryl group on the piperidine
  - sulfonamide region
  - heteroaryl-carboxamide region

Why not make it seed number 1:

- It does not yet have local analog support in this dataset.
- The measured duplicates are inconsistent, so the SAR confidence is weaker than the raw potency suggests.

### `CHEMBL5410606`

- pIC50 `5.10`
- MW `429.9`
- cLogP `2.39`
- TPSA `110.5`

Why keep it:

- Potency and property window are both acceptable.
- The scaffold looks more tractable than the very large or highly fused alternatives.
- Enumeration vectors are clear enough:
  - anisyl / benzylamine region
  - acidic side chain length and replacement
  - aryl halogen / heteroaryl substitutions

Why it stays secondary:

- It is still a singleton in the current data, so you do not have real internal SAR yet.

## Conditional seed only

### `CHEMBL1923233`

- pIC50 `5.52`
- MW `384.3`
- cLogP `4.41`
- TPSA `65.8`

This is the most potent measured scaffold and the duplicate measurements are consistent.

Do **not** make it the default base series unless you intentionally want to pursue a potentially reactive motif. The `NC(=O)/C(C#N)=C/Ar` region looks like an **electrophile / Michael-acceptor-like liability** and should be treated carefully unless that chemistry is deliberate.

Use it only if:

- you want a focused liability-managed series
- you can confirm the mechanism justifies that motif
- you are comfortable running a reactivity / selectivity screen early

## Seeds to park for now

### `CHEMBL2012938`

- pIC50 `5.23`, but MW `562.6` and cLogP `6.14`
- Polyphenolic / oversized profile
- More liability risk and less clean combinatorial tractability

### `CHEMBL4635160`

- pIC50 `5.00`, but MW `524.1` and cLogP `5.44`
- Structurally complex and less attractive for broad enumeration

### `CHEMBL5267683` and `CHEMBL4129140`

- Both are weaker (`pIC50 3.7`)
- Both are already lipophilic
- Neither is the best place to spend the first enumeration cycle

### Label-only large conjugates and biotin-linked molecules

- Do not use these as initial enumeration seeds
- They are better treated as follow-up probes or special-purpose chemistry, not core lead-like starting points

## Practical enumeration plan

### Round 1

Enumerate around the **`CHEMBL5276142 / CHEMBL1410015 / CHEMBL1493046`** chemotype first.

Keep the core fixed and vary:

- terminal aryl electronics and substitution pattern
- cyclic amine ring size and heteroatom pattern
- polarity on the amine ring
- small linker edits that do not blow up molecular weight

### Round 2

Run a smaller focused library around:

- `CHEMBL5278336`
- `CHEMBL5410606`

### Optional side program

Only if intended, run a separate liability-aware mini-library around:

- `CHEMBL1923233`

Do not mix this one into the main lead-like enumeration funnel unless you are explicitly accepting the reactive-motif risk.

## Narrowing filters that make sense

Your rough list can be cleaned up into this practical funnel:

### Keep as real filters

- **cLogP:** target roughly `1` to `5`
- **TPSA / surface area:** prefer roughly `40` to `120` for the first pass
- **Synthetic feasibility:** use as a hard practical gate
- **PAINS / structural alert filter:** use as an exclusion or manual-review gate
- **Obvious reactivity flags:** especially Michael acceptors, strongly redox-active polyphenols, unstable acylating motifs
- **Molecular weight:** prefer `<= 500`, ideally lower for early rounds
- **Rotatable bonds:** keep controlled; Veber-style discipline is helpful

### Important, but not first-pass library filters

- **hERG risk:** important later as a liability screen, not the first rule for defining the enumeration base
- **Ames mutagenicity risk:** useful later as a liability screen
- **Kd / binding affinity:** important assay output, but not a pre-enumeration structural filter

### Not enough to define the seed by themselves

- **GenAI:** useful for proposing vectors once the seed is chosen
- **Fragmentation:** useful for ideation, but not a substitute for choosing the seed scaffold
- **SAR:** this is the analysis you do around a chosen series, not a separate filter

## Terminology cleanup

Some of the shorthand in the notes looks like this:

- `Ferber Law` -> probably **Veber rule**
- `PINS Filter` / `Pain Assay Interface` -> **PAINS filter** (`Pan-Assay Interference Compounds`)
- `HARG Binding` -> probably **hERG binding**
- `AIMS Mutation Check` -> probably **Ames mutagenicity**

## Recommendation in one sentence

Use **`CHEMBL5276142` as the first base molecule**, treat **`CHEMBL1410015` and `CHEMBL1493046` as its immediate analog reference set**, keep **`CHEMBL5278336` and `CHEMBL5410606` as secondary focused-enumeration seeds**, and keep **`CHEMBL1923233` separate unless you intentionally want to explore a reactive motif.**
