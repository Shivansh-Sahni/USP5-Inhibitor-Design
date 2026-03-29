# Similarity-Window Validation

## Result

A new similarity-constrained enumeration library was generated in:

- `outputs/enumeration_similarity_window_0p55_0p95.csv`

This library contains:

- `2255` unique products
- mean nearest-original Morgan Tanimoto similarity: `0.6827`
- median nearest-original Morgan Tanimoto similarity: `0.6724`
- enforced similarity window: `0.55` to `0.95`

This is a much better fit for controlled analog exploration than the earlier broad library, because it stays in a moderate-to-high similarity band relative to the 13 original positive parents.

## Why `0.55-0.95` is defensible

There is **no universal Tanimoto cutoff** that is correct for every fingerprint and every task. However, the literature supports the idea that:

- high similarity cutoffs around `0.85` are commonly used as a marker of close analogs
- meaningful activity-relevant thresholds for circular fingerprints can be much lower than `0.85`
- clustering and neighborhood analyses often use thresholds in the rough `0.5-0.7+` range for useful structure grouping with ECFP-like fingerprints

That makes `0.55-0.95` a reasonable project-specific window for:

- staying meaningfully related to known actives
- allowing controlled deviation
- avoiding exact duplicates and trivial near-identity molecules at the top end

## Articles supporting this choice

### 1. Jasial et al., 2016

Article:

- [Activity-relevant similarity values for fingerprints and implications for similarity searching](https://pmc.ncbi.nlm.nih.gov/articles/PMC4830209/)

Why it matters:

- The paper explicitly says there is **no generally applicable Tc threshold** for activity similarity.
- It also shows that activity-relevant similarity values depend strongly on the fingerprint.
- For ECFP4-like fingerprints, useful activity-relevant values can be substantially lower than the classic `0.85` rule-of-thumb used for older fingerprints.

Why this supports the project window:

- It argues against insisting on `0.85+` as the only “valid” region.
- It supports using a lower similarity floor for circular fingerprints when controlled exploration is desired.

### 2. Zahoránszky-Kőhalmi et al., 2016

Article:

- [Impact of similarity threshold on the topology of molecular similarity networks and clustering outcomes](https://link.springer.com/article/10.1186/s13321-016-0127-5)

Why it matters:

- This paper studies how different Tanimoto thresholds affect molecular similarity networks and clustering outcomes.
- For ECFP4/Tanimoto analyses, useful threshold behavior appears in the moderate range, and example clustering maxima are reported around `0.68` to `0.72` in large datasets.

Why this supports the project window:

- It supports the idea that the `0.55-0.70+` region is not arbitrary noise; it is still chemically meaningful for neighborhood and cluster structure.
- It is consistent with treating `0.55-0.70` as controlled exploration and `0.70+` as close analog space.

### 3. TCMSID database paper, 2022

Article:

- [TCMSID: a simplified integrated database for drug discovery from traditional chinese medicine](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00670-z)

Why it matters:

- The authors explicitly use `Tc = 0.85` as a threshold for **high** similarity and `Tc = 0.5` as a threshold for **medium** similarity.

Why this supports the project window:

- It gives a published example where the literature-style interpretation is:
  - `0.85` high similarity
  - `0.5` medium similarity
- Our chosen band of `0.55-0.95` sits squarely between medium similarity and very high similarity.

### 4. Martin et al., 2002 / related medicinal-chemistry practice

Article:

- [2D similarity searching: how much 3D similarity is enough?](https://pubs.acs.org/doi/abs/10.1021/jm020155c)

Why it matters:

- This work is one of the sources often cited for the practical use of `0.85` as a strong similarity cutoff in medicinal chemistry and library design contexts.

Why this supports the project window:

- It reinforces using the upper part of the window, especially `0.70-0.95`, as near-analog territory rather than random drift.

## Practical interpretation for this project

For Morgan fingerprint similarity to the 13 original positive parents:

- `0.85-0.95`: very close analogs
- `0.70-0.85`: close analogs suitable for focused SAR
- `0.55-0.70`: controlled exploration
- `< 0.55`: broader chemical exploration or scaffold drift

Accordingly, this new library should be treated as:

- a **focused-to-moderately exploratory analog library**
- not a pure scaffold-hopping library
- much better aligned with lead-finding and early SAR than the earlier unconstrained `3264`-compound broad library
