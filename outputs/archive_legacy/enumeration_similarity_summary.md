# Enumeration Similarity To Original Positives

Similarity metric: Morgan fingerprint Tanimoto, radius 2, 2048 bits.
Reference set: 13 original positive molecules from the training/modeling dataset.

## Overall unique-product similarity

- Unique products: 3264
- Mean max similarity to any original: 0.3085
- Median max similarity to any original: 0.2818
- Min / max similarity: 0.1429 / 1.0
- Very high (>= 0.85): 36
- High (0.70-0.85): 82
- Medium (0.55-0.70): 30
- Low (< 0.55): 3116

## Per-method mean max similarity

- `ring_size_and_ring_system_scans`: mean 0.9333, median 1.0, min 0.8, max 1.0
- `bioisosteric_swaps`: mean 0.854, median 0.8571, min 0.7778, max 0.9082
- `matched_molecular_pair_expansion`: mean 0.8207, median 0.825, min 0.6316, max 1.0
- `reaction_based_enumeration`: mean 0.8032, median 0.807, min 0.6, max 0.939
- `heteroatom_walks`: mean 0.7808, median 0.7857, min 0.6122, max 0.907
- `linker_scans`: mean 0.7338, median 0.7273, min 0.6875, max 0.7846
- `r_group_enumeration`: mean 0.3038, median 0.2921, min 0.1429, max 0.8261
- `scaffold_hopping`: mean 0.28, median 0.2576, min 0.1489, max 0.6667
- `reagent_pool_combinatorics`: mean 0.2762, median 0.2772, min 0.18, max 0.4561
- `fragment_growing`: mean 0.2682, median 0.2419, min 0.1489, max 0.6667