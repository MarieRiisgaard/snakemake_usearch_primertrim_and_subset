# Changelog

## [1.2.1](https://github.com/MarieRiisgaard/snakemake_usearch_primertrim_and_subset/compare/v1.2.0...v1.2.1) (2025-10-28)


### Bug Fixes

* corrected read summary parsing (NR==1), increased memory for concat_subsample, and improved workflow stability ([6ac388c](https://github.com/MarieRiisgaard/snakemake_usearch_primertrim_and_subset/commit/6ac388ca4e4581342c4bada6077969be2d5ca366))
* remove zero-length reads after primer trimming ([ddb4bd9](https://github.com/MarieRiisgaard/snakemake_usearch_primertrim_and_subset/commit/ddb4bd97f025380338bd4f73059df7e682f624fe))

## [1.1.1](https://github.com/KasperSkytte/snakemake_usearch/compare/v1.1.0...v1.1.1) (2025-09-03)


### Bug Fixes

* ignore "unclassified" folder correctly+add description of steps in readme ([8448c90](https://github.com/KasperSkytte/snakemake_usearch/commit/8448c90f68ce314698d27eb7481f50e53136fa35))
* typo ([75bf7e0](https://github.com/KasperSkytte/snakemake_usearch/commit/75bf7e0c4e6261b6bdfa1e07775c373f7ba0391d))

## [1.1.0](https://github.com/KasperSkytte/snakemake_usearch/compare/v1.0.0...v1.1.0) (2025-08-28)


### Features

* make otutab rarefaction optional (requires usearch 11) ([df737b9](https://github.com/KasperSkytte/snakemake_usearch/commit/df737b9f2daf5ec1ea4fbf2e190e84fd5485a042))


### Bug Fixes

* Docker container doesn't load conda env when run using apptainer ([a91f2fe](https://github.com/KasperSkytte/snakemake_usearch/commit/a91f2fe2685ce7d53d31ba145003666b08dd2c6b))
* env was unsolvable with strict channel prio ([859af74](https://github.com/KasperSkytte/snakemake_usearch/commit/859af74578ac48810be1707fbd6c3f44b9e01a50))
* exit a task if output file is empty ([b6548e8](https://github.com/KasperSkytte/snakemake_usearch/commit/b6548e8f44e0cc950203d14db6addc496605082a))
* update conda env hash ([8259f64](https://github.com/KasperSkytte/snakemake_usearch/commit/8259f649169095311f3f0a1ebcee79d2292c0dc4))

## 1.0.0 (2025-08-21)


### Bug Fixes

* lint action failed ([9824c46](https://github.com/KasperSkytte/snakemake_usearch/commit/9824c46691e46385890ef6f33cf0463a548d84dc))
