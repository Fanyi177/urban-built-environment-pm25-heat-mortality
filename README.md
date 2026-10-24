# Multidimensional built environment, PM2.5 and extreme heat-related mortality in China

## Overview

This repository contains R code associated with a nationwide ecological study of multidimensional urban built-environment (BE) characteristics, long-term PM2.5 pollution, cumulative extreme heat, and natural-cause mortality in 361 cities in mainland China from 2013 to 2020. The study drew on records of 50,891,259 natural-cause deaths.

The analytical framework considers an overall BE score and four dimension-specific scores:

- **LUUF:** land use and urban form
- **UFSD:** urban functional service density
- **TFD:** transport facility density
- **RNA:** road network accessibility

The analyses investigate two complementary processes: (1) statistical pathways involving changes in PM2.5 and cumulative extreme-heat exposure and (2) modification of the respective exposure–mortality associations. Counterfactual scenarios explore potential mortality changes under historical trends and alternative urban-planning strategies.

## Repository contents

| File | Purpose |
| --- | --- |
| `code_urban_environment_PM2.5_extreme_heat.R` | Main analysis script for built-environment characterization, mortality associations, statistical pathway decomposition, effect-modification analyses, and associated tables or figures. |
| `run_urban_resilience_scenarios_final_v4_covariates.R` | Scenario-analysis functions for projected built-environment changes, alternative urban-planning profiles, exposure prediction, mortality-impact modeling, and bootstrap uncertainty estimation. |

The repository contains analysis code, not individual-level mortality records or the complete underlying datasets. Additional preprocessing inputs, fixed PCA parameters, scenario projections, and external data may be required to reproduce the full study.

## Study design and methods

- **Study setting:** 361 cities in mainland China, 2013–2020.
- **Outcome:** annual natural-cause mortality.
- **Environmental exposures:** annual mean PM2.5 concentration and annual cumulative extreme-heat intensity.
- **Built environment:** an overall PCA-derived BE score and four dimension-specific scores (LUUF, UFSD, TFD and RNA).
- **Association analyses:** regression and generalized additive model approaches.
- **Statistical pathway decomposition:** exposure-related pathway components and remaining associations, interpreted as decompositions of observational associations rather than causal mediation effects.
- **Effect modification:** assessment of how BE characteristics modify PM2.5–mortality and extreme heat–mortality associations.
- **Planning scenarios:** historical-trend projections for 2030 and 2050 and policy-oriented profiles representing compact-city, complete-community, and active-mobility-supportive development.
- **Uncertainty:** bootstrap resampling and empirical confidence intervals, as configured in the analysis scripts.

Scenario estimates are model-based projections. They are not observed effects of implemented urban interventions. PM2.5-specific and heat-specific modeled mortality impacts may overlap and should not be summed.

## Software requirements

The code is written in **R**. It uses packages including `dplyr`, `tidyr`, `purrr`, `data.table`, `mgcv`, `randomForest`, `ggplot2`, and additional packages for spatial processing, statistical modeling, and visualization. See the beginning of each script for its package requirements.

## Running the code

1. Obtain access to the necessary input datasets through their respective data providers and any required institutional approvals.
2. Review each R script and replace computer-specific absolute file paths with paths valid for your environment. Check all expected input column names, units, and precomputed PCA or projection parameters.
3. Prepare the original city-year analysis datasets and related auxiliary inputs before running the analyses.
4. Run the main analysis script in the intended section order. The main script is not a standalone, one-command pipeline and may require files that are not distributed here.
5. Source the scenario function definitions from `run_urban_resilience_scenarios_final_v4_covariates.R`, then call the appropriate runner with your prepared data, exposure variables, BE scores, projections, and fixed PCA parameters.
6. Check scenario settings before execution. The script's configurable defaults may include years or scenario types beyond those presented in the manuscript. Align the bootstrap replicate count, covariate definitions and scenario settings with the final study specification.

**Important:** The scenario engine expects the urbanization-rate quintile to be defined consistently across the full city-year sample before stratification or bootstrap resampling, alongside the appropriate male and older-population proportions in the mortality models. Do not recompute quintile cutoffs separately within each subgroup.

## Data availability and privacy

Individual mortality records are subject to data-access restrictions and cannot be redistributed through this public repository. Other datasets may be governed by their respective providers' licenses or access conditions. Users are responsible for obtaining the appropriate approvals and source data.

No individual mortality data, restricted datasets, credentials, or nonpublic input files should be committed to this repository. Review all scripts and files before publishing, especially absolute paths, local filenames, and derived outputs that may be subject to data-use agreements.

## Reproducibility notes

This repository documents the analytical code but is not, by itself, a fully executable reproducibility package. Exact reproduction depends on the underlying datasets, harmonized variable definitions, external configuration or projection objects, and the software environment. The published manuscript and its supplementary methods should be used to identify the final analytical specifications.

## Citation

If you use this repository, please cite the associated research article once bibliographic details are available. A formal article citation and DOI can be added here after publication.

## Contact

For scientific questions, please use the repository's GitHub Issues function or the corresponding-author contact provided in the associated manuscript.
