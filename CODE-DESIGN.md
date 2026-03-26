# Code Design

EpiEstim is an R package for estimating the time-varying instantaneous reproduction number from epidemic curves. This document describes the design and architecture of the codebase.

## Core Workflow

The typical user workflow is:

1. Get Incidence data (e.g., incidence & incidence2)
2. Estimate R - uses configuration (`make_config()`)
    - `estimate_R()` (also wrappers `estimate_R_agg()`, `estimate_advantage()`
    - `wallinga_teunis()`
3. Plot results (`plot()`)
4. Send to other packages (e.g., projections)

## Be Careful

- [Any common pitfalls a developer working with this package might encounter?]

## Style

- Use R version >= 3.3.0 (i.e., no `|>` and no `\(x)`)
- No use of pipes base or magrittr in code base (can be used in vignettes)
- pkgdown for documentation website
- ggplot2 for plots and patchwork for combining plots
- Validate inputs at the start of functions
- Use base messages (i.e., `stop()`, `warning()`, `message()`)
- Prefer `pkg::foo()` usage of dependencies, avoid full package imports
- Write documentation with roxygen2 and markdown (https://roxygen2.r-lib.org/articles/rd-formatting.html)
- Re-use documentation with `@inheritXXX` (https://roxygen2.r-lib.org/articles/reuse.html)


## Testing

- Tests use testthat (not v3)
- Plot tests use [`vdiffr`](https://vdiffr.r-lib.org/) for visual snapshot tests.


## Data

Built-in datasets include:

- **Flu2009**: 1918 pandemic influenza data from England with serial interval distribution
- **SARS2003**: SARS outbreak data from Hong Kong
- **Measles1861**: Historical measles outbreak in England
- **Smallpox1972**: Last smallpox case cluster in England
- **Flu1918**: 1918 influenza pandemic data
- **MockRotavirus**: Simulated rotavirus outbreak data
- **MERS2014_15**: MERS outbreak data from 2014-2015
- **covid_deaths_2020_uk**: COVID-19 deaths in UK (2020)
- **flu_2009_NYC_school**: 2009 influenza in NYC school

Accessed via `data("DatasetName")` after loading package.
