# EpiEstim

<!-- badges: start -->
[![R-CMD-check](https://github.com/mrc-ide/EpiEstim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mrc-ide/EpiEstim/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/mrc-ide/EpiEstim/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/EpiEstim)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3871387.svg)](https://doi.org/10.5281/zenodo.3871387)
<!-- badges: end -->

EpiEstim is a tool to estimate the time-varying instantaneous reproduction number during epidemics.
To install the latest version, use:
```r
install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev', 'https://cloud.r-project.org'))
```

## Functionalities

EpiEstim offers multiple functionalities around reproduction number estimation, 
which are summarised below. 

### Instantaneous reproduction number estimation: `estimate_R()`

The function `estimate_R` estimates the time-varying case instantaneous 
reproduction number using the methodology described in Cori et al. AJE 2013. 

Among other features, this supports the following: 

- Within the R estimation, nesting an estimation of the serial interval from 
double censored data using the `coarseDataTools` package, as described in 
Thompson et al. Epidemics 2019 
(EpiEstim versions >= 1.1-0, using `method = "si_from_data"` or 
`method = "si_from_sample"`)

- Allowing user-specified incidence split between imported and local cases, as 
described in Thompson et al. Epidemics 2019 
(v >= 1.1-0, by modifying the `incid` argument) 

- Allowing input incidence aggregated at time steps longer than 1, as described 
in Nash et al. PLoS Comp Biol 2023 (v >= 2.4, by modifying `incid` and 
specifying the `dt` argument) 

- Enabling backcalculation of early incidence (before first observation) to 
improve early R estimation as described in Brizzi et al. CID 2022. (v >= 2.5, 
with the `backimputation_window` argument)

### Estimation of the transmission advantage of pathogen variants: `estimate_advantage()`

The function `estimate_advantage` jointly estimates the instantaneous 
reproduction number for a reference variant (or pathogen or strain) and the 
relative transmissibility of a "new" variant (or pathogen or strain), 
using the methodology described in Bhatia et al. Epidemics 2023.

### Case reproduction number estimation: `wallinga_teunis()`

The function `wallinga_teunis` estimates the time-varying case reproduction 
number using the methodology described in Wallinga and Teunis et al. AJE 2004.
Note that EpiEstim currently does not implement any additional features such as
correcting for right censoring as proposed by Cauchemez et al. AJE 2006.

## Vignettes

Please see https://mrc-ide.github.io/EpiEstim/ for vignettes with worked examples, 
FAQs and details about how EpiEstim can be used alongside some other R packages 
in an outbreak analysis workflow.

## Cite our papers

The methodology underlying EpiEstim is detailed in the following papers:

Anne Cori, Neil M. Ferguson, Christophe Fraser, Simon Cauchemez, [A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics](https://doi.org/10.1093/aje/kwt133), American Journal of Epidemiology, Volume 178, Issue 9, 1 November 2013, Pages 1505–1512. 

Thompson RN, Stockwin JE, van Gaalen RD, Polonsky JA, Kamvar ZN, Demarsh PA, et al. [Improved inference of time-varying reproduction numbers during infectious disease outbreaks](https://doi.org/10.1016/j.epidem.2019.100356), Epidemics, Volume 29, 1 December 2019, 100356.

Nash RK, Nouvellet P, Cori A. [Real-time estimation of the epidemic reproduction number: Scoping review of the applications and challenges](https://doi.org/10.1371/journal.pdig.0000052), PLOS Digital Health, Volume 1, Issue 6, 27 June 2022, e0000052.

Bhatia S, Wardle J, Nash RK, Nouvellet P, Cori A. [Extending EpiEstim to estimate the transmission advantage of pathogen variants in real-time: SARS-CoV-2 as a case-study](https://doi.org/10.1016/j.epidem.2023.100692), Epidemics, 21 June 2023, 100692.

Nash RK, Cori A, Nouvellet P. [Estimating the epidemic reproduction number from temporally aggregated incidence data: a statistical modelling approach and software tool](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439), PLOS Computational Biology, Volume 19, Issue 8, 28 August 2023, e1011439.

Brizzi A, O'Driscoll M, Dorigatti I., [Refining Reproduction Number Estimates to Account for Unobserved Generations of Infection in Emerging Epidemics](https://doi.org/10.1093/cid/ciac138), Clinical Infectious Diseases, Volume 75, Issue 1, 1 July 2022, Pages e114–e121.

You can download a formatted bibtex file containing all our papers [here](inst/epiestimpapers.bib).

## Citing this code resource

We kindly request that you cite this codebase as follows (BibTeX format):

```console
  @Manual{,
    title = {EpiEstim: Estimate Time Varying Reproduction Numbers from Epidemic Curves},
    author = {Anne Cori and Rebecca Nash and Thibaut Jombart and Zhian N. Kamvar and Jake Stockwin and Robin Thompson and Sangeeta Bhatia and Andrea Brizzi and Elisabeth Dahlqwist and Rich FitzJohn and Hugo Gruson and Jack Wardle and Rolina {van Gaalen} and Tim Pollington and Pietro Monticone and Ewout {ter Hoeven} and Jonathan A. Polonsky and Shikun Li and Justin Lessler and P. Alex Demarsh and Neil M. Ferguson and Christophe Fraser and Simon Cauchemez},
    year = {2026},
    note = {R package version 2.4},
    url = {https://github.com/mrc-ide/EpiEstim},
  }

```

You can also retrieve the above entry using 

```
cite("EpiEstim")
```
