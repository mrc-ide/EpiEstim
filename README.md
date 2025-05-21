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

### Vignettes
Please see https://mrc-ide.github.io/EpiEstim/ for vignettes with worked examples, 
FAQs and details about how EpiEstim can be used alongside some other R packages 
in an outbreak analysis workflow.

### Cite our papers
Anne Cori, Neil M. Ferguson, Christophe Fraser, Simon Cauchemez, [A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics](https://doi.org/10.1093/aje/kwt133), American Journal of Epidemiology, Volume 178, Issue 9, 1 November 2013, Pages 1505–1512. 

```console
 @article{Cori2013,
 author={Cori, A and Ferguson, NM and Fraser, C and Cauchemez, S},  
 year={2013},  
 title={{A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics}},  
 journal={Am. J. Epidemiol.},  
 doi={10.1093/aje/kwt133},  
}
```

Thompson RN, Stockwin JE, van Gaalen RD, Polonsky JA, Kamvar ZN, Demarsh PA, et al. [Improved inference of time-varying reproduction numbers during infectious disease outbreaks](https://doi.org/10.1016/j.epidem.2019.100356), Epidemics, Volume 29, 1 December 2019, 100356.

Nash RK, Nouvellet P, Cori A. [Real-time estimation of the epidemic reproduction number: Scoping review of the applications and challenges](https://doi.org/10.1371/journal.pdig.0000052), PLOS Digital Health, Volume 1, Issue 6, 27 June 2022, e0000052.

Bhatia S, Wardle J, Nash RK, Nouvellet P, Cori A. [Extending EpiEstim to estimate the transmission advantage of pathogen variants in real-time: SARS-CoV-2 as a case-study](https://doi.org/10.1016/j.epidem.2023.100692), Epidemics, 21 June 2023, 100692.

Nash RK, Cori A, Nouvellet P. [Estimating the epidemic reproduction number from temporally aggregated incidence data: a statistical modelling approach and software tool](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439), PLOS Computational Biology, Volume 19, Issue 8, 28 August 2023, e1011439.

Brizzi A, O'Driscoll M, Dorigatti I., [Refining Reproduction Number Estimates to Account for Unobserved Generations of Infection in Emerging Epidemics](https://doi.org/10.1093/cid/ciac138), Clinical Infectious Diseases, Volume 75, Issue 1, 1 July 2022, Pages e114–e121.


### Citing this code resource
We kindly request that you cite this codebase as follows (BibTeX format):

```console
@misc{Cori2022,
 author={Cori, A and Kamvar, ZN and Stockwin, J and Jombart, T and Dahlqwist, E and FitzJohn, R and Thompson, R and Nash, RK and Wardle, J and Bhatia, S},  
 year={2022},  
 title={{EpiEstim v2.2-4: A tool to estimate time varying instantaneous reproduction number during epidemics}},  
 publisher={GitHub},
 journal={GitHub repository},  
 howpublished = {\url{https://github.com/mrc-ide/EpiEstim}},  
}
```
