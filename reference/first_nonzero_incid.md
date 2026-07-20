# First day of non-zero incidence

Get the first day of non-zero incidence across all variants and
locations.

## Usage

``` r
first_nonzero_incid(incid)
```

## Arguments

- incid:

  a multidimensional array containing values of the incidence for each
  time step (1st dimension), location (2nd dimension) and
  pathogen/strain/variant (3rd dimension)

## Value

integer

## Details

For each variant, find the first day of non-zero incidence. The maximum
of these is the smallest possible point at which estimation can begin.

## Author

Sangeeta Bhatia
