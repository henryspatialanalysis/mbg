# Run variance inflation factor (VIF) selection on input covariates

Run variance inflation factor (VIF) selection on input covariates

## Usage

``` r
vif_covariate_select(dataset, vif_cutoff = 5)
```

## Arguments

- dataset:

  data.frame-like object with named columns containing all covariates to
  consider in the VIF analysis.

- vif_cutoff:

  (`numeric(1)`) Cutoff for maximum variance inflation factor in
  `dataset`

## Value

data.table listing each variable, VIF in most recent round, and whether
the indicator should be included or not
