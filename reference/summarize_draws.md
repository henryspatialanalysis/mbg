# Summarize draws

Helper function to summarize a matrix or data.frame of predictive draws

## Usage

``` r
summarize_draws(
  draws,
  id_fields = NULL,
  draw_fields = NULL,
  ui_width = 0.95,
  na.rm = TRUE
)
```

## Arguments

- draws:

  A `matrix`, `data.frame`, or
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  of predictive draws.

- id_fields:

  (default NULL) Only considered for data.frame-like `draws`. What
  identifier fields in the data should be kept in the summary table and
  not included among the draw fields?

- draw_fields:

  (default NULL) Only considered for data.frame-like `draws`. What
  fields represent actual draws, as opposed to identifier fields or
  other metadata like population? If `NULL`, the default, automatically
  determines the draw fields as all columns not included in the
  `id_fields`.

- ui_width:

  (`numeric`, default 0.95) Size of the uncertainty interval width when
  calculating the upper and lower summary rasters

- na.rm:

  (`logical`, default TRUE) Should NA values be removed when calculating
  summaries across draws?

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
containing at least the following fields:

- The `id_fields`, if passed

- "mean": Mean across predictive draws

- "lower": Lower bound of the (X%) uncertainty interval

- "upper": Upper bound of the (X%) uncertainty interval

- "ui_width": "upper" - "lower"

## Examples

``` r
# Summarize a draws matrix
draws_matrix <- matrix(rnorm(200), nrow = 10)
summary_table_a <- summarize_draws(draws_matrix)
head(summary_table_a)
#>           mean     lower    upper ui_width
#>          <num>     <num>    <num>    <num>
#> 1:  0.04883910 -1.696132 1.984732 3.680864
#> 2: -0.03271432 -1.279289 1.186630 2.465919
#> 3:  0.12182437 -2.075344 1.975478 4.050822
#> 4:  0.06324496 -1.476200 1.577184 3.053384
#> 5:  0.24955491 -1.774292 1.829265 3.603557
#> 6:  0.54884338 -1.441610 2.456848 3.898458

# Summarize a draws data.table with location IDs
draws_table <- matrix(c(1:10, rnorm(200)), nrow = 10) |>
  data.table::as.data.table() |>
  data.table::setnames(c('location_id', paste0('draw_', 1:20)))
summary_table_b <- summarize_draws(draws_table, id_fields = 'location_id')
head(summary_table_b)
#>    location_id        mean     lower    upper ui_width
#>          <num>       <num>     <num>    <num>    <num>
#> 1:           1 -0.13098397 -1.034371 1.033082 2.067454
#> 2:           2 -0.26414518 -2.083516 1.741731 3.825247
#> 3:           3  0.35835144 -1.106188 1.667994 2.774182
#> 4:           4 -0.10309306 -1.697354 1.770278 3.467632
#> 5:           5  0.07679248 -1.964361 1.878477 3.842837
#> 6:           6  0.37618266 -1.823267 1.959878 3.783145
```
