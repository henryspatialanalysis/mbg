# Dissolve sf object by attribute

Dissolve an SF object by attribute

## Usage

``` r
dissolve_sf_by_attribute(x, by = character(0))
```

## Arguments

- x:

  ([sf::sf](https://r-spatial.github.io/sf/reference/sf.html) object) SF
  object to dissolve

- by:

  (`character(N)`, default character(0)) Attributes to dissolve by

## Value

Dissolved [sf::sf](https://r-spatial.github.io/sf/reference/sf.html)
object

## Details

Inspired by spatialEco::sf_dissolve

## Examples

``` r
if (FALSE) { # \dontrun{
  communes_sf <- sf::st_read(system.file("extdata/Benin_communes.gpkg", package = "mbg"))
  departments_sf <- mbg::dissolve_sf_by_attribute(
    x = communes_sf,
    by = c('department', 'department_code')
  )
} # }
```
