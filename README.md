# Model-Based Geostatistics

`mbg` is an R package for modern model-based geostatistics. It combines `sf`, `terra`, and
`data.table` functions for spatial data processing, `caret` for spatial ML models, and
`R-INLA` for geostatistical models.

## Partner packages

The `mbg` packages works well with two other partner R packages, although they are not
dependencies:

- [`versioning`](https://github.com/henryspatialanalysis/versioning): Work with config files and versioned input/output folders
- [`pixel2poly`](https://github.com/henryspatialanalysis/pixel2poly): Translate between raster surfaces and administrative unit boundaries

## Package workflow

The general MBG workflow using this package is as follows:

- Load point data on outcomes, raster covariate surfaces, and a raster population surface
- (Optional): run machine learning models relating the input covariate surfaces to the outcome, producing predictive raster surfaces from a variety of methods
- Prepare geostatistical model inputs. This includes the outcomes point data, model specifications, a spatial 2-D mesh, and either the input covariate surfaces or the ML predictive surfaces
- Run the geostatistical model using `R-INLA`. This model predicts the outcome as a linear combination of the raster surfaces and a SPDE approximation to a Gaussian process over space.
- Using the model fit, predict results across the entire study area. Uncertainty is captured by generating 250 "draws" at each pixel location.
- Summarize draws as raster surfaces by taking the mean, median, and 95% uncertainty interval bounds of draws at each pixel location
- (Optional): aggregate from pixels to administrative boundaries, preserving uncertainty

## Using the package

You can install the R package from Github:

```devtools::install_github("henryspatialanalysis/mbg")```

Or clone the package to your computer and load it dynamically (this is the preferred option while the package is under active development):

```devtools::load_all("/local/path/to/mbg/")```

After loading the package, you can access package documentation by running `help(mbg)`, or get documentation for a specific function by running e.g. `help(fit_inla_model)`.
