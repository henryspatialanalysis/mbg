# Model-Based Geostatistics

`mbg` is an R package for modern model-based geostatistics. `mbg` provides a simple interface to run spatial machine learning models and geostatistical models that estimate a continuous (raster) surface from point-referenced outcomes and, optionally, a set of raster covariates. The package also includes functions to summarize raster outcomes by (polygon) region while preserving uncertainty.

The `mbg` package combines features from the [`sf`](https://r-spatial.github.io/sf/), [`terra`](https://rspatial.github.io/terra/), and [`data.table`](https://CRAN.R-project.org/package=data.table) packages for spatial data processing; [`caret`](https://topepo.github.io/caret/) for spatial ML models; and [`R-INLA`](https://www.r-inla.org/) for geostatistical models.

---

## Using the package

You can install the R package from Github:

```devtools::install_github("henryspatialanalysis/mbg")```

Or clone the package to your computer and load it dynamically:

```devtools::load_all("/local/path/to/mbg/")```

After loading the package, you can access package documentation by running `help(mbg)`, or get documentation for a specific function by running e.g. `help(MbgModelRunner)`.

---

## Package workflow

The general MBG workflow using this package is as follows:

1. Load point data on outcomes, raster covariate surfaces, and a raster population surface
2. _(Optional):_ Run machine learning models relating the input covariate surfaces to the outcome, producing predictive raster surfaces from a variety of methods
3. Prepare geostatistical model inputs. This includes the outcomes point data, model specifications, a spatial 2-D mesh, and either the input covariate surfaces or the ML predictive surfaces
4. Run the geostatistical model. This model predicts the outcome as a linear combination of the raster surfaces and a SPDE approximation to a Gaussian process over space.
5. Using the model fit, predict results across the entire study area. Uncertainty is captured by generating 250 posterior predictive draws at each pixel location.
6. Summarize draws as raster surfaces by taking the mean, median, and 95% uncertainty interval bounds of draws at each pixel location
7. _(Optional):_ Aggregate from pixels to administrative boundaries, preserving uncertainty

For example workflows, see the [introductory vignette](https://henryspatialanalysis.github.io/mbg/articles/mbg.html).
