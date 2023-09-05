## #######################################################################################
##
## Example SPDE model with covariates in INLA
##
## Gambia example based on Paula Moraga's Geospatial Health Data book:
## https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldataexamplespatial.html
##
## Adapted to mainly use the `terra` and `sf` packages, plus full INLA predictive draws
##
## #######################################################################################

load_packages <- c(
  'assertthat', 'data.table', 'geodata', 'ggplot2', 'INLA', 'rgeoboundaries', 'sf', 'terra'
)
invisible(lapply(load_packages, library, character.only = T))
# Use pixel2poly to summarize posterior draws
devtools::load_all('~/repos/pixel2poly')

options(geodata_default_path = '~/temp_data/geodata')

# Save results and visualizations to a folder
save_dir <- '~/temp_data/geostats/gambia_example'


## Prepare input data: country boundary, data, and covariate ---------------------------->

# Load country boundaries
ad0_sf <- rgeoboundaries::geoboundaries('Gambia', adm_lvl = 0)
ad1_sf <- rgeoboundaries::geoboundaries('Gambia', adm_lvl = 1)
ad2_sf <- rgeoboundaries::geoboundaries('Gambia', adm_lvl = 2)
ad2_sf$adm2_id <- seq_len(nrow(ad2_sf))

# Load input data: malaria cases
individual_data <- get(data(gambia, package = 'geoR')) |> data.table::as.data.table()
group_data <- (
  individual_data
  [, .(pos = sum(pos), samplesize = .N), by = .(x, y)]
  [, prevalence := pos / samplesize]
)
# Convert from UTM zone 28 to lat-long
latlong_dt <- (
  sf::st_as_sf(
    group_data,
    coords = c('x', 'y'),
    crs = sf::st_crs('+proj=utm +zone=28')
  ) |>
  sf::st_transform(crs = sf::st_crs(4326)) |>
  sf::st_coordinates() |>
  as.data.table()
)
group_data[, c('x', 'y') := NULL ]
group_data$x <- latlong_dt$X
group_data$y <- latlong_dt$Y

# Get mean annual precipitation by grid cell across the gambia
covar_raw <- geodata::worldclim_country(country = 'GMB', var = 'prec') |>
  terra::app(sum)
# Aggregate to a more coarse resolution
covar <- terra::aggregate(covar_raw, fact = 5, fun = "mean", na.rm=T) |>
  terra::extend(y = ad2_sf, snap = 'out')
# Normalize
covar <- (
  (covar - terra::global(covar, mean, na.rm=T)[1, 1]) /
    terra::global(covar, sd, na.rm=T)[1, 1]
)
group_data$prec <- terra::extract(
  x = covar,
  y = terra::vect(group_data, geom = c('x', 'y')),
  ID = F
)[, 1]
group_data$intercept <- 1

# The ID raster (prediction locations) has the same extent as the covariate
id_raster <- pixel2poly::build_id_raster(
  polygons = terra::vect(ad2_sf),
  template_raster = covar
)
covar <- terra::crop(x = covar, y = id_raster)
# Create an aggregation table
agg_table <- pixel2poly::build_aggregation_table(
  polygons = terra::vect(ad2_sf),
  id_raster = id_raster,
  polygon_id_field = 'adm2_id'
)
# Turn into a table of prediction locations
prediction_locations <- data.table::as.data.table(id_raster, xy = T) |> na.omit()
names(prediction_locations)[3] <- 'pixel_id'
prediction_locations$intercept <- 1
prediction_locations$prec <- terra::extract(
  covar,
  y = terra::vect(prediction_locations, geom = c('x', 'y')),
  ID = FALSE
)


## Build INLA data stack ---------------------------------------------------------------->

# Prediction mesh
mesh <- INLA::inla.mesh.2d(
  loc.domain = prediction_locations[, .(x, y)],
  max.edge = c(0.1, 5),
  cutoff = 0.01
)

# SPDE object
spde <- INLA::inla.spde2.matern(
  mesh = mesh,
  alpha = 2,
  constr = TRUE
)

# Spatial index
index_s <- INLA::inla.spde.make.index(name = "s", n.spde = spde$n.spde)

# Projection matrix: mesh to data
A_proj_data <- INLA::inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(group_data[, .(x, y)])
)

# Projection matrix: mesh to all prediction locations
A_proj_predictions <- INLA::inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(prediction_locations[, .(x, y)])
)

# Data stack for estimation
stack_estimation <- INLA::inla.stack(
  tag = 'est',
  data = list(y = group_data$pos, samplesize = group_data$samplesize),
  A = list(1, A_proj_data),
  effects = list(group_data[, .(intercept, prec)], s = index_s)
)

# INLA model formula
formula_string <- 'y ~ 0 + intercept + prec + f(s, model = spde)'


## Run INLA model ----------------------------------------------------------------------->

fitted_model <- INLA::inla(
  stats::as.formula(formula_string),
  family = 'binomial',
  Ntrials = samplesize,
  control.family = list(link = 'logit'),
  data = inla.stack.data(stack_estimation),
  control.predictor = list(
    compute = FALSE,
    link = 1,
    A = inla.stack.A(stack_estimation)
  ),
  control.compute = list(
    config = TRUE
  )
)

## Generate predictive draws across the entire country ---------------------------------->

posterior_samples <- INLA::inla.posterior.sample(
  n = 250,
  result = fitted_model,
  add.names = FALSE
)

fes <- c('intercept', 'prec')
spatial_res <- c('s')

latent_matrix <- lapply(posterior_samples, function(x) x$latent) |> do.call(what = cbind)
rownames(latent_matrix) <- vapply(
  strsplit(rownames(posterior_samples[[1]]$latent), split = ':'),
  `[`, 1, FUN.VALUE = character(1)
)

# Split into fixed effect coeffients and spatial mesh effects
fe_rows <- sapply(fes, function(fe){
  which_row <- which(rownames(latent_matrix) == fe)
  # There should be exactly one row per covariate effect
  assertthat::assert_that(length(which_row) == 1)
  return(which_row)
})
fe_coefficients <- latent_matrix[fe_rows, ]

spatial_mesh_effects <- latent_matrix[rownames(latent_matrix) == 's', ]
assertthat::assert_that(nrow(spatial_mesh_effects) == ncol(A_proj_predictions))

# Project to all grid locations
fe_draws <- as.matrix(prediction_locations[, .(intercept, prec)]) %*% fe_coefficients
re_draws <- as.matrix(A_proj_predictions %*% spatial_mesh_effects)

assertthat::assert_that(all.equal(dim(fe_draws), dim(re_draws)))
predictive_draws <- plogis(fe_draws + re_draws)

# Insert into raster
mean_raster <- pixel2poly::values_to_raster(
  Matrix::rowMeans(predictive_draws),
  id_raster = id_raster
)
lower_raster <- pixel2poly::values_to_raster(
  matrixStats::rowQuantiles(predictive_draws, probs = 0.025),
  id_raster = id_raster
)
upper_raster <- pixel2poly::values_to_raster(
  matrixStats::rowQuantiles(predictive_draws, probs = 0.975),
  id_raster = id_raster
)

# Aggregate pixel draws to admin draws
# NOT population weighted for example - normally this would be pop-weighted
ad2_draws <- pixel2poly::aggregate_draws_to_polygons(
  draws_matrix = predictive_draws,
  aggregation_table = agg_table,
  aggregation_cols = c('polygon_id', 'shapeName'),
  method = 'mean',
  na.rm = T
)
draw_cols <- grep('^draw_', colnames(ad2_draws), value = T)
draws_mat <- as.matrix(ad2_draws[, ..draw_cols])
ad2_summary <- data.table::data.table(
  polygon_id = ad2_draws$polygon_id,
  shapeName = ad2_draws$shapeName,
  mean = rowMeans(draws_mat),
  lower = matrixStats::rowQuantiles(draws_mat, probs = 0.025),
  upper = matrixStats::rowQuantiles(draws_mat, probs = 0.975)
)

## Save and plot some basic results ----------------------------------------------------->

data.table::fwrite(predictive_draws, file = file.path(save_dir, 'pred_draws.csv'))
data.table::fwrite(ad2_draws, file = file.path(save_dir, 'pred_draws_ad2.csv'))

terra::writeRaster(x = mean_raster, file = file.path(save_dir, 'pred_mean_raster.tif'), overwrite = T)
terra::writeRaster(x = lower_raster, file = file.path(save_dir, 'pred_lower_raster.tif'), overwrite = T)
terra::writeRaster(x = upper_raster, file = file.path(save_dir, 'pred_upper_raster.tif'), overwrite = T)
data.table::fwrite(ad2_summary, file = file.path(save_dir, 'pred_summary_ad2.csv'))

# Plot the mesh outlines
pdf(file.path(save_dir, 'mesh_s.pdf'), height = 6, width = 6)
plot(mesh)
lines(ad0_sf, col = 'red')
dev.off()

# Plot the summary rasters
rplot <- function(rr, title){
  rr_dt <- as.data.table(rr, xy = T)
  colnames(rr_dt)[3] <- 'val'
  fig <- ggplot() +
    geom_raster(data = rr_dt, aes(x=x, y=y, fill = val)) +
    geom_sf(data = ad0_sf, fill = NA, color = '#222222') +
    labs(title = title, fill = 'Estimate') +
    scale_fill_gradientn(
      colors = viridisLite::viridis(100), limits = c(0, 1), labels = scales::percent
    ) +
    theme_minimal()
  return(fig)
}
pdf(file.path(save_dir, 'raster_summaries.pdf'), height = 6, width = 6)
print(rplot(mean_raster, title = 'Mean pixel estimates'))
print(rplot(lower_raster, title = 'Lower (2.5th %ile) pixel estimates'))
print(rplot(upper_raster, title = 'Upper (97.5th %ile) pixel estimates'))
dev.off()

# Plot the admin summaries
merge_vars <- c('polygon_id', 'shapeName')
ad2_sf$polygon_id <- ad2_sf$adm2_id 
ad2_summary_sf <- merge(
  x = ad2_sf,
  y = melt(ad2_summary, id.vars = merge_vars),
  by = merge_vars,
  allow.cartesian = T
)
poly_fig <- ggplot() + 
  facet_wrap('variable', ncol = 1) +
  geom_sf(data = ad2_summary_sf, aes(fill = value), linewidth = .1) +
  geom_sf(data = ad0_sf, fill = NA, color = '#222222', linewidth = .25) +
  labs(title = 'Admin2 summaries', fill = 'Prevalence') +
  scale_fill_gradientn(
    colors = viridisLite::viridis(100), limits = c(0, 1), labels = scales::percent
  ) +
  theme_minimal()
pdf(file.path(save_dir, 'ad2_summaries.pdf'), height = 10, width = 6)
print(poly_fig)
dev.off()
