## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

download_data <- FALSE
create_custom_data <- FALSE

## ----initialis_pastclim, echo=FALSE, results="hide", eval=!download_data------
library(pastclim)
set_data_path(on_CRAN = TRUE)

## -----------------------------------------------------------------------------
library(pastclim)
tavg_vars <- c(paste0("temperature_0", 1:9), paste0("temperature_", 10:12))
time_steps <- get_time_bp_steps(dataset = "Example")
n_europe_ext <- c(-10, 15, 45, 60)

## ----eval=download_data-------------------------------------------------------
# download_dataset(dataset = "Beyer2020", bio_variables = tavg_vars)
# tavg_series <- region_series(
#   bio_variables = tavg_vars,
#   time_bp = time_steps,
#   dataset = "Beyer2020",
#   ext = n_europe_ext
# )

## ----echo=FALSE, results="hide", eval=create_custom_data----------------------
# terra::saveRDS(tavg_series,
#   file = "../inst/extdata/delta/tavg_series.RDS"
# )

## ----echo=FALSE, results="hide", eval=!download_data--------------------------
library(pastclim)
set_data_path(on_CRAN = TRUE)
tavg_series <- terra::readRDS(system.file("extdata/delta/tavg_series.RDS",
  package = "pastclim"
))
# get back the time units that are lost when saving the rds
old_names <- names(tavg_series) # there is a bug in terra
terra::time(tavg_series, tstep = "years") <- terra::time(tavg_series)
names(tavg_series) <- old_names
rm(old_names)

## -----------------------------------------------------------------------------
tavg_model_lres_rast <- tavg_series$temperature_01
tavg_model_lres_rast

## ----fig.width=6, fig.height=5------------------------------------------------
plot(tavg_model_lres_rast, main = time_bp(tavg_model_lres_rast))

## ----eval=download_data-------------------------------------------------------
# download_dataset(dataset = "WorldClim_2.1_10m", bio_variables = tavg_vars)
# tavg_obs_hres_all <- region_series(
#   bio_variables = tavg_vars,
#   time_ce = 1985,
#   dataset = "WorldClim_2.1_10m",
#   ext = n_europe_ext
# )

## ----echo=FALSE, results="hide", eval=create_custom_data----------------------
# terra::saveRDS(tavg_obs_hres_all,
#   file = "../inst/extdata/delta/tavg_obs_hres_all.RDS"
# )

## ----echo=FALSE, results="hide", eval=!download_data--------------------------
tavg_obs_hres_all <- terra::readRDS(
  system.file("extdata/delta/tavg_obs_hres_all.RDS",
    package = "pastclim"
  )
)

## -----------------------------------------------------------------------------
tavg_obs_range <- range(
  unlist(
    lapply(tavg_obs_hres_all, minmax, compute = TRUE)
  )
)
tavg_obs_range

## ----fig.width=4, fig.height=4------------------------------------------------
tavg_obs_hres_all <- terra::crop(tavg_obs_hres_all, n_europe_ext)
# extract the January raster
tavg_obs_hres_rast <- tavg_obs_hres_all[[1]]
plot(tavg_obs_hres_rast)

## -----------------------------------------------------------------------------
ext(tavg_obs_hres_rast) == ext(tavg_model_lres_rast)

## ----eval=download_data-------------------------------------------------------
# download_etopo()
# relief_rast <- load_etopo()
# relief_rast <- terra::resample(relief_rast, tavg_obs_hres_rast)

## ----echo=FALSE, results="hide", eval=create_custom_data----------------------
# terra::saveRDS(relief_rast,
#   file = "../inst/extdata/delta/relief_rast.RDS"
# )

## ----echo=FALSE, results="hide", eval=!download_data--------------------------
relief_rast <- terra::readRDS(system.file("extdata/delta/relief_rast.RDS",
  package = "pastclim"
))

## ----fig.width=6, fig.height=5------------------------------------------------
land_mask_high_res <- make_land_mask(
  relief_rast = relief_rast,
  time_bp = time_bp(tavg_model_lres_rast)
)
plot(land_mask_high_res, main = time_bp(land_mask_high_res))

## ----eval=download_data-------------------------------------------------------
# ice_mask_low_res <- get_ice_mask(time_bp = time_steps, dataset = "Beyer2020")
# ice_mask_high_res <- downscale_ice_mask(
#   ice_mask_low_res = ice_mask_low_res,
#   land_mask_high_res = land_mask_high_res
# )
# plot(ice_mask_high_res)

## ----echo=FALSE, results="hide", eval=create_custom_data----------------------
# terra::saveRDS(ice_mask_low_res,
#   file = "../inst/extdata/delta/ice_mask_low_res.RDS"
# )

## ----echo=FALSE, eval=!download_data------------------------------------------
ice_mask_low_res <- terra::readRDS(
  system.file("extdata/delta/ice_mask_low_res.RDS",
    package = "pastclim"
  )
)
ice_mask_high_res <- downscale_ice_mask(
  ice_mask_low_res = ice_mask_low_res,
  land_mask_high_res = land_mask_high_res
)
plot(ice_mask_high_res)

## -----------------------------------------------------------------------------
land_mask_high_res <- mask(land_mask_high_res,
  ice_mask_high_res,
  inverse = TRUE
)
plot(land_mask_high_res)

## ----eval=FALSE---------------------------------------------------------------
# internal_seas <- readRDS(system.file("extdata/internal_seas.RDS",
#   package = "pastclim"
# ))
# land_mask_high <- mask(land_mask_high_res,
#   internal_seas,
#   inverse = TRUE
# )

## -----------------------------------------------------------------------------
delta_rast <- delta_compute(
  x = tavg_model_lres_rast, ref_time = 0,
  obs = tavg_obs_hres_rast
)
model_downscaled <- delta_downscale(
  x = tavg_model_lres_rast,
  delta_rast = delta_rast,
  x_landmask_high = land_mask_high_res,
  range_limits = tavg_obs_range
)
model_downscaled

## ----fig.width=6, fig.height=5------------------------------------------------
panel(model_downscaled, main = time_bp(model_downscaled))

## ----fig.width=6, fig.height=5------------------------------------------------
panel(tavg_model_lres_rast, main = time_bp(tavg_model_lres_rast))

## -----------------------------------------------------------------------------
tavg_downscaled_list <- list()
for (i in 1:12) {
  delta_rast <- delta_compute(
    x = tavg_series[[i]], ref_time = 0,
    obs = tavg_obs_hres_all[[i]]
  )
  tavg_downscaled_list[[i]] <- delta_downscale(
    x = tavg_series[[i]],
    delta_rast = delta_rast,
    x_landmask_high = land_mask_high_res,
    range_limits = tavg_obs_range
  )
}
tavg_downscaled <- terra::sds(tavg_downscaled_list)

## -----------------------------------------------------------------------------
tavg_downscaled

## ----eval=download_data-------------------------------------------------------
# prec_vars <- c(paste0("precipitation_0", 1:9), paste0("precipitation_", 10:12))
# prec_series <- region_series(
#   bio_variables = prec_vars,
#   time_bp = time_steps,
#   dataset = "Beyer2020",
#   ext = n_europe_ext
# )

## ----echo=FALSE, results="hide", eval=create_custom_data----------------------
# terra::saveRDS(prec_series,
#   file = "../inst/extdata/delta/prec_series.RDS"
# )

## ----echo=FALSE, results="hide", eval=!download_data--------------------------
prec_vars <- c(paste0("precipitation_0", 1:9), paste0("precipitation_", 10:12))
prec_series <- terra::readRDS(system.file("extdata/delta/prec_series.RDS",
  package = "pastclim"
))
# get back the time units that are lost when saving the rds
old_names <- names(prec_series) # there is a bug in terra
terra::time(prec_series, tstep = "years") <- terra::time(prec_series)
names(prec_series) <- old_names
rm(old_names)

## ----eval=download_data-------------------------------------------------------
# download_dataset(dataset = "WorldClim_2.1_10m", bio_variables = prec_vars)
# prec_obs_hres_all <- region_series(
#   bio_variables = prec_vars,
#   time_ce = 1985,
#   dataset = "WorldClim_2.1_10m",
#   ext = n_europe_ext
# )

## ----echo=FALSE, results="hide", eval=create_custom_data----------------------
# terra::saveRDS(prec_obs_hres_all,
#   file = "../inst/extdata/delta/prec_obs_hres_all.RDS"
# )

## ----echo=FALSE, results="hide", eval=!download_data--------------------------
prec_obs_hres_all <- terra::readRDS(
  system.file("extdata/delta/prec_obs_hres_all.RDS",
    package = "pastclim"
  )
)

## -----------------------------------------------------------------------------
prec_obs_range <- range(
  unlist(
    lapply(prec_obs_hres_all, minmax,
      compute = TRUE
    )
  )
)
prec_obs_range

## -----------------------------------------------------------------------------
prec_downscaled_list <- list()
for (i in 1:12) {
  delta_rast <- delta_compute(
    x = prec_series[[i]], ref_time = 0,
    obs = prec_obs_hres_all[[i]]
  )
  prec_downscaled_list[[i]] <- delta_downscale(
    x = prec_series[[i]],
    delta_rast = delta_rast,
    x_landmask_high = land_mask_high_res,
    range_limits = prec_obs_range
  )
}
prec_downscaled <- terra::sds(prec_downscaled_list)

## -----------------------------------------------------------------------------
bioclim_downscaled <- bioclim_vars(
  tavg = tavg_downscaled,
  prec = prec_downscaled
)

## -----------------------------------------------------------------------------
bioclim_downscaled

## -----------------------------------------------------------------------------
panel(bioclim_downscaled[[1]], main = time_bp(bioclim_downscaled[[1]]))

## -----------------------------------------------------------------------------
terra::writeCDF(bioclim_downscaled,
  paste0(tempdir(), "/EA_bioclim_downscaled.nc"),
  overwrite = TRUE
)

## -----------------------------------------------------------------------------
custom_data <- region_series(
  bio_variables = c("bio01", "bio04", "bio19"),
  dataset = "custom",
  path_to_nc = paste0(tempdir(), "/EA_bioclim_downscaled.nc")
)

## -----------------------------------------------------------------------------
custom_data

## -----------------------------------------------------------------------------
panel(custom_data$bio01, main = time_bp(custom_data$bio01))

