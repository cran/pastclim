## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
#  devtools::install_github("EvolEcolGroup/pastclim")

## ----install_vignette, eval=FALSE---------------------------------------------
#  devtools::install_github("EvolEcolGroup/pastclim", build_vignettes = TRUE)

## ----vignette, eval=FALSE-----------------------------------------------------
#  vignette("pastclim_overview", package = "pastclim")

## ----eval=FALSE---------------------------------------------------------------
#  vignette("custom_datasets", package = "pastclim")

## ----eval=FALSE---------------------------------------------------------------
#  vignette("available_datasets", package = "pastclim")

## -----------------------------------------------------------------------------
citation("pastclim")

## ----eval=FALSE---------------------------------------------------------------
#  help("Beyer2020")

## ----echo=FALSE---------------------------------------------------------------
pastclim:::get_dataset_info(dataset="Beyer2020")

## ----eval=FALSE---------------------------------------------------------------
#  library(pastclim)
#  set_data_path()

## ----echo=FALSE, results='hide'-----------------------------------------------
library(pastclim)
set_data_path_for_CRAN()

## ----eval=FALSE---------------------------------------------------------------
#  set_data_path(path_to_nc = "~/my_reconstructions")

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset = "Example")

## -----------------------------------------------------------------------------
get_time_steps(dataset = "Example")

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset = "Beyer2020")

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset = "Krapp2021")

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset="Example", details = TRUE)

## -----------------------------------------------------------------------------
get_downloaded_datasets()

## ----eval=FALSE---------------------------------------------------------------
#  download_dataset(dataset = "Beyer2020", bio_variables = c("bio01", "bio05"))

## -----------------------------------------------------------------------------
locations <- data.frame(
  name = c("Iho Eleru","La Riera", "Chalki", "Oronsay","Atlantis"), 
  longitude = c(5,-4, 27, -6, -24), latitude = c(7, 44, 36, 56, 31),
  time_bp = c(-11200, -18738,-10227, -10200, -11600)
)
locations

## -----------------------------------------------------------------------------
location_slice(
  x = locations, bio_variables = c("bio01", "bio12"),
  dataset = "Example", nn_interpol = FALSE
)

## -----------------------------------------------------------------------------
location_slice(
  x = locations, bio_variables = c("bio01", "bio12"),
  dataset = "Example", nn_interpol = TRUE)

## -----------------------------------------------------------------------------
locations_ts <- location_series(
  x = locations,
  bio_variables = c("bio01", "bio12"),
  dataset = "Example")

## -----------------------------------------------------------------------------
subset(locations_ts, name == "Iho Eleru")

## -----------------------------------------------------------------------------
subset(locations_ts, name == "Oronsay")

## ---- warning=TRUE, fig.width=4, fig.height=3---------------------------------
library(ggplot2)
ggplot(data=locations_ts, aes(x=time_bp, y=bio01, group=name)) +
  geom_line(aes(col=name))+
  geom_point(aes(col=name))

## ---- warning=TRUE, fig.width=4, fig.height=3---------------------------------
library(ggplot2)
ggplot(data=locations_ts, aes(x=time_bp, y=bio01, group=name)) +
  geom_line(aes(col=name))+
  geom_point(aes(col=name))+
  labs(y = var_labels("bio01", dataset="Example", abbreviated=TRUE),
       x = "time BP (yr)")

## -----------------------------------------------------------------------------
climate_20k <- region_slice(
  time_bp = -20000,
  bio_variables = c("bio01", "bio10", "bio12"),
  dataset = "Example"
)

## -----------------------------------------------------------------------------
library(terra)
climate_20k

## ---- fig.width=6, fig.height=5-----------------------------------------------
terra::plot(climate_20k, main = names(climate_20k))

## ---- fig.width=6, fig.height=5-----------------------------------------------
terra::plot(climate_20k, main = var_labels(climate_20k, dataset = "Example"))

## -----------------------------------------------------------------------------
climate_region <- region_series(
  time_bp = list(min = -15000, max = 0), 
  bio_variables = c("bio01", "bio10", "bio12"),
  dataset = "Example")
climate_region

## -----------------------------------------------------------------------------
climate_region$bio01

## ---- fig.width=6, fig.height=5-----------------------------------------------
terra::plot(climate_region$bio01, main=time_bp(climate_region$bio01))

## ---- fig.width=6, fig.height=5-----------------------------------------------
slice_10k <- slice_region_series(climate_region, time_bp = -10000)
terra::plot(slice_10k, main = names(slice_10k))

## -----------------------------------------------------------------------------
mis1_steps <- get_mis_time_steps(mis = 1, dataset = "Example")
mis1_steps

## -----------------------------------------------------------------------------
climate_mis1 <- region_series(
  time_bp = mis1_steps, 
  bio_variables = c("bio01", "bio10", "bio12"),
  dataset = "Example"
)
climate_mis1

## -----------------------------------------------------------------------------
names(region_extent)

## -----------------------------------------------------------------------------
region_extent$Europe

## ---- fig.width=6, fig.height=5-----------------------------------------------
europe_climate_20k <- region_slice(
  time_bp = -20000,
  bio_variables = c("bio01", "bio10", "bio12"),
  dataset = "Example",
  ext = region_extent$Europe
)
terra::plot(europe_climate_20k,
            main = names(europe_climate_20k))

## -----------------------------------------------------------------------------
names(region_outline)

## ---- fig.width=6, fig.height=5-----------------------------------------------
europe_climate_20k <- region_slice(
  time_bp = -20000,
  c("bio01", "bio10", "bio12"),
  dataset = "Example",
  crop = region_outline$Europe
)
terra::plot(europe_climate_20k,
            main = names(europe_climate_20k))

## ---- fig.width=6, fig.height=5-----------------------------------------------
library(sf)
afr_eurasia <- sf::st_union(region_outline$Africa, region_outline$Eurasia)
climate_20k_afr_eurasia <- region_slice(
  time_bp = -20000,
  bio_variables  = c("bio01", "bio10", "bio12"),
  dataset = "Example",
  crop = afr_eurasia)
terra::plot(climate_20k_afr_eurasia,
            main = names(climate_20k_afr_eurasia))

## ---- fig.width=6, fig.height=5-----------------------------------------------
custom_vec <- terra::vect("POLYGON ((0 70, 25 70, 50 80, 170 80, 170 10,
                              119 2.4, 119 0.8, 116 -7.6, 114 -12, 100 -40,
                              -25 -40, -25 64, 0 70))")
climate_20k_custom <- region_slice(
  time_bp = -20000,
  bio_variables = c("bio01", "bio10", "bio12"),
  dataset = "Example",
  crop = custom_vec)
terra::plot(climate_20k_custom, 
            main = names(climate_20k_custom))

## -----------------------------------------------------------------------------
get_biome_classes("Example")

## ---- fig.width=6, fig.height=5-----------------------------------------------
climate_20k <- region_slice(
  time_bp = -20000,
  bio_variables = c("bio01", "bio10", "biome"),
  dataset = "Example"
)
climate_20k$desert <- climate_20k$biome
climate_20k$desert[climate_20k$desert != 21] <- FALSE
climate_20k$desert[climate_20k$desert == 21] <- TRUE
terra::plot(climate_20k, main = names(climate_20k))

## ---- fig.width=6, fig.height=5-----------------------------------------------
climate_20k <- region_slice(
  time_bp = -20000,
  bio_variables = c("bio01", "bio10"),
  dataset = "Example"
)
climate_20k$ice_mask <- get_ice_mask(-20000, dataset = "Example")
climate_20k$land_mask <- get_land_mask(-20000, dataset = "Example")
terra::plot(climate_20k, main = names(climate_20k))

## ---- fig.width=4, fig.height=4-----------------------------------------------
bio_vars <- c("bio01", "bio10", "bio12")
climate_10k <- region_slice(-10000,
  bio_variables = bio_vars,
  dataset = "Example"
)
climate_values_10k <- df_from_region_slice(climate_10k)
climate_10k_pca <- prcomp(climate_values_10k[, bio_vars], 
                          scale = TRUE, center = TRUE)
plot(climate_10k_pca$x[, 2] ~ climate_10k_pca$x[, 1],
  pch = 20, col = "lightgray",
  xlab = "PC1", ylab = "PC2"
)

## -----------------------------------------------------------------------------
locations_10k <- data.frame(
  longitude = c(0, 90, 20, 5), latitude = c(20, 45, 50, 47),
  time_bp = c(-9932, -9753, -10084, -10249)
)
climate_loc_10k <- location_slice(
  x = locations_10k[, c("longitude", "latitude")],
  time_bp = locations_10k$time_bp, bio_variables = bio_vars,
  dataset = "Example"
)
locations_10k_pca_scores <- predict(climate_10k_pca, 
                                    newdata = climate_loc_10k[, bio_vars])

## ---- fig.width=4, fig.height=4-----------------------------------------------
plot(climate_10k_pca$x[, 2] ~ climate_10k_pca$x[, 1],
  pch = 20, col = "lightgray",
  xlab = "PC1", ylab = "PC2"
)
points(locations_10k_pca_scores, pch = 20, col = "red")

## -----------------------------------------------------------------------------
climate_20k <- region_slice(
  time_bp = -20000,
  bio_variables = c("bio01", "bio10"),
  dataset = "Example"
)
this_sample <- sample_region_slice(climate_20k, size = 100)
head(this_sample)

## -----------------------------------------------------------------------------
climate_ts <- region_series(
  time_bp = c(-20000,-10000),
  bio_variables = c("bio01", "bio10", "bio12"),
  dataset = "Example",
  ext = terra::ext(region_extent$Europe)
)
sampled_climate <- sample_region_series(climate_ts, size = c(3,5))
sampled_climate

## ---- fig.width=6, fig.height=4-----------------------------------------------
Namerica_10k <- region_slice(dataset="Example", c("bio01"),
                               time_bp=-10000, ext=region_extent$N_America)
terra::plot(Namerica_10k)

## ---- fig.width=6, fig.height=4-----------------------------------------------
Namerica_ds <- terra::disagg(Namerica_10k, fact=25, method='bilinear')
terra::plot(Namerica_ds)


