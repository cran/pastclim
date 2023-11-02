## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(pastclim)
get_vars_for_dataset("WorldClim_2.1_10m")

## -----------------------------------------------------------------------------
get_vars_for_dataset("WorldClim_2.1_10m", monthly = TRUE, annual = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  download_dataset(
#    dataset = "WorldClim_2.1_10m",
#    bio_variables = c("bio01", "bio02", "altitude")
#  )

## ----eval=FALSE---------------------------------------------------------------
#  climate_present <- region_slice(
#    time_ce = 1985,
#    bio_variables = c("bio01", "bio02", "altitude"),
#    dataset = "WorldClim_2.1_10m"
#  )

## -----------------------------------------------------------------------------
list_available_datasets()

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m")

## ----eval=FALSE---------------------------------------------------------------
#  download_dataset(
#    dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=FALSE---------------------------------------------------------------
#  future_slice <- region_slice(
#    time_ce = 2030,
#    dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=FALSE---------------------------------------------------------------
#  future_series <- region_series(
#    dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    bio_variables = c("bio01", "bio02")
#  )

