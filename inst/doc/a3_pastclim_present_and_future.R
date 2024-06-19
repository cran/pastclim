## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# boolean to turn off chunks that download large datasets
eval_chunk <- FALSE 

## -----------------------------------------------------------------------------
library(pastclim)
get_vars_for_dataset("WorldClim_2.1_10m")

## -----------------------------------------------------------------------------
get_vars_for_dataset("WorldClim_2.1_10m", monthly = TRUE, annual = FALSE)

## ----eval=eval_chunk----------------------------------------------------------
#  download_dataset(
#    dataset = "WorldClim_2.1_10m",
#    bio_variables = c("bio01", "bio02", "altitude")
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  climate_present <- region_slice(
#    time_ce = 1985,
#    bio_variables = c("bio01", "bio02", "altitude"),
#    dataset = "WorldClim_2.1_10m"
#  )

## -----------------------------------------------------------------------------
library(pastclim)
get_vars_for_dataset("CHELSA_2.1_0.5m")

## -----------------------------------------------------------------------------
get_vars_for_dataset("CHELSA_2.1_0.5m", monthly = TRUE, annual = FALSE)

## ----eval=eval_chunk----------------------------------------------------------
#  download_dataset(
#    dataset = "CHELSA_2.1_0.5m",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  climate_present <- region_slice(
#    time_ce = 1990,
#    bio_variables = c("bio01", "bio02"),
#    dataset = "CHELSA_2.1_0.5m"
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  download_dataset(dataset = "CHELSA_2.1_0.5m_vsi",
#                   bio_variables = c("bio12","temperature_01"))

## ----eval=eval_chunk----------------------------------------------------------
#  climate_present <- region_slice(
#    time_ce = 1990,
#    bio_variables = c("bio12","temperature_01"),
#    dataset = "CHELSA_2.1_0.5m_vsi"
#  )

## -----------------------------------------------------------------------------
list_available_datasets()[grepl("WorldClim_2.1",list_available_datasets())]

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m")

## ----eval=eval_chunk----------------------------------------------------------
#  download_dataset(
#    dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  future_slice <- region_slice(
#    time_ce = 2030,
#    dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  future_series <- region_series(
#    dataset = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    bio_variables = c("bio01", "bio02")
#  )

## -----------------------------------------------------------------------------
list_available_datasets()[grepl("CHELSA_2.1",list_available_datasets())]

## -----------------------------------------------------------------------------
get_vars_for_dataset(dataset = "CHELSA_2.1_UKESM1-0-LL_ssp585_0.5m", monthly=TRUE)

## ----eval=eval_chunk----------------------------------------------------------
#  download_dataset(
#    dataset = "CHELSA_2.1_UKESM1-0-LL_ssp585_0.5m_vsi",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  future_slice <- region_slice(
#    time_ce = 2025,
#    dataset = "CHELSA_2.1_UKESM1-0-LL_ssp585_0.5m_vsi",
#    bio_variables = c("bio01", "bio02")
#  )

## ----eval=eval_chunk----------------------------------------------------------
#  future_series <- region_series(
#    dataset = "CHELSA_2.1_UKESM1-0-LL_ssp585_0.5m_vsi",
#    bio_variables = c("bio01", "bio02")
#  )

