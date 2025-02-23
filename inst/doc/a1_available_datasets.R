## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(pastclim)

## ----echo=FALSE, results='hide'-----------------------------------------------
data_path <- file.path(tempdir(), "pastclim_data")
# clear it in case it exists already
unlink(data_path, recursive = TRUE)
# set data path
set_data_path(
  path_to_nc = data_path,
  ask = FALSE,
  write_config = FALSE,
  copy_example = TRUE
)

## -----------------------------------------------------------------------------
get_available_datasets()

## -----------------------------------------------------------------------------
list_available_datasets()

## ----eval=FALSE---------------------------------------------------------------
# help("Example")

## ----echo=FALSE---------------------------------------------------------------
pastclim:::get_dataset_info(dataset = "Example")

## ----echo=FALSE---------------------------------------------------------------
list_datasets <- list_available_datasets()
# replace datasets with multiple versions with a single string
list_datasets <- c(
  list_datasets[!grepl("WorldClim_2.1", list_datasets)],
  "WorldClim_2.1"
)
list_datasets <- c(
  list_datasets[!grepl("paleoclim_1.0", list_datasets)],
  "paleoclim_1.0"
)
list_datasets <- c(
  list_datasets[!grepl("CHELSA_2.1", list_datasets)],
  "CHELSA_2.1"
)
for (i in sort(list_datasets)) {
  pastclim:::get_dataset_info(i)
  cat("\n#######################################################\n")
}

