## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
tiffs_path <- system.file("extdata/CHELSA_bio01",package="pastclim")
list_of_tiffs <- file.path(tiffs_path,dir(tiffs_path))
bio01 <- terra::rast(list_of_tiffs)

## -----------------------------------------------------------------------------
terra::time(bio01)<-c(0,-100,-200)
names(bio01)<-paste("bio01",terra::time(bio01),sep="_")

## -----------------------------------------------------------------------------
nc_name <- file.path(tempdir(),"CHELSA_TraCE21k_bio01.nc")
terra::writeCDF(bio01, filename = nc_name, varname = "bio01", overwrite=TRUE)

## -----------------------------------------------------------------------------
nc_in <- ncdf4::nc_open(nc_name, write=TRUE)
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "units",
                 attval = "years since 1950-01-01 00:00:00.0")
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "long_name",
                 attval = "years BP")
ncdf4::ncatt_put(nc_in, varid="time", attname="axis", attval = "T")
ncdf4::nc_close(nc_in)


## -----------------------------------------------------------------------------
library(pastclim)

## ----echo=FALSE, results='hide'-----------------------------------------------
data_path <- file.path(tempdir(),"pastclim_data")
# clear it in case it exists already
unlink(data_path, recursive = TRUE) 
# set data path
set_data_path(path_to_nc = data_path,
              ask = FALSE,
              write_config = FALSE,
              copy_example = TRUE)

## -----------------------------------------------------------------------------
custom_series <- region_series(bio_variables = "bio01",
                                dataset = "custom",
                                path_to_nc = nc_name
)
custom_series

## -----------------------------------------------------------------------------
get_time_steps(dataset="custom", path_to_nc = nc_name)

## -----------------------------------------------------------------------------
climate_100<-slice_region_series(custom_series, time_bp = -100)
terra::plot(climate_100)

## -----------------------------------------------------------------------------
head(read.csv(system.file("extdata/dataset_list_included.csv",
                     package="pastclim")), n=2)

