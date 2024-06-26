# this file tests the setting and getting of metadata from a vrt file
# we work in the temp directory
skip_on_cran()
# buildvrt from gdal is not properly installed on the linux machine on CRAN
test_that("setting and getting vrt meta", {
  vrt_path <- file.path(tempdir(),"test.vrt")
  tif_files <- list.files(system.file("extdata/CHELSA_bio01", package="pastclim"),
             full.names = TRUE)
  # create the file
  # vrt_path <- terra::vrt(x = tif_files,
  #                        filename = vrt_path,
  #                        options="-separate", overwrite=TRUE, return_filename=TRUE)
  sf::gdal_utils(
    util = "buildvrt",
    source = tif_files,
    destination = vrt_path,
    options = c("-separate","-overwrite")
  )
  
  
  description <- "band_name_1"
  time_vector <- c(0,-10,-1000)
  expect_true(vrt_set_meta(vrt_path, description, time_vector))
  # check we have the correct description in the file
  vrt_rast <- terra::rast(vrt_path)
  expect_true(identical(names(vrt_rast),paste(description,time_vector, sep="_")))
  vrt_meta <- vrt_get_meta(vrt_path = vrt_path)
  expect_true(identical(vrt_meta$time_bp,time_vector))
  expect_true(identical(vrt_meta$description,description))
  # check that get_time_bp_steps works with a vrt file
  expect_true(identical(
    get_time_bp_steps(dataset="custom", path_to_nc = vrt_path),time_vector))
  
  
  # expect a warning if we try to set the metadata a second time
  expect_warning(vrt_res <- vrt_set_meta(vrt_path, description, time_vector),
               "metadata for pastclim is already present")
  expect_false(vrt_res)
  
  # expect a warning if we pass an incorrect number of times
  # vrt_path <- terra::vrt(x = tif_files,
  #                        filename = vrt_path,
  #                        options="-separate", overwrite=TRUE, return_filename=TRUE)
  sf::gdal_utils(
    util = "buildvrt",
    source = tif_files,
    destination = vrt_path,
    options = c("-separate","-overwrite")
  )
  expect_warning(vrt_res <- vrt_set_meta(vrt_path, description, c(time_vector,4)),
                 "the vrt has a different number of")
  expect_false(vrt_res)
  
}
)
