#' Create a land mask
#'
#' Create a land mask for a given time step. The land mask is based on the simple
#' logic of moving the ocean up and down given the current relief profile (
#' topography+bathymetry, i.e. the elevation both above and below sea level).
#' Note that this approach ignores any rebound due to changing mass and distribution of
#' ice sheets.
#' LIMITATIONS: The land mask will show internal lakes/seas as land, as their level
#' is unrelated to the general sea level. If you have specific reconstructions
#' of internal lakes (or want to simply reuse their current extents), you will
#' have to add them onto the masks generated by this function. Also note that the
#' land mask does not include ice sheets. This means that some areas that are permanently
#' covered by ice at the two poles will show up as sea. This means that, for any
#' reconstruction including Greenland or Antarctica, the resulting land mask will
#' need to be modified to include the appropriate ice sheets.
#'
#' @param relief_rast a [`terra::SpatRaster`] with relief
#' @param time_bp the time of interest
#' @param sea_level sea level at the time of interest (if left to NULL, this is
#' computed using Spratt 2016)
#' @returns a [`terra::SpatRaster`] of the land masks (with land as 1's and sea
#' as NAs), where the layers are different times
#'
#' @keywords internal

make_land_mask <- function(relief_rast, time_bp, sea_level = NULL) {
  message("This function is still under development; do not use it for real analysis")
  if (is.null(sea_level)) {
    sea_level <- get_sea_level(time_bp = time_bp)
  } else { # check that we have as many sea level estimates as times
    if (length(time_bp) != length(sea_level)) {
      stop("time_bp and sea_level should have the same number of elements")
    }
  }
  land_mask <- NULL
  for (i in seq_along(time_bp)) {
    # create binary relief map for areas above and below the relevant sea level
    relief_bin <- relief_rast
    relief_bin[relief_bin > sea_level[i]] <- NA
    relief_bin[!is.na(relief_bin)] <- 1
    sea_patches <- patches(relief_bin, directions = 8)
    # get mode of a vector (removing any NAs)
    modal_vector <- function(x) {
      x <- x[!is.na(x)]
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    patch_to_get <- modal_vector(values(sea_patches))
    sea_patches[is.na(sea_patches)] <- -1
    sea_patches[sea_patches != patch_to_get] <- -1
    sea_patches[sea_patches > 0] <- NA
    sea_patches <- -sea_patches
    names(sea_patches) <- "mask"
    if (is.null(land_mask)) {
      land_mask <- sea_patches
    } else {
      add(land_mask) <- sea_patches
    }
  }
  terra::time(land_mask, tstep = "years") <- (time_bp + 1950)
  return(land_mask)
}
