#' Credible regions of location surfaces
#'
#' @description
#' For each animal and bin, the smallest set of cells that holds `level` of
#' the probability of a location surface: the cells where the animal most
#' probably was, and how large that area is.
#'
#' @details
#' Cells enter in order of probability until their sum reaches `level`; ties
#' go to the cell listed first. A surface stores only the cells that hold its
#' kept mass (0.999 by default), so `level` cannot exceed it. A surface cut
#' down by subsetting may store less than `level` in some bins; their regions
#' then hold every stored cell, and a warning gives their number.
#'
#' The area is the number of cells times the cell area of the surface's grid,
#' in square metres; on a longitude-latitude grid each cell's own area is
#' summed.
#'
#' @param surface A location surface from [location_surface()].
#' @param level Share of the probability the region holds, above 0 and at most
#'   the surface's kept mass.
#' @param raster If `TRUE`, return a `terra` raster with one layer per animal
#'   and bin, holding each cell's probability inside the region and `NA`
#'   outside it, instead of the table.
#' @return A data frame of class `"credible_region"` with the columns
#'   `animal_id`, `bin_start`, `cell` and `prob`, one row per cell in a
#'   region, and the attribute `summary`: one row per animal and bin with
#'   `n_cells`, `area_m2` and `mass`, the probability the region holds. With
#'   `raster = TRUE`, a `SpatRaster`.
#' @seealso [location_surface()]
#' @export
#' @examples
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Beyond 2 km the Stoney curve is below 0.001 at every receiver; leaving
#' # those pairs out keeps the example quick.
#' dist <- water_distances(grid, stoney_track$receivers,
#'                         max_distance = 2000)
#' model <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' field <- de_field(dist, de_model(model, sigma = list(station_hour = 0)))
#' counts <- bin_detections(stoney_track$detections, stoney_track$receivers,
#'                          stoney_track$tags)
#' surface <- location_surface(counts, field)
#' region <- credible_region(surface, level = 0.95)
#' attr(region, "summary")
#' # As a raster, one layer per animal and hour. Each cell in the region holds
#' # its probability; cells outside are NA.
#' r <- credible_region(surface, raster = TRUE)
#' terra::nlyr(r)
#' names(r)[1:3]
credible_region <- function(surface, level = 0.95, raster = FALSE) {
  sl_check_surface(surface)
  kept <- attr(surface, "kept")$mass
  if (!is.numeric(level) || length(level) != 1 || !(level > 0) || level > kept) {
    stop(sprintf("`level` must be above 0 and at most the kept mass of the surface, %s.",
                 format(kept)), call. = FALSE)
  }
  s <- data.table::data.table(animal_id = surface$animal_id, bin_start = surface$bin_start,
                              cell = surface$cell, prob = surface$prob)
  before <- prob <- NULL
  data.table::setorderv(s, c("animal_id", "bin_start", "prob", "cell"), c(1, 1, -1, 1))
  s[, before := cumsum(prob) - prob, by = c("animal_id", "bin_start")]
  cr_check_mass(s, level)
  region <- s[before < level * (1 - 1e-12)]
  region[, before := NULL]
  grid <- sl_grid(surface)
  if (isTRUE(raster)) return(cr_raster(region, grid))
  cr_table(region, grid, level)
}

# A bin that stores less than the level, as after subsetting, cannot give a
# region holding the level.
cr_check_mass <- function(s, level) {
  prob <- NULL
  stored <- s[, list(mass = sum(prob)), by = c("animal_id", "bin_start")]$mass
  short <- stored < level * (1 - 1e-9)
  if (any(short)) {
    warning(sprintf(paste0("%s of %s bins store less than the level %s (lowest %s), as after ",
                           "subsetting a surface; their regions hold every stored cell and ",
                           "less than the level."), fmt_n(sum(short)),
                    fmt_n(length(stored)), format(level),
                    format(signif(min(stored), 3))), call. = FALSE)
  }
}

cr_table <- function(region, grid, level) {
  cell_area <- cr_cell_area(grid, region$cell)
  prob <- area <- NULL
  region[, area := cell_area]
  summ <- region[, list(n_cells = .N, area_m2 = sum(area), mass = sum(prob)),
                 by = c("animal_id", "bin_start")]
  out <- data.frame(animal_id = region$animal_id, bin_start = region$bin_start,
                    cell = region$cell, prob = region$prob, stringsAsFactors = FALSE)
  structure(out, class = c("credible_region", "data.frame"), level = level,
            summary = as.data.frame(summ))
}

# Cell area in square metres: one value on a projected grid, each cell's own
# on a longitude-latitude grid, NA without a grid.
cr_cell_area <- function(grid, cells) {
  if (is.null(grid)) return(rep(NA_real_, length(cells)))
  if (!terra::is.lonlat(grid)) return(rep(prod(terra::res(grid)), length(cells)))
  terra::cellSize(grid, unit = "m")[cells][[1]]
}

cr_raster <- function(region, grid) {
  if (is.null(grid)) stop("The surface holds no grid, so no raster can be made.", call. = FALSE)
  groups <- unique(region[, c("animal_id", "bin_start")])
  template <- terra::rast(grid, nlyrs = 1)
  layers <- lapply(seq_len(nrow(groups)), function(i) {
    r <- region[region$animal_id == groups$animal_id[i] & region$bin_start == groups$bin_start[i]]
    v <- rep(NA_real_, terra::ncell(template))
    v[r$cell] <- r$prob
    terra::setValues(template, v)
  })
  out <- terra::rast(layers)
  names(out) <- paste(groups$animal_id, ls_time(groups$bin_start))
  out
}

#' @export
print.credible_region <- function(x, ...) {
  s <- attr(x, "summary")
  cat(sprintf("Credible regions at level %s: %s, median %s cells (%s m^2)\n",
              format(attr(x, "level")), plural(nrow(s), "animal-bin"),
              fmt_n(stats::median(s$n_cells)),
              fmt_n(stats::median(s$area_m2))))
  print(utils::head(as.data.frame(unclass(x)[c("animal_id", "bin_start", "cell", "prob")]),
                    n = 6))
  invisible(x)
}
