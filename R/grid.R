# ---- Water grid helpers ----
# The grid travels inside package objects as a PackedSpatRaster, because a
# SpatRaster points to C++ memory that saveRDS() cannot keep.

is_one_number <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x)

is_one_string <- function(x) is.character(x) && length(x) == 1 && !is.na(x)

as_water_grid <- function(grid, call = NULL) {
  if (inherits(grid, "PackedSpatRaster")) grid <- terra::unwrap(grid)
  if (!inherits(grid, "SpatRaster")) {
    grid <- tryCatch(terra::rast(grid), error = function(e) {
      detfield_abort(paste("`grid` must be a SpatRaster, a raster file path or an object",
                           "terra::rast() reads."),
                     "grid", call = call)
    })
  }
  if (terra::nlyr(grid) > 1) grid <- grid[[1]]
  if (!nzchar(terra::crs(grid))) {
    detfield_abort("`grid` has no CRS. Set a projected CRS in metres with terra::crs().",
                   "grid", call = call)
  }
  if (terra::is.lonlat(grid)) {
    detfield_abort(wrap_text(paste("`grid` is in longitude and latitude.",
                                   "Project it to a CRS in metres first.")),
                   "grid", call = call)
  }
  unit_m <- terra::linearUnits(grid)
  if (!isTRUE(all.equal(unit_m, 1))) {
    detfield_abort(wrap_text(sprintf(paste("`grid` must have a CRS in metres; one unit of its",
                                           "CRS is %s m. Project it to a CRS in metres first."),
                                     format(signif(unit_m, 7)))),
                   "grid", call = call)
  }
  check_scale(grid, call)
  grid
}

# ---- CRS scale ----
# The largest departure from 1 that the scale of a projected CRS may have at
# the raster centre. UTM stays within 0.1 percent across a zone.
max_scale_error <- 0.01

# Refuses a projected CRS whose scale at the raster centre differs from 1 by
# more than max_scale_error on either axis or either diagonal. Distances are
# read as metres on the ground, so a grid in Web Mercator gives distances too
# long by a factor of about 1 / cos(latitude). At 45 N a sinusoidal CRS
# keeps both axes within 1 percent up to about 11 degrees from its central
# meridian, and departs on the diagonals. A CRS with no transformation to
# longitude and latitude, such as a local engineering CRS, has no scale to
# measure. Nothing then shows that its metres are metres on the ground, so it
# is refused too.
check_scale <- function(grid, call = NULL) {
  s <- crs_scale(grid)
  if (anyNA(s$scale)) {
    detfield_abort(wrap_text(paste(
      "The scale of this CRS at the raster centre cannot be measured, so give the",
      "raster in a projected CRS with a geographic datum, such as its UTM zone."
    )), "grid", call = call)
  }
  if (all(abs(s$scale - 1) <= max_scale_error)) return(invisible())
  shown <- sprintf("%.3f", s$scale)
  scale_text <- if (length(unique(shown)) == 1L) {
    paste(shown[1], "on its two axes and its diagonals")
  } else {
    paste0(shown[1], " east-west, ", shown[2], " north-south, ", shown[3],
           " south-west to north-east and ", shown[4], " north-west to south-east")
  }
  epsg <- utm_epsg(s$lon, s$lat)
  detfield_abort(wrap_text(paste0(
    "`grid` has a CRS whose scale at the raster centre is ", scale_text,
    ", more than 1 percent from 1, so its metres are not metres on the ground. ",
    "Project it to the UTM zone of the centre, ", epsg %% 100,
    if (epsg < 32700) "N" else "S", " (EPSG:", epsg, ")."
  )), "grid", call = call)
}

# The scale of a grid's CRS at the raster centre: the grid length of four 1 km
# segments centred there (east-west, north-south and the two diagonals), each
# over its geodesic length. Also returns the centre in longitude and latitude.
crs_scale <- function(grid, half = 500) {
  e <- as.vector(terra::ext(grid))
  cx <- (e[["xmin"]] + e[["xmax"]]) / 2
  cy <- (e[["ymin"]] + e[["ymax"]]) / 2
  d <- half / sqrt(2)
  from <- rbind(c(cx - half, cy), c(cx, cy - half), c(cx - d, cy - d), c(cx - d, cy + d))
  to <- rbind(c(cx + half, cy), c(cx, cy + half), c(cx + d, cy + d), c(cx + d, cy - d))
  # When no transformation exists, terra warns and may then stop with an error.
  # The warning is muffled and the error gives a matrix of NA. An NA is the
  # answer, and check_scale() raises the one condition the user sees.
  ll <- tryCatch(suppressWarnings(terra::project(rbind(c(cx, cy), from, to),
                                                 from = terra::crs(grid), to = "EPSG:4326")),
                 error = function(e) matrix(NA_real_, 9, 2))
  if (anyNA(ll) || any(!is.finite(ll))) {
    return(list(scale = rep(NA_real_, 4), lon = NA_real_, lat = NA_real_))
  }
  geo <- terra::distance(ll[2:5, , drop = FALSE], ll[6:9, , drop = FALSE],
                         lonlat = TRUE, pairwise = TRUE)
  list(scale = 2 * half / geo, lon = ll[1, 1], lat = ll[1, 2])
}

# EPSG code of the WGS 84 UTM zone that holds a point: 326zz north of the
# equator, 327zz south. A longitude of exactly 180 falls in zone 60.
utm_epsg <- function(lon, lat) {
  zone <- min(floor((lon + 180) / 6) + 1, 60)
  as.integer(if (lat >= 0) 32600 + zone else 32700 + zone)
}

grid_info <- function(grid) {
  e <- as.vector(terra::ext(grid))
  list(
    xmin = e[["xmin"]], ymax = e[["ymax"]], res = terra::res(grid),
    ncol = terra::ncol(grid), nrow = terra::nrow(grid),
    water = !is.na(terra::values(grid, mat = FALSE))
  )
}

# Only NA cells count as land. A bathymetry raster that codes land as 0, or
# by the sign of an elevation, would otherwise give a lake with no shore.
check_land_coding <- function(grid) {
  v <- terra::values(grid, mat = FALSE)
  n_na <- sum(is.na(v))
  if (n_na == 0) {
    warning(sprintf("`grid` has no NA cells, so all %s cells count as water. %s",
                    fmt_n(length(v)), "Set land cells to NA."), call. = FALSE)
    return(invisible(FALSE))
  }
  low <- sum(v <= 0, na.rm = TRUE)
  high <- sum(v > 0, na.rm = TRUE)
  if (low > 0 && high > 0) {
    warning(sprintf(paste("%s of %s water cells hold values of 0 or less and %s hold values",
                          "above 0. Only NA cells count as land; if 0 or a sign marks land,",
                          "set those cells to NA."),
                    fmt_n(low), fmt_n(low + high),
                    fmt_n(high)), call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# The EPSG code of a grid's CRS; NA when there is no grid or the CRS has no
# numeric EPSG code.
grid_epsg <- function(grid) {
  if (!inherits(grid, "SpatRaster")) return(NA_integer_)
  code <- terra::crs(grid, describe = TRUE)$code
  if (length(code) && !is.na(code) && grepl("^[0-9]+$", code)) as.integer(code) else NA_integer_
}

# Cell numbers follow terra: row-major from the top-left corner, starting at 1.
# A point outside the grid gets NA.
cell_of_xy <- function(info, x, y) {
  col <- floor((x - info$xmin) / info$res[1]) + 1
  row <- floor((info$ymax - y) / info$res[2]) + 1
  out <- (row - 1) * info$ncol + col
  out[col < 1 | col > info$ncol | row < 1 | row > info$nrow] <- NA
  out
}

xy_of_cell <- function(info, cell) {
  row <- (cell - 1) %/% info$ncol + 1
  col <- (cell - 1) %% info$ncol + 1
  cbind(
    x = info$xmin + (col - 0.5) * info$res[1],
    y = info$ymax - (row - 0.5) * info$res[2]
  )
}

# A raster of `values` on the water cells of `grid`; other cells are NA.
cells_to_raster <- function(grid, cells, values, names = NULL) {
  grid <- as_water_grid(grid)
  values <- as.matrix(values)
  out <- terra::rast(grid, nlyrs = ncol(values))
  m <- matrix(NA_real_, terra::ncell(grid), ncol(values))
  m[cells, ] <- values
  terra::values(out) <- m
  if (!is.null(names)) names(out) <- names
  out
}
