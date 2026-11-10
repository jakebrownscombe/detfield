# ---- Receiver positions on the water grid ----

# Coordinates of each row of `receivers` in the CRS of the grid, as columns x
# and y. Projected x and y win over latitude and longitude when both exist,
# because they carry no datum shift.
receiver_xy <- function(receivers, grid) {
  if (inherits(receivers, "sf")) {
    pts <- sf::st_transform(sf::st_geometry(receivers), terra::crs(grid))
    xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
    receivers <- sf::st_drop_geometry(receivers)
    receivers$x <- xy[, 1]
    receivers$y <- xy[, 2]
    return(as.data.frame(receivers))
  }
  if (!is.data.frame(receivers)) {
    stop("`receivers` must be a data frame or an sf object.", call. = FALSE)
  }
  receivers <- as.data.frame(receivers)
  if (all(c("x", "y") %in% names(receivers))) {
    check_epsg(attr(receivers, "epsg"), grid_epsg(grid))
    check_numeric(receivers, c("x", "y"), "receivers")
    return(receivers)
  }
  if (!all(c("deploy_lat", "deploy_lon") %in% names(receivers))) {
    stop("`receivers` needs columns x and y in the CRS of `grid`, or deploy_lat and deploy_lon.",
         call. = FALSE)
  }
  check_numeric(receivers, c("deploy_lat", "deploy_lon"), "receivers")
  xy <- terra::project(cbind(receivers$deploy_lon, receivers$deploy_lat),
                       from = "EPSG:4326", to = terra::crs(grid))
  receivers$x <- xy[, 1]
  receivers$y <- xy[, 2]
  receivers
}

check_epsg <- function(table_epsg, grid_epsg) {
  if (is.null(table_epsg) || is.na(grid_epsg)) return(invisible(TRUE))
  if (!identical(as.integer(table_epsg), grid_epsg)) {
    stop(sprintf("`receivers` has epsg %s and `grid` has EPSG %s; they must match.",
                 table_epsg, grid_epsg), call. = FALSE)
  }
  invisible(TRUE)
}

check_numeric <- function(df, cols, what) {
  for (col in cols) {
    if (!is.numeric(df[[col]]) || anyNA(df[[col]])) {
      stop(sprintf("Column `%s` of `%s` must be numeric with no missing values.", col, what),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

receiver_key <- function(receivers) {
  key <- intersect(c("station_id", "receiver_serial"), names(receivers))[1]
  if (is.na(key)) stop("`receivers` needs a station_id or a receiver_serial column.", call. = FALSE)
  if (!is.character(receivers[[key]]) || anyNA(receivers[[key]])) {
    stop(sprintf("Column `%s` of `receivers` must be character with no missing values.", key),
         call. = FALSE)
  }
  key
}

# One position per key. A station's deployments are averaged; a spread wider
# than a cell means the receivers really moved, which the caller should know.
# A serial is not a place: one that moved between stations would be averaged
# to a point where it never sat, so that case stops.
collapse_receivers <- function(receivers, key, res) {
  id <- receivers[[key]]
  ids <- unique(id)
  g <- match(id, ids)
  n <- tabulate(g, length(ids))
  x <- rowsum(receivers$x, g, reorder = FALSE)[, 1] / n
  y <- rowsum(receivers$y, g, reorder = FALSE)[, 1] / n
  spread <- vapply(seq_along(ids), function(i) {
    max(sqrt((receivers$x[g == i] - x[i])^2 + (receivers$y[g == i] - y[i])^2))
  }, numeric(1))
  out <- data.frame(id = ids, x = unname(x), y = unname(y), n_rows = n, spread_m = spread,
                    stringsAsFactors = FALSE)
  names(out)[1] <- key
  if ("deploy_z" %in% names(receivers) && is.numeric(receivers$deploy_z)) {
    out$deploy_z <- unname(rowsum(receivers$deploy_z, g, reorder = FALSE)[, 1] / n)
  }
  wide <- out$spread_m > max(res)
  if (any(wide) && key == "receiver_serial") {
    stop(sprintf(
      "%s %s at positions more than one cell (%s m) apart: %s. %s",
      plural(sum(wide), "receiver_serial value"), agree(sum(wide), "sits", "sit"), fmt_n(max(res)),
      paste(utils::head(ids[wide], 10), collapse = ", "),
      "Add a station_id column so that each position has its own key."
    ), call. = FALSE)
  }
  if (any(wide)) {
    warning(sprintf(
      "%s %s rows more than one cell (%s m) from %s mean position: %s. %s",
      plural(sum(wide), paste(key, "value")), agree(sum(wide), "has", "have"), fmt_n(max(res)),
      agree(sum(wide), "its", "their"), paste(utils::head(ids[wide], 10), collapse = ", "),
      "The mean is used; give each position its own id if the receivers moved."
    ), call. = FALSE)
  }
  out
}

receiver_lookup <- function(receivers) {
  cols <- intersect(c("receiver_serial", "station_id", "deploy_datetime", "recover_datetime"),
                    names(receivers))
  if (!all(c("receiver_serial", "station_id") %in% cols)) return(NULL)
  for (col in intersect(cols, c("deploy_datetime", "recover_datetime"))) {
    v <- receivers[[col]]
    if (!inherits(v, "POSIXct") || !identical(attr(v, "tzone"), "UTC")) {
      stop(sprintf("Column `%s` of `receivers` must be POSIXct in UTC.", col), call. = FALSE)
    }
  }
  out <- unique(receivers[, cols, drop = FALSE])
  rownames(out) <- NULL
  out
}
