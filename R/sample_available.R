#' Draw availability points from the detection surface
#'
#' @description
#' For each animal, detected bin and presence draw, draws `n` cells from the
#' whole water body that the array could detect the animal in: in proportion
#' to the probability that the bin is detected from each cell, or uniformly
#' over the cells at or above a minimum probability. The draws are the
#' available half of a used-available design whose used half is
#' [sample_locations()]. positionR 1.x drew availability from the detection
#' surface of the array too (`sample_points_from_system_de()`); the
#' differences are listed under Details.
#'
#' @details
#' **Weighted (`uniform_above = NULL`).** Cells are the rows of `coverage` on
#' the bin's UTC date with `p` above 0, drawn with replacement in proportion to
#' `p`. Presence points exist only for detected bins, and a bin is detected
#' more often where the array hears well. If an animal uses cell x in
#' proportion to u(x), and its bin is detected from x with probability d(x),
#' the presence draws follow u(x) d(x). Availability drawn in proportion to
#' d(x) carries the same factor, so it cancels in the ratio of used to
#' available.
#'
#' **Which coverage.** Pass the probability that a bin is detected:
#' [system_de()] with `transmissions` set to the transmissions a tag sends in
#' one bin (24 an hour at a 150 s mean delay). An hour is a presence when the
#' animal was heard at least once in it, so that is the filter the presence
#' points went through. The per-transmission table weights by a different
#' event. On Stoney Lake, 20 simulated animals with no depth preference on 70
#' receivers gave, with the per-transmission weight, -0.091 per m of bottom
#' depth (95% CI -0.161 to -0.020) on their true tracks and -0.065 (-0.123 to
#' -0.008) from position draws; with the hourly weight, -0.032 (-0.090 to
#' 0.026) and -0.010 (-0.057 to 0.038). These are the means over the animals
#' of block 1 at beta 0 in `stoney_selection$trial`, with t intervals. On
#' that array the hourly weight is close to uniform: 77.7% of cell-days are
#' detected with probability 0.95 or more in an hour.
#'
#' **Uniform (`uniform_above = 0.1`, say).** Cells are the rows of `coverage`
#' on the bin's date with `p` at or above `uniform_above`, drawn uniformly. With
#' the per-transmission table of [system_de()] and 0.1 this follows the
#' example of positionR 1.x, which kept the cells strictly above 0.1. It cuts
#' away water the array barely hears, so used points there have no available
#' counterpart. `uniform_above = 0` takes every row of `coverage`, including
#' cells no receiver can hear; leave out the rows of
#' cells you do not count as available, such as water cut off from the
#' receivers, before the call.
#'
#' **Differences from positionR 1.x.** `sample_points_from_system_de()`
#' weighted each cell by the chance, for one transmission, that one or more
#' receivers detect it (`cumulative_prob`), kept cells strictly above a
#' minimum of 0.001 and strictly below an optional maximum. With a surface
#' carrying fish and time columns it drew one set per fish and time period;
#' with a template of positions and a surface without those columns, it
#' sampled once and copied that one set to every fish and time period. This
#' function weights by the coverage
#' table it is given (hourly, with 24 transmissions, in the article), has no
#' minimum or maximum, keeps cells with `p >=` `uniform_above` in the uniform
#' mode, and draws one set per presence draw.
#'
#' **Why the whole water body.** Availability drawn from a disc round each
#' bin's position can show selection the animals do not have on Stoney Lake:
#' cells with more water within reach are deeper offshore, so a walk with no
#' preference looks like avoidance of deep water. Northrup et al. (2013) found
#' the same kind of bias in simulations when the spatial extent of
#' availability was set wrongly. The coefficient from the whole water body is
#' a lake-scale one, of where in the water body the animal spent its detected
#' bins against where the array could have heard it.
#' `vignette("habitat-selection")` sets it against the true track.
#'
#' Each presence draw gets its own `n` cells. A stratum (one animal, bin and
#' draw) is therefore one matched set of one used and `n` available points,
#' and the strata of one `draw` index form a full data set, as multiple
#' imputation over the presence draws needs. A stratum with no cell to draw,
#' because `coverage` has no row on its date or none that qualifies, is
#' dropped, counted in a warning and listed in the attribute `dropped`.
#'
#' Left to the user: the coverage table, the mode and the model.
#' [used_available()] stacks the used and available points with covariates
#' read at each cell, and fits nothing.
#'
#' @param presence The used points: a data frame with one row per animal, bin
#'   and draw and the columns `animal_id` (character), `bin_start` (POSIXct,
#'   UTC), `draw` (whole numbers), `x` and `y` (in the units of `grid`), such as
#'   the output of [sample_locations()] or true positions sampled by hand. A
#'   `cell` column, if given, must hold the cell of `x` and `y` on `grid`.
#' @param coverage A data frame with the columns `date` (`Date`), `cell` (cells
#'   of `grid`) and `p` (0 to 1), one row per date and cell, such as the table
#'   of [system_de()].
#' @param grid The `SpatRaster` whose cells `coverage` indexes, such as the
#'   lake raster the field was built on.
#' @param n Number of available cells per stratum (animal, bin and draw).
#' @param uniform_above `NULL` to draw cells in proportion to `p`; or one
#'   number in 0 to 1 to draw uniformly among the cells with `p` at or above
#'   it.
#' @param seed An integer seed for the draws, or `NULL` to use the current
#'   random number stream. The caller's stream is restored afterwards when a
#'   seed is given.
#' @return A data frame with `n` rows per kept stratum: `animal_id`,
#'   `bin_start`, `draw` (the presence draw the row is matched to), `cell`, and
#'   `x` and `y`, the cell centre. Attributes: `"estimand"`; `"uniform_above"`
#'   and `"n"` as given; `"epsg"`, the EPSG code of `grid`; and `"dropped"`,
#'   one row per dropped stratum with `animal_id`, `bin_start`, `draw` and
#'   `reason`.
#' @references
#' Northrup, J. M., Hooten, M. B., Anderson, C. R. and Wittemyer, G. (2013).
#' Practical guidance on characterizing availability in resource selection
#' functions under a use-availability design. *Ecology*, 94, 1456-1463.
#' \doi{10.1890/12-1688.1}
#' @seealso [sample_locations()] for the used points, [system_de()] for the
#'   weight, [used_available()] to stack the two.
#' @export
#' @examples
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Beyond 2 km the Stoney curve is below 0.001 at every receiver; leaving
#' # those pairs out keeps the example quick.
#' dist <- water_distances(grid, stoney_track$receivers, max_distance = 2000)
#' model <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' field <- de_field(dist, de_model(model, sigma = list(station_hour = 0)))
#' counts <- bin_detections(stoney_track$detections, stoney_track$receivers,
#'                          stoney_track$tags)
#' surface <- location_surface(counts, field)
#' presence <- sample_locations(surface, n = 20, seed = 1)
#' # Every receiver listens on each day of the study
#' days <- unique(as.Date(presence$bin_start, tz = "UTC"))
#' listening <- expand.grid(station_id = stoney_track$receivers$station_id,
#'                          date = days, stringsAsFactors = FALSE)
#' # Weighted by the probability that an hour is detected: 24 transmissions
#' hourly <- system_de(field, listening, transmissions = 24)
#' available <- sample_available(presence, hourly, grid, n = 20, seed = 1)
#' head(available)
#' # Uniform over the cells heard with probability 0.1 or more per transmission
#' per_transmission <- system_de(field, listening)
#' uniform <- sample_available(presence, per_transmission, grid, n = 20,
#'                             uniform_above = 0.1, seed = 1)
#' nrow(uniform) / nrow(presence)
sample_available <- function(presence, coverage, grid, n, uniform_above = NULL, seed = NULL) {
  # Every step runs under the seed, so a seeded call leaves the session's
  # random stream as it found it.
  with_seed_or_not(seed, sa_run(presence, coverage, grid, n, uniform_above))
}

sa_run <- function(presence, coverage, grid, n, uniform_above) {
  grid <- sa_grid(grid)
  presence <- sa_check_points(presence, grid, "presence")
  coverage <- sa_check_coverage(coverage, grid)
  if (!ls_whole(n) || n < 1) stop("`n` must be one whole number of at least 1.", call. = FALSE)
  if (!is.null(uniform_above)) {
    ok <- is.numeric(uniform_above) && length(uniform_above) == 1 &&
      isTRUE(uniform_above >= 0 && uniform_above <= 1)
    if (!ok) {
      stop("`uniform_above` must be NULL or one number in 0 to 1.", call. = FALSE)
    }
    coverage <- coverage[coverage$p >= uniform_above, , drop = FALSE]
    coverage$p <- rep(1, nrow(coverage))
  }
  strata <- unique(presence[c("animal_id", "bin_start", "draw")])
  strata <- strata[do.call(order, unname(as.list(strata))), , drop = FALSE]
  domains <- sa_domains(strata, coverage, grid)
  key <- domains$key
  dropped <- data.frame(strata[is.na(key), , drop = FALSE], reason = domains$reason[is.na(key)],
                        row.names = NULL, stringsAsFactors = FALSE)
  sa_warn_dropped(dropped, nrow(strata), uniform_above)
  kept <- strata[!is.na(key), , drop = FALSE]
  rows <- sa_draw(domains$cells, key[!is.na(key)], as.integer(n))
  out <- data.frame(animal_id = rep(kept$animal_id, each = n),
                    bin_start = rep(kept$bin_start, each = n),
                    draw = rep(kept$draw, each = n), cell = as.integer(domains$cells$cell[rows]),
                    x = domains$cells$x[rows], y = domains$cells$y[rows],
                    stringsAsFactors = FALSE)
  attr(out, "estimand") <- if (is.null(uniform_above)) {
    paste("cells of the water body available in each detected bin, weighted by the",
          "probability that the bin is detected from the cell")
  } else {
    paste("cells of the water body available in each detected bin, uniform over the cells",
          "with p at or above", format(uniform_above))
  }
  attr(out, "uniform_above") <- uniform_above
  attr(out, "n") <- as.integer(n)
  attr(out, "epsg") <- grid_epsg(grid)
  attr(out, "dropped") <- dropped
  out
}

# ---- Checks ----
sa_grid <- function(grid) {
  if (inherits(grid, "PackedSpatRaster")) grid <- terra::rast(grid)
  if (!inherits(grid, "SpatRaster")) {
    stop("`grid` must be a SpatRaster, the grid whose cells `coverage` indexes.",
         call. = FALSE)
  }
  grid
}

# The used points of sample_available() and used_available(): one row per
# animal, bin and draw, with finite coordinates on the CRS of `grid`. The
# checked columns are returned as a plain data frame.
sa_check_points <- function(points, grid, what) {
  bd_check_columns(points, c(animal_id = "character", bin_start = "POSIXct", draw = "numeric",
                             x = "numeric", y = "numeric"), what)
  bd_check_utc(points, "bin_start", what)
  if (!nrow(points)) stop(sprintf("`%s` has no rows.", what), call. = FALSE)
  if (anyNA(points$animal_id) || anyNA(points$bin_start)) {
    stop(sprintf("`%s$animal_id` and `%s$bin_start` must have no NA.", what, what),
         call. = FALSE)
  }
  if (!all(is.finite(points$draw) & points$draw == round(points$draw) & points$draw >= 1)) {
    stop(sprintf("`%s$draw` must hold whole numbers of at least 1.", what), call. = FALSE)
  }
  if (!all(is.finite(points$x) & is.finite(points$y))) {
    stop(sprintf("`%s$x` and `%s$y` must be finite.", what, what), call. = FALSE)
  }
  repeats <- sum(duplicated(points[c("animal_id", "bin_start", "draw")]))
  if (repeats) {
    stop(wrap_text(sprintf("`%s` must have one row per animal, bin and draw; %s %s one.",
                           what, plural(repeats, "row"), agree(repeats, "repeats", "repeat"))),
         call. = FALSE)
  }
  sa_check_epsg(attr(points, "epsg", exact = TRUE), grid, what)
  keep <- c("animal_id", "bin_start", "draw", "x", "y", if ("cell" %in% names(points)) "cell")
  out <- as.data.frame(points)[keep]
  if ("cell" %in% keep) sa_check_cells(out, grid, what)
  out
}

# CONTRACT.md, section 3: a table combined with a raster stops when their
# EPSG codes differ.
sa_check_epsg <- function(epsg, grid, what) {
  code <- grid_epsg(grid)
  if (length(epsg) == 1 && !is.na(epsg) && !is.na(code) && !identical(as.integer(epsg), code)) {
    stop(sprintf("`%s` is in EPSG:%d and the grid in EPSG:%d.", what, as.integer(epsg), code),
         call. = FALSE)
  }
  invisible(NULL)
}

# A `cell` column must hold the cell of x and y, or the points belong to
# another grid.
sa_check_cells <- function(points, grid, what) {
  at <- terra::cellFromXY(grid, cbind(points$x, points$y))
  bad <- is.na(at) | is.na(points$cell) | at != points$cell
  if (any(bad)) {
    stop(wrap_text(sprintf(paste("%s of %s rows of `%s` have a `cell` that does not hold",
                                 "their x and y on the grid; give the grid the points were",
                                 "drawn on."), fmt_n(sum(bad)),
                           fmt_n(nrow(points)), what)), call. = FALSE)
  }
  invisible(NULL)
}

sa_check_coverage <- function(coverage, grid) {
  bd_check_columns(coverage, c(date = "Date", cell = "numeric", p = "numeric"), "coverage")
  if (!nrow(coverage)) stop("`coverage` has no rows.", call. = FALSE)
  cell <- coverage$cell
  if (anyNA(coverage$date)) stop("`coverage$date` must have no NA.", call. = FALSE)
  if (!all(is.finite(cell) & cell == round(cell) & cell >= 1 & cell <= terra::ncell(grid))) {
    stop(sprintf("`coverage$cell` must hold cells of `grid`, 1 to %s.",
                 fmt_n(terra::ncell(grid))), call. = FALSE)
  }
  if (!all(is.finite(coverage$p) & coverage$p >= 0 & coverage$p <= 1)) {
    stop("`coverage$p` must lie in 0 to 1, with no NA.", call. = FALSE)
  }
  repeats <- sum(duplicated(coverage[c("date", "cell")]))
  if (repeats) {
    stop(sprintf("`coverage` must have one row per date and cell; %s %s one.",
                 plural(repeats, "row"), agree(repeats, "repeats", "repeat")), call. = FALSE)
  }
  data.frame(date = coverage$date, cell = coverage$cell, p = coverage$p)
}

# ---- Domains ----
# Strata draw from the coverage rows of their bin's UTC date with p above 0,
# and strata of one date share them. Keys 1, 2, ... number the dates in the
# order of the strata; a stratum whose date has no such row has key NA and a
# reason. `cells` holds the rows of each key in turn.
sa_domains <- function(strata, coverage, grid) {
  bin_date <- format(as.Date(strata$bin_start, tz = "UTC"))
  have_date <- bin_date %in% format(coverage$date)
  cov <- coverage[coverage$p > 0, , drop = FALSE]
  xy <- terra::xyFromCell(grid, cov$cell)
  cov$x <- xy[, 1]
  cov$y <- xy[, 2]
  by_date <- split(seq_len(nrow(cov)), format(cov$date))
  dates <- unique(bin_date[bin_date %in% names(by_date)])
  rows <- unname(by_date[dates])
  key <- match(bin_date, dates)
  reason <- ifelse(!is.na(key), NA_character_,
                   ifelse(have_date, "no qualifying cell on the date", "no coverage on the date"))
  cells <- cov[unlist(rows), c("cell", "p", "x", "y")]
  cells$key <- rep(seq_along(rows), lengths(rows))
  list(key = key, reason = reason, cells = cells)
}

# ---- Draws ----
# Draws n rows of `cells` for each stratum, in proportion to p within the
# stratum's domain, in one pass: domain k occupies the interval (k - 1, k] of
# one cumulative scale, as in sample_locations(), and a uniform draw in
# (k - 1, k) falls in one of its rows. Rounding at an edge is clamped to the
# domain's own rows.
sa_draw <- function(cells, key, n) {
  w <- cells$p / stats::ave(cells$p, cells$key, FUN = sum)
  edge <- cells$key - 1 + stats::ave(w, cells$key, FUN = cumsum)
  last <- !duplicated(cells$key, fromLast = TRUE)
  edge[last] <- cells$key[last]
  first_row <- match(seq_len(max(cells$key, 0L)), cells$key)
  last_row <- which(last)
  k <- rep(key, each = n)
  u <- k - 1 + stats::runif(length(k))
  rows <- findInterval(u, edge) + 1L
  pmin(pmax(rows, first_row[k]), last_row[k])
}

sa_warn_dropped <- function(dropped, total, uniform_above) {
  if (!nrow(dropped)) return(invisible(NULL))
  no_date <- sum(dropped$reason == "no coverage on the date")
  rule <- if (is.null(uniform_above)) "p above 0" else
    paste("p at or above", format(uniform_above))
  msg <- paste("%s of %s strata (animal, bin and draw) have no cell to draw and are dropped:",
               "%s with no coverage on their date, %s with no cell of %s on it.",
               "attr(, \"dropped\") lists them.")
  warning(wrap_text(sprintf(msg, fmt_n(nrow(dropped)), fmt_n(total), fmt_n(no_date),
                            fmt_n(nrow(dropped) - no_date), rule)),
          call. = FALSE)
}
