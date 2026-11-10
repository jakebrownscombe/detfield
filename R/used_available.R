#' Stack used and available points for a habitat-selection model
#'
#' @description
#' Stacks the used points of [sample_locations()] and the available points of
#' [sample_available()] into the long table a conditional logistic regression
#' takes: one stratum per animal, bin and draw, holding its one used point
#' (`case = 1`) and its available points (`case = 0`), with one column per
#' layer of `layers` read at each point's cell. This is the use-availability
#' design, which compares the resources an animal used with those available
#' to it (Johnson 1980).
#'
#' @details
#' Each point is read at the cell of `layers` that holds its `x` and `y`. A
#' `cell` column in `presence` or `available` must name that same cell, or the
#' points were drawn on another grid and the call stops.
#'
#' A row with a missing layer value is dropped, and the attribute `dropped`
#' counts the rows by reason. A stratum whose used point has a missing value is
#' dropped whole, as is a stratum left with no available point; a used point
#' whose stratum has no available points in `available` (a stratum that
#' [sample_available()] dropped) is dropped too. Every kept stratum therefore
#' holds one used point and at least one available point, and without missing
#' values exactly `attr(available, "n")` of them. A warning gives the counts
#' whenever a row is dropped.
#'
#' The function fits no model. A conditional logit per animal on one draw's
#' strata, such as
#' `survival::clogit(case ~ depth_m + strata(stratum), data = d[d$draw == 1, ])`,
#' gives one imputation; repeating it for each draw and pooling with Rubin's
#' rules carries the position uncertainty of the surfaces. Fitting every draw
#' at once treats one position as many and narrows the interval. Columns of
#' `presence` other than the ones named below are not carried.
#'
#' @param presence The used points: a data frame with `animal_id`,
#'   `bin_start` (POSIXct, UTC), `draw`, `x` and `y`, one row per animal, bin
#'   and draw, such as the output of [sample_locations()].
#' @param available The available points from [sample_available()] for the
#'   same `presence`.
#' @param layers A `SpatRaster` of covariates on the grid of the points, one
#'   layer per covariate, each with a name that becomes a column.
#' @return A data frame with one row per point: `stratum` (an integer per
#'   animal, bin and draw), `case` (1 used, 0 available), `animal_id`,
#'   `bin_start`, `draw`, `cell`, `x`, `y`, and one column per layer, ordered
#'   by stratum with the used point first. Attributes: `"estimand"` and `"n"`
#'   from `available`, `"epsg"`, and `"dropped"`, a data frame of rows dropped
#'   by `case` and `reason`.
#' @references
#' Johnson, D. H. (1980). The comparison of usage and availability
#' measurements for evaluating resource preference. *Ecology*, 61, 65-71.
#' \doi{10.2307/1937156}
#' @seealso [sample_locations()], [sample_available()]
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
#' presence <- sample_locations(location_surface(counts, field), n = 5, seed = 1)
#' listening <- expand.grid(station_id = stoney_track$receivers$station_id,
#'                          date = unique(as.Date(presence$bin_start, tz = "UTC")),
#'                          stringsAsFactors = FALSE)
#' hourly <- system_de(field, listening, transmissions = 24)
#' available <- sample_available(presence, hourly, grid, n = 20, seed = 1)
#' # Bottom depth in metres, positive down, as the covariate
#' depth <- abs(grid)
#' names(depth) <- "depth_m"
#' stacked <- used_available(presence, available, depth)
#' head(stacked)
#' table(case = stacked$case)
used_available <- function(presence, available, layers) {
  layers <- ua_layers(layers)
  carried <- list(estimand = attr(available, "estimand", exact = TRUE),
                  n = attr(available, "n", exact = TRUE))
  presence <- sa_check_points(presence, layers, "presence")
  available <- ua_check_available(available, layers)
  used <- ua_read(presence, layers, 1L)
  avail <- ua_read(available, layers, 0L)
  key_used <- ua_key(used)
  key_avail <- ua_key(avail)
  if (!all(key_avail %in% key_used)) {
    n_orphan <- sum(!key_avail %in% key_used)
    stop(wrap_text(sprintf(paste("%s of `available` %s to no row of `presence`;",
                                 "give the presence points `available` was drawn for."),
                           plural(n_orphan, "row"), agree(n_orphan, "belongs", "belong"))),
         call. = FALSE)
  }
  vars <- names(layers)
  used_na <- !stats::complete.cases(used[vars])
  avail_na <- !stats::complete.cases(avail[vars])
  no_avail <- !key_used %in% key_avail
  lost_used <- key_avail %in% key_used[used_na]
  keep_avail <- !avail_na & !lost_used
  left <- !key_used %in% key_avail[keep_avail]
  empty <- left & !used_na & !no_avail
  keep_used <- !(used_na | no_avail | empty)
  reasons <- c("stratum has no available points", "missing layer value",
               "every available point has a missing layer value", "missing layer value",
               "used point has a missing layer value")
  counts <- c(sum(no_avail), sum(used_na & !no_avail), sum(empty),
              sum(avail_na & !lost_used), sum(lost_used))
  dropped <- data.frame(case = c(1L, 1L, 1L, 0L, 0L), reason = reasons, rows = counts,
                        stringsAsFactors = FALSE)
  ua_warn(dropped)
  out <- rbind(used[keep_used, ], avail[keep_avail, ])
  out <- out[order(out$animal_id, out$bin_start, out$draw, -out$case), ]
  out$stratum <- match(ua_key(out), unique(ua_key(out)))
  out <- out[c("stratum", "case", "animal_id", "bin_start", "draw", "cell", "x", "y", vars)]
  rownames(out) <- NULL
  attr(out, "estimand") <- carried$estimand
  attr(out, "n") <- carried$n
  attr(out, "epsg") <- grid_epsg(layers)
  attr(out, "dropped") <- dropped
  out
}

# ---- Helpers ----
ua_reserved <- c("stratum", "case", "animal_id", "bin_start", "draw", "cell", "x", "y")

ua_layers <- function(layers) {
  if (inherits(layers, "PackedSpatRaster")) layers <- terra::rast(layers)
  if (!inherits(layers, "SpatRaster")) {
    stop("`layers` must be a SpatRaster of covariates.", call. = FALSE)
  }
  vars <- names(layers)
  bad <- anyNA(vars) || any(!nzchar(vars)) || anyDuplicated(vars) > 0 ||
    any(vars != make.names(vars))
  if (bad) {
    stop("The layers of `layers` must have unique names that are valid column names.",
         call. = FALSE)
  }
  if (any(vars %in% ua_reserved)) {
    stop(wrap_text(sprintf("`layers` has a layer named %s, a column the result already holds.",
                           paste0("`", vars[vars %in% ua_reserved], "`", collapse = ", "))),
         call. = FALSE)
  }
  layers
}

# The available points: the columns of sample_available(), one stratum per
# animal, bin and draw.
ua_check_available <- function(available, grid) {
  bd_check_columns(available, c(animal_id = "character", bin_start = "POSIXct",
                                draw = "numeric", cell = "numeric", x = "numeric",
                                y = "numeric"), "available")
  if (!nrow(available)) stop("`available` has no rows.", call. = FALSE)
  bd_check_utc(available, "bin_start", "available")
  if (!all(is.finite(available$x) & is.finite(available$y))) {
    stop("`available$x` and `available$y` must be finite.", call. = FALSE)
  }
  sa_check_epsg(attr(available, "epsg", exact = TRUE), grid, "available")
  out <- as.data.frame(available)[c("animal_id", "bin_start", "draw", "cell", "x", "y")]
  sa_check_cells(out, grid, "available")
  out
}

# Points with their cell on `layers` and the layer values there.
ua_read <- function(points, layers, case) {
  # CONTRACT.md types `cell` as integer; cellFromXY() returns a double.
  cell <- as.integer(terra::cellFromXY(layers, cbind(points$x, points$y)))
  # Read every row at once, a point outside the layers at cell 1, then blank
  # those rows, so a categorical layer keeps its levels.
  outside <- is.na(cell)
  values <- terra::extract(layers, ifelse(outside, 1, cell))
  values[outside, ] <- NA
  data.frame(points[c("animal_id", "bin_start", "draw")], cell = cell, x = points$x,
             y = points$y, case = case, values, stringsAsFactors = FALSE,
             check.names = FALSE)
}

ua_key <- function(points) {
  paste(points$animal_id, as.numeric(points$bin_start), points$draw, sep = "\r")
}

ua_warn <- function(dropped) {
  if (!sum(dropped$rows)) return(invisible(NULL))
  used <- sum(dropped$rows[dropped$case == 1L])
  avail <- sum(dropped$rows[dropped$case == 0L])
  msg <- paste("%s used and %s available rows dropped: %s with no available",
               "points, %s with a missing layer value, %s whose available points all have",
               "one; %s with a missing layer value and %s in a stratum",
               "whose used point has one. attr(, \"dropped\") counts them.")
  r <- dropped$rows
  warning(wrap_text(sprintf(msg, fmt_n(used), fmt_n(avail), plural(r[1], "used point"),
                            fmt_n(r[2]), fmt_n(r[3]), plural(r[4], "available point"),
                            fmt_n(r[5]))),
          call. = FALSE)
}
