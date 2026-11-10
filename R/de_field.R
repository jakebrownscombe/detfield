#' Detection probability from every water cell to every receiver
#'
#' Evaluates a detection model at the distance through water between each
#' receiver and each water cell, for one tag model and one date. The
#' result is the per-transmission detection probability that the location
#' likelihood uses, stored on the logit scale.
#'
#' A pair that no water path connects, or that [water_distances()] left out
#' beyond `max_distance`, has a logit of `-Inf`: the receiver cannot hear
#' that cell. A pair whose straight line meets land gets the model's
#' `blocked` offset added to its logit; see [de_model()]. With `blocked = 0`,
#' a message counts the pairs with land between them that the curve still
#' gives a detection probability of 0.05 or more.
#'
#' @param distances A [water_distances()] table.
#' @param model A [de_model()].
#' @param depth Receiver depth in metres, the `depth_m` the model reads: one
#'   number for every receiver, or a vector named by receiver key. By default
#'   the `deploy_z` column of the receivers given to [water_distances()], or
#'   `NA` when there was none.
#' @param tag_model The tag model, for a model that reads `tag_model`.
#' @param date One date, a `Date` or a POSIXct in UTC, at which the model's
#'   modulators are evaluated. Needed only when a modulator is a function.
#'
#' @return An object of class `de_field`, a list with:
#'   * `eta`: a matrix of logits, one row per water cell and one column per
#'     receiver. `field$p` returns the same matrix on the probability scale,
#'     computed when asked for.
#'   * `cells`: the terra cell numbers of the rows.
#'   * `receivers`: the receiver keys of the columns.
#'   * `key`: `"station_id"` or `"receiver_serial"`, the column the keys come
#'     from.
#'   * `lookup`: receiver serials, stations and deployment times when
#'     [water_distances()] had them, else `NULL`.
#'   * `receiver_xy`: the key and the position (`x`, `y`, in the CRS of the
#'     grid) used for each receiver.
#'   * `grid`, `model`, `distances`, and `inputs`, a list holding `depth`,
#'     `tag_model`, `date` and the modulator values used.
#'
#'   `as.data.frame()` gives one row per pair in `distances`, with columns
#'   `cell`, the receiver key, `distance_m`, `in_water`, `eta` and `p`.
#'   `field$p` and `field[["p"]]` give the probability matrix.
#'
#' @seealso [system_de()] for the probability that the array detects a
#'   transmission from each cell; [plot.de_field()].
#' @export
#' @examples
#' # The Stoney example: Stoney Lake (100 m cells, land NA) and 70 receivers.
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Pairs beyond 2 km are left out: there every receiver hears less than 1e-4.
#' d <- water_distances(grid, stoney_track$receivers, max_distance = 2000)
#' # The curve of the Stoney example; each receiver reads it at its own depth,
#' # deploy_z.
#' curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' field <- de_field(d, curve)
#' field
#' range(field$p[, "S005"])
#' head(as.data.frame(field))
de_field <- function(distances, model, depth = NULL, tag_model = NULL, date = NULL) {
  check_distances(distances)
  if (!inherits(model, "de_model")) stop("`model` must be a de_model.", call. = FALSE)
  key <- attr(distances, "key")
  rx <- attr(distances, "receivers")
  ids <- rx[[key]]
  if ("tag_model" %in% model$uses && !is_one_string(tag_model)) {
    stop("This detection model reads `tag_model`; give one tag model.", call. = FALSE)
  }
  depth_by_rx <- field_depth(depth, rx, ids, "depth_m" %in% model$uses)
  date <- check_date(date)
  mods <- eval_modulators(model$modulators, date, ids)
  col <- match(distances[[key]], ids)
  args <- c(list(distance_m = distances$distance_m, depth_m = depth_by_rx[col],
                 tag_model = if (is.null(tag_model)) NA_character_ else tag_model),
            lapply(mods, function(v) v[col]))
  eta_pairs <- do.call(model$predict, args)
  blocked <- if (is.null(model$blocked)) 0 else model$blocked
  if (blocked != 0) eta_pairs[!distances$in_water] <- eta_pairs[!distances$in_water] + blocked
  if (blocked == 0) note_shadowed(distances$in_water, eta_pairs)
  cells <- attr(distances, "cells")
  eta <- matrix(-Inf, length(cells), length(ids), dimnames = list(NULL, ids))
  eta[cbind(match(distances$cell, cells), col)] <- eta_pairs
  receiver_xy <- rx[, c(key, "x", "y")]
  rownames(receiver_xy) <- NULL
  structure(list(
    eta = eta, cells = cells, receivers = ids, key = key, lookup = attr(distances, "lookup"),
    receiver_xy = receiver_xy, grid = attr(distances, "grid"), model = model,
    distances = distances,
    inputs = list(depth = depth_by_rx, tag_model = tag_model, date = date, modulators = mods,
                  blocked = blocked)
  ), class = "de_field")
}

# With blocked = 0 the curve, usually fitted on clear lines of sight, is read
# on pairs with land between them. Reference tags can detect at a small
# share of the clear curve on such pairs, so a message gives the count.
note_shadowed <- function(in_water, eta) {
  heard <- eta >= stats::qlogis(0.05)
  shadowed <- !in_water & heard
  if (any(shadowed)) {
    message(wrap_text(sprintf(
      paste0("%s of %s receiver-cell pairs with p >= 0.05 have land on the ",
             "straight line between them; with blocked = 0 the curve is read on ",
             "them as if the line were clear. See `blocked` in ?de_model."),
      fmt_n(sum(shadowed)), fmt_n(sum(heard))
    )))
  }
  invisible(sum(shadowed))
}

check_distances <- function(distances) {
  need <- c("grid", "cells", "key", "receivers")
  has <- vapply(need, function(a) !is.null(attr(distances, a)), logical(1))
  if (!inherits(distances, "water_distances") || !all(has)) {
    stop("`distances` must be a table from water_distances(), with its attributes.",
         call. = FALSE)
  }
  invisible(TRUE)
}

field_depth <- function(depth, rx, ids, needed) {
  if (is.null(depth)) {
    depth <- if ("deploy_z" %in% names(rx)) rx$deploy_z else rep(NA_real_, length(ids))
    names(depth) <- ids
  } else if (!is.numeric(depth)) {
    stop("`depth` must be numeric.", call. = FALSE)
  } else if (length(depth) == 1 && is.null(names(depth))) {
    depth <- stats::setNames(rep(depth, length(ids)), ids)
  } else {
    missing <- setdiff(ids, names(depth))
    if (is.null(names(depth)) || length(missing)) {
      stop("`depth` must be one number or a vector named by receiver key, with every receiver.",
           call. = FALSE)
    }
    depth <- depth[ids]
  }
  if (needed && anyNA(depth)) {
    stop(sprintf("The detection model reads depth_m, but %s %s no depth.",
                 plural(sum(is.na(depth)), "receiver"), agree(sum(is.na(depth)), "has", "have")),
         call. = FALSE)
  }
  depth
}

check_date <- function(date) {
  if (is.null(date)) return(NULL)
  if (length(date) != 1 || is.na(date)) stop("`date` must be one date.", call. = FALSE)
  if (inherits(date, "POSIXct")) {
    if (!identical(attr(date, "tzone"), "UTC")) {
      stop("`date` must be in UTC; set it with `attr(date, \"tzone\") <- \"UTC\"`.",
           call. = FALSE)
    }
    return(date)
  }
  if (!inherits(date, "Date")) stop("`date` must be a Date or a POSIXct in UTC.", call. = FALSE)
  date
}

# One value per receiver for each modulator, in the order of `ids`.
eval_modulators <- function(modulators, date, ids) {
  lapply(stats::setNames(names(modulators), names(modulators)), function(nm) {
    m <- modulators[[nm]]
    if (is.function(m)) {
      if (is.null(date)) stop(sprintf("Modulator `%s` is a function of date; give `date`.", nm),
                              call. = FALSE)
      m <- m(date)
    }
    if (!is.null(names(m)) && length(m) > 1) {
      missing <- setdiff(ids, names(m))
      if (length(missing)) {
        stop(sprintf("Modulator `%s` has no value for %s %s.", nm,
                     agree(length(missing), "receiver", "receivers"),
                     paste(utils::head(missing, 5), collapse = ", ")), call. = FALSE)
      }
      return(unname(m[ids]))
    }
    if (length(m) != 1) {
      stop(sprintf("Modulator `%s` must give one value or one named value per receiver.", nm),
           call. = FALSE)
    }
    rep(m, length(ids))
  })
}

# A field built by hand may store `p` instead of `eta`; that `p` is returned.
field_element <- function(x, name) {
  x <- unclass(x)
  if (identical(name, "p") && is.null(x[["p"]]) && !is.null(x[["eta"]])) {
    return(stats::plogis(x[["eta"]]))
  }
  x[[name]]
}

#' @export
`$.de_field` <- function(x, name) field_element(x, name)

#' @export
`[[.de_field` <- function(x, i, ...) field_element(x, i)

# The formals repeat the generic's, row.names included.
# nolint start: object_name_linter.
#' @export
as.data.frame.de_field <- function(x, row.names = NULL, optional = FALSE, ...) {
  # nolint end
  d <- x$distances
  col <- match(d[[x$key]], x$receivers)
  eta <- x$eta[cbind(match(d$cell, x$cells), col)]
  out <- data.frame(cell = d$cell, key = d[[x$key]], distance_m = d$distance_m,
                    in_water = d$in_water, eta = eta, p = stats::plogis(eta),
                    stringsAsFactors = FALSE)
  names(out)[2] <- x$key
  out
}

#' @export
print.de_field <- function(x, ...) {
  cat_wrapped(sprintf("<de_field> %s x %s (by %s)\n", plural(length(x$cells), "water cell"),
                      plural(length(x$receivers), "receiver"), x$key))
  cat_wrapped(paste("  model:", x$model$source))
  if (!is.null(x$inputs$tag_model)) cat("  tag_model:", x$inputs$tag_model, "\n")
  if (!is.null(x$inputs$date)) cat("  date:", format(x$inputs$date), "\n")
  d <- x$distances
  cat_wrapped(sprintf("  blocked line of sight: %s of %s pairs, %s\n",
                      fmt_n(sum(!d$in_water)), fmt_n(nrow(d)),
                      blocked_label(x$inputs$blocked)))
  p <- x$p
  heard <- colSums(p >= 0.05)
  cat_wrapped(sprintf("  cells with p >= 0.05 per receiver: median %s, range %s to %s\n",
                      fmt_n(stats::median(heard)),
                      fmt_n(min(heard)), fmt_n(max(heard))))
  invisible(x)
}
