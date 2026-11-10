#' Detection probability of the whole array, per listening day
#'
#' For each water cell and each day, the probability that at least `k`
#' listening receivers detect one transmission from the cell, treating the
#' receivers as independent given the cell. With `k = 1` this is one minus
#' the product over listening receivers of one minus their detection
#' probability. The station effects in `sigma` are not integrated here: the
#' result is the array's coverage under the curve as given. Pedersen et al.
#' (2014) define the coverage of a grid cell for a sensor network the same
#' way, from the product over sensors of one minus each sensor's detection
#' probability.
#'
#' With `transmissions` above 1, the same probability for a bin in which the
#' tag sent that many transmissions: each receiver hears at least one of them
#' with probability `1 - (1 - p)^transmissions`, the transmissions being
#' independent given the cell, and at least `k` receivers must each hear one.
#' This is the probability that the bin is detected from the cell, which is
#' the weight [sample_available()] needs for availability in hourly bins: an
#' hour holds `3600 / mean delay` transmissions, 24 at the 150 s mean delay of
#' the Stoney tags. A spread of a receiver's detection probability between
#' hours (`sigma$station_hour`) makes the transmissions of one hour fail
#' together, and the probability lower than the value given here; it is not
#' integrated. `transmissions = 1` gives exactly the per-transmission result.
#'
#' Days that share the same set of listening receivers are computed once. A
#' day with fewer listening receivers than `k` gives 0 in every cell, and
#' every value lies in 0 to 1: the sum behind it is clamped, so rounding never
#' gives a probability below 0.
#'
#' **Days with no listening receiver.** A day on which every row has
#' `listening = FALSE` and a `status` of `"off"` or `"deaf"` gives 0 in every
#' cell and `n_listening` 0: the receivers' own logs show that none of them
#' could hear. Any other day without a listening row is left out of the
#' result: a `status` of `"not_downloaded"`, `"reference_silent"` or
#' `"silent_run"`, a missing `status` column, or `listening = NA`. On such a
#' day a receiver may have been listening with its data missing, and a 0
#' would misstate the coverage. A date with no row in `listening` is also
#' left out. On a day with at least one listening row, the rows that are not
#' listening are dropped whatever their `status`, so a `"not_downloaded"`
#' receiver is left out and the coverage comes from the receivers that
#' listened.
#'
#' @param field A [de_field()].
#' @param listening A data frame with one row per receiver and day: a
#'   `date` column (`Date`), and the field's key column (`station_id` or
#'   `receiver_serial`) or, for a field keyed by station, `receiver_serial`
#'   together with the deployment lookup of the field. A day on which a
#'   serial moves goes to its later deployment; a recovery day with no later
#'   deployment of that serial goes to the station it left. An optional logical
#'   `listening` column drops the rows that are not `TRUE`, and an optional
#'   `status` column, as `detqc::audit_listening()` writes it, says why; see
#'   Details for a day with no listening row. The list `list(deployments = , days = )`
#'   that `detqc::as_detfield_input()` returns as `listening` is also taken:
#'   its `days` table is used.
#' @param k The number of receivers that must detect the transmission; at
#'   most the number of receivers the field holds.
#' @param transmissions The number of transmissions the tag sends in a bin,
#'   one whole number of at least 1. The default, 1, gives the probability for
#'   one transmission.
#' @param format `"table"` for a data frame with columns `date`, `cell`,
#'   `p` and `n_listening`; `"raster"` for a [terra::SpatRaster] with one
#'   layer per day, named by date.
#'
#' @return A data frame or a SpatRaster, as `format` says. `p` is the
#'   probability for one transmission, or for a bin of `transmissions`
#'   transmissions. A receiver in `listening` that the field does not hold
#'   is an error on any row, listening or not, including the rows of a day
#'   that is left out. A field built for one date (`de_field(date = )`)
#'   scored against a listening table with no row on that date gives a
#'   warning, because the result takes its dates from `listening`.
#'
#' @references
#' Pedersen, M. W., Burgess, G. and Weng, K. C. (2014). A quantitative
#' approach to static sensor network design. *Methods in Ecology and
#' Evolution*, 5, 1043-1051. \doi{10.1111/2041-210X.12255}
#' @seealso [de_field()]; [sample_available()], which takes the per-bin table.
#' @export
#' @examples
#' # The Stoney example: Stoney Lake (100 m cells, land NA) and 70 receivers.
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' # Pairs beyond 2 km are left out: there every receiver hears less than 1e-4.
#' dist <- water_distances(grid, stoney_track$receivers, max_distance = 2000)
#' field <- de_field(dist, curve)
#' # All 70 stations listen on 1 June; S005 to S009 are silent on 2 June.
#' stations <- stoney_track$receivers$station_id
#' listening <- data.frame(
#'   station_id = c(stations, setdiff(stations, sprintf("S%03d", 5:9))),
#'   date = as.Date(rep(c("2024-06-01", "2024-06-02"), c(70, 65)))
#' )
#' cover <- system_de(field, listening)
#' tapply(cover$p, cover$date, mean)
#' system_de(field, listening, format = "raster")
#' # Per hour: 24 transmissions at the tags' 150 s mean delay
#' hourly <- system_de(field, listening, transmissions = 24)
#' tapply(hourly$p >= 0.95, hourly$date, mean)
system_de <- function(field, listening, k = 1, transmissions = 1,
                      format = c("table", "raster")) {
  if (!inherits(field, "de_field")) stop("`field` must be a de_field.", call. = FALSE)
  format <- match.arg(format)
  check_system_de_counts(k, transmissions, length(field$receivers))
  read <- listening_field_columns(field, listening_days(listening))
  col <- read$cols
  days <- sort(unique(c(col$date, read$deaf_days)))
  check_field_date(field, days)
  # A day with no listening receiver has the empty set.
  sets <- vapply(days, function(d) paste(sort(unique(col$col[col$date == d])), collapse = ","), "")
  unique_sets <- unique(sets)
  p <- field$p
  # A receiver hears a bin when it hears at least one of its transmissions.
  # With one transmission p is used as it is, so the result is unchanged.
  if (transmissions > 1) p <- 1 - (1 - p)^transmissions
  cover <- vapply(unique_sets, function(s) {
    at_least_k(p[, set_columns(s), drop = FALSE], k)
  }, numeric(length(field$cells)))
  cover <- matrix(cover, nrow = length(field$cells))[, match(sets, unique_sets), drop = FALSE]
  n_listening <- vapply(sets, function(s) length(set_columns(s)), integer(1), USE.NAMES = FALSE)
  if (format == "raster") {
    return(cells_to_raster(field$grid, field$cells, cover, names = format(days)))
  }
  data.frame(date = rep(days, each = length(field$cells)), cell = field$cells,
             p = as.vector(cover), n_listening = rep(n_listening, each = length(field$cells)))
}

# `k` and `transmissions` of system_de(), against the field's `n_receivers`.
check_system_de_counts <- function(k, transmissions, n_receivers) {
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != round(k)) {
    stop("`k` must be one whole number of at least 1.", call. = FALSE)
  }
  if (k > n_receivers) {
    stop(sprintf("`k` is %s, but the field holds %s, so no cell can be heard by %s.",
                 fmt_n(as.integer(k)), plural(n_receivers, "receiver"), fmt_n(as.integer(k))),
         call. = FALSE)
  }
  if (!ls_whole(transmissions) || transmissions < 1) {
    stop("`transmissions` must be one whole number of at least 1.", call. = FALSE)
  }
  invisible(NULL)
}

# A field built for one date and scored against listening rows of other dates
# would take their dates silently.
check_field_date <- function(field, days) {
  built_for <- field$inputs$date
  if (is.null(built_for)) return(invisible(NULL))
  built_for <- as.Date(built_for, tz = "UTC")
  if (!built_for %in% days) {
    warning(sprintf(paste0("The field was built for %s, but `listening` has no row on that ",
                           "date, so the result is dated %s. Give the listening rows of %s."),
                    format(built_for), paste(format(days), collapse = ", "), format(built_for)),
            call. = FALSE)
  }
  invisible(NULL)
}

# The hand-off of detqc::as_detfield_input() carries listening as a list of
# deployments and days, which bin_detections() takes whole; here the days
# table is the one that counts.
listening_days <- function(listening) {
  if (is.list(listening) && !is.data.frame(listening) && is.data.frame(listening$days)) {
    return(listening$days)
  }
  listening
}

set_columns <- function(s) {
  if (nzchar(s)) as.integer(strsplit(s, ",", fixed = TRUE)[[1]]) else integer(0)
}

# Statuses of detqc::audit_listening() that come from the receiver's own log
# and say it could not hear. Every other not-listening status is inferred
# from missing data, so the receiver may have been listening.
deaf_statuses <- c("off", "deaf")

# Rows of `listening` as field columns and dates (`cols`), and the dates on
# which no receiver listened and every row is off or deaf (`deaf_days`).
listening_field_columns <- function(field, listening) {
  if (!is.data.frame(listening) || !inherits(listening[["date"]], "Date")) {
    stop("`listening` must be a data frame with a `date` column of class Date.", call. = FALSE)
  }
  listening <- as.data.frame(listening)
  deaf_days <- listening$date[0]
  if ("listening" %in% names(listening)) {
    if (!is.logical(listening$listening)) stop("`listening$listening` must be logical.",
                                               call. = FALSE)
    deaf_days <- deaf_only_days(listening)
    # A receiver the field does not hold is a wrong input whether or not it
    # listened that day.
    field_columns(field, listening)
    listening <- listening[listening$listening %in% TRUE, , drop = FALSE]
  }
  if (!nrow(listening) && !length(deaf_days)) {
    stop("`listening` has no listening rows and no day on which every receiver is off or deaf.",
         call. = FALSE)
  }
  if (!nrow(listening)) return(list(cols = data.frame(date = deaf_days[0], col = integer(0)),
                                    deaf_days = deaf_days))
  list(cols = data.frame(date = listening$date, col = field_columns(field, listening)),
       deaf_days = deaf_days)
}

# The field column of each row of `listening`; a receiver the field does not
# hold is an error.
field_columns <- function(field, listening) {
  if (!nrow(listening)) return(integer(0))
  id <- listening_ids(field, listening)
  col <- match(id, field$receivers)
  if (anyNA(col)) {
    stop(sprintf("%s of `listening` %s the field does not hold: %s.",
                 plural(sum(is.na(col)), "row"),
                 agree(sum(is.na(col)), "names a receiver", "name receivers"),
                 paste(utils::head(unique(id[is.na(col)]), 5), collapse = ", ")),
         call. = FALSE)
  }
  col
}

# The field's key of each row of `listening`: its own key column, or the
# station of its serial through the field's deployment lookup.
listening_ids <- function(field, listening) {
  if (field$key %in% names(listening)) return(listening[[field$key]])
  if ("receiver_serial" %in% names(listening) && !is.null(field$lookup)) {
    return(stations_of_serials(field$lookup, listening$receiver_serial, listening$date))
  }
  stop(sprintf("`listening` needs a `%s` column.", field$key), call. = FALSE)
}

# Dates with no listening row on which every row is FALSE with an off or deaf
# status. Without a status column no day qualifies.
deaf_only_days <- function(listening) {
  if (!"status" %in% names(listening)) return(listening$date[0])
  deaf <- listening$listening %in% FALSE & listening$status %in% deaf_statuses
  all_deaf <- tapply(deaf, as.character(listening$date), all)
  sort(as.Date(names(all_deaf)[all_deaf]))
}

stations_of_serials <- function(lookup, serial, date) {
  out <- serial_to_station(lookup, serial, date)
  lost <- is.na(out)
  if (any(lost)) {
    pairs <- unique(paste(serial[lost], "on", format(date[lost])))
    stop(sprintf("%s of `listening` %s a receiver_serial with no deployment on that date: %s.",
                 plural(sum(lost), "row"), agree(sum(lost), "names", "name"),
                 paste(utils::head(pairs, 5), collapse = ", ")),
         call. = FALSE)
  }
  out
}

# The station each serial sat at on each date. CONTRACT.md, section 1, gives
# a deployment the times deploy_datetime <= t < recover_datetime, so a day
# holds every deployment that listened during it. When a serial moves on a
# day, the day goes to the later deployment; a recovery day with no later
# deployment of that serial belongs to the deployment that ended.
serial_to_station <- function(lookup, serial, date) {
  if (!all(c("deploy_datetime", "recover_datetime") %in% names(lookup))) {
    one <- tapply(lookup$station_id, lookup$receiver_serial, function(s) length(unique(s)))
    if (any(one > 1)) {
      stop("A receiver_serial sits at more than one station and the lookup has no deployment ",
           "times; give `listening` the station_id column.", call. = FALSE)
    }
    return(lookup$station_id[match(serial, lookup$receiver_serial)])
  }
  start <- as.Date(lookup$deploy_datetime, tz = "UTC")
  end <- as.Date(lookup$recover_datetime, tz = "UTC")
  end[is.na(end)] <- as.Date(Inf)
  out <- rep(NA_character_, length(serial))
  for (i in seq_len(nrow(lookup))) {
    at <- which(serial == lookup$receiver_serial[i] & date == end[i] & is.na(out))
    out[at] <- lookup$station_id[i]
  }
  for (i in seq_len(nrow(lookup))) {
    at <- which(serial == lookup$receiver_serial[i] & date >= start[i] & date < end[i])
    out[at] <- lookup$station_id[i]
  }
  out
}

# P(at least k successes) for independent Bernoulli columns, row by row. The
# vector q holds P(exactly j successes so far) for j = 0, ..., k - 1. One
# minus their sum rounds below 0 where the array cannot hear a cell, as at
# k = 3 on a large grid, so it is clamped; with
# fewer columns than k the answer is exactly 0.
at_least_k <- function(p, k) {
  if (ncol(p) < k) return(rep(0, nrow(p)))
  q <- matrix(0, nrow(p), k)
  q[, 1] <- 1
  for (s in seq_len(ncol(p))) {
    ps <- p[, s]
    if (k > 1) q[, k:2] <- q[, k:2] * (1 - ps) + q[, (k - 1):1] * ps
    q[, 1] <- q[, 1] * (1 - ps)
  }
  pmin(1, pmax(0, 1 - rowSums(q)))
}
