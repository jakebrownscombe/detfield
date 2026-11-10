# ---- Count a range test ----

#' Count range-test detections per receiver, range-test tag and hour
#'
#' @description
#' Turns the detections of a range test into the counts a range model is
#' fitted to. For every receiver and range-test tag deployed together for a
#' whole bin, the result has one row: the transmissions the tag sent
#' (`sent`), the detections the receiver logged (`heard`, zero when it logged
#' none) and the distance between them (`distance_m`).
#'
#' @details
#' **Every pair, every bin.** A bin in which a receiver heard nothing from a
#' tag is kept as a row with `heard = 0`. These rows matter most at long
#' range, where silence is the common outcome: without them the far end of
#' the curve rests only on the bins that happened to be heard, and the fitted
#' range comes out too long. A receiver and a tag form a pair in a bin when
#' both deployments cover the whole bin.
#'
#' **Transmissions.** `sent` is the bin length over the tag's mean delay,
#' rounded, as in [bin_detections()]. Where a receiver logged more
#' detections than that in a bin, `sent` is raised to the count, and a
#' message gives the number of rows raised.
#'
#' **Distance.** From the deployment positions of receiver and tag
#' (`deploy_lat` and `deploy_lon`, WGS 84), measured on an azimuthal
#' equidistant projection centred on the test and rounded to 0.1 m.
#'
#' **Receiver depth.** When `receivers` has `deploy_z`, the result carries
#' it as `depth_m`, the name [de_model()] passes to a range model.
#'
#' A deployment with no `recover_datetime` ends at the last detection in
#' `detections`. Detections of other transmitters or receivers, and
#' detections in bins that a deployment covers only in part, are not counted;
#' a message gives their number.
#'
#' @param detections Detections of the range test, one row each:
#'   `datetime` (POSIXct, UTC), `receiver_serial` and `transmitter`.
#' @param receivers Receiver deployments of the range test, one row each:
#'   `receiver_serial`, `station_id`, `deploy_datetime`, `recover_datetime`,
#'   `deploy_lat`, `deploy_lon` and, optionally, `deploy_z`.
#' @param tags Range-test tags, one row per transmitter: `transmitter`,
#'   `deploy_datetime`, `recover_datetime`, `deploy_lat`, `deploy_lon`, and
#'   either `delay` or both `min_delay` and `max_delay`, in seconds. Other
#'   columns, such as `tag_model`, are carried into the result.
#' @param bin Bin length: `"hour"`, `"day"` or a number of seconds.
#' @return A data frame with one row per receiver, transmitter and bin:
#'   `datetime` (start of the bin, UTC), `receiver_serial`, `station_id`,
#'   `depth_m` when `receivers` has `deploy_z`, `transmitter`, the other
#'   columns of `tags`, `distance_m`, `sent` and `heard`.
#' @seealso [de_model()], which wraps the range model fitted to these counts,
#'   and `vignette("range-model")`, which fits one with tag model, receiver
#'   depth and wind.
#' @export
#' @examples
#' # The Stoney range test: two lines, at 1 m and 25 m, each of three receivers
#' # and five range-test tags, for two days. Serial numbers are read as text.
#' in_file <- function(name) system.file("extdata", name, package = "detfield")
#' read_utc <- function(name, cols, serial = TRUE) {
#'   x <- utils::read.csv(in_file(name),
#'                        colClasses = if (serial) {
#'                          c(receiver_serial = "character")
#'                        } else NA)
#'   for (col in cols) x[[col]] <- as.POSIXct(x[[col]], tz = "UTC")
#'   x
#' }
#' times <- c("deploy_datetime", "recover_datetime")
#' detections <- read_utc("stoney_range_test_detections.csv", "datetime")
#' receivers <- read_utc("stoney_range_test_receivers.csv", times)
#' tags <- read_utc("stoney_range_test_tags.csv", times, serial = FALSE)
#' # The messages count the detections in part-hours of a deployment, which are
#' # not counted, and the rows in which a receiver heard more transmissions than
#' # the mean delay gives, where `sent` is raised.
#' counts <- range_test_counts(detections, receivers, tags)
#' head(counts)
#' # Hours in which a receiver heard nothing from a tag stay in, with heard = 0.
#' c(rows = nrow(counts), heard_nothing = sum(counts$heard == 0))
range_test_counts <- function(detections, receivers, tags, bin = "hour") {
  bin_s <- bd_bin_seconds(bin)
  rt_check(detections, receivers, tags)
  delay <- rt_delays(tags)
  ends_at <- max(detections$datetime, na.rm = TRUE)
  pairs <- rt_pairs(receivers, tags, ends_at)
  rows <- rt_bins(pairs, bin_s)
  if (!nrow(rows)) {
    stop("No receiver and tag were deployed together for a whole bin.", call. = FALSE)
  }
  start <- floor(as.numeric(detections$datetime) / bin_s) * bin_s
  hit <- match(paste(detections$receiver_serial, detections$transmitter, start),
               paste(rows$receiver_serial, rows$transmitter, as.numeric(rows$datetime)))
  if (anyNA(hit)) {
    message(wrap_text(sprintf(
      paste0("%s of %s detections are not counted: other transmitters or ",
             "receivers, or bins a deployment covers only in part."),
      fmt_n(sum(is.na(hit))), fmt_n(length(hit))
    )))
  }
  rows$heard <- tabulate(hit[!is.na(hit)], nbins = nrow(rows))
  rows$sent <- as.integer(round(bin_s / delay[rows$transmitter]))
  rt_finish(rows)
}

rt_check <- function(detections, receivers, tags) {
  bd_check_columns(detections, c(datetime = "POSIXct", receiver_serial = "character",
                                 transmitter = "character"), "detections")
  bd_check_utc(detections, "datetime", "detections")
  placed <- c(deploy_datetime = "POSIXct", recover_datetime = "POSIXct",
              deploy_lat = "numeric", deploy_lon = "numeric")
  bd_check_columns(receivers, c(receiver_serial = "character", station_id = "character",
                                placed), "receivers")
  bd_check_columns(tags, c(transmitter = "character", placed), "tags")
  rt_check_placed(receivers, "receivers")
  rt_check_placed(tags, "tags")
  if (!nrow(detections)) stop("`detections` has no rows.", call. = FALSE)
  if (anyDuplicated(tags$transmitter)) {
    stop("`tags` must have one row per transmitter.", call. = FALSE)
  }
  invisible(TRUE)
}

rt_check_placed <- function(tab, what) {
  for (col in c("deploy_datetime", "recover_datetime")) bd_check_utc(tab, col, what)
  if (anyNA(tab$deploy_datetime) || anyNA(tab$deploy_lat) || anyNA(tab$deploy_lon)) {
    stop(sprintf("`%s` needs `deploy_datetime`, `deploy_lat` and `deploy_lon` in every row.",
                 what), call. = FALSE)
  }
}

# Mean delay per transmitter, named by transmitter.
rt_delays <- function(tags) {
  keep <- intersect(names(tags), c("transmitter", "delay", "min_delay", "max_delay"))
  tab <- bd_delay_table(tags[, keep, drop = FALSE])
  missing <- setdiff(tags$transmitter, tab$transmitter)
  if (length(missing)) {
    stop(sprintf("`tags` gives no delay for %s.", bd_first(missing)), call. = FALSE)
  }
  stats::setNames(tab$delay, tab$transmitter)
}

# One row per receiver deployment and tag: the time both were in the water
# and the distance between them.
rt_pairs <- function(receivers, tags, ends_at) {
  rx <- receivers
  rx$recover_datetime[is.na(rx$recover_datetime)] <- ends_at
  tg <- tags
  tg$recover_datetime[is.na(tg$recover_datetime)] <- ends_at
  ij <- expand.grid(r = seq_len(nrow(rx)), t = seq_len(nrow(tg)))
  pairs <- data.frame(
    receiver_serial = rx$receiver_serial[ij$r], station_id = rx$station_id[ij$r],
    transmitter = tg$transmitter[ij$t],
    from = pmax(as.numeric(rx$deploy_datetime[ij$r]), as.numeric(tg$deploy_datetime[ij$t])),
    to = pmin(as.numeric(rx$recover_datetime[ij$r]), as.numeric(tg$recover_datetime[ij$t])),
    distance_m = rt_distance(rx, tg)[cbind(ij$r, ij$t)],
    stringsAsFactors = FALSE
  )
  if ("deploy_z" %in% names(rx)) pairs$depth_m <- rx$deploy_z[ij$r]
  extra <- setdiff(names(tg), c("transmitter", "deploy_datetime", "recover_datetime",
                                "deploy_lat", "deploy_lon", "delay", "min_delay", "max_delay"))
  for (col in extra) pairs[[col]] <- tg[[col]][ij$t]
  attr(pairs, "extra") <- extra
  pairs
}

# Distances in metres between every receiver (rows) and tag (columns).
rt_distance <- function(rx, tg) {
  lon <- c(rx$deploy_lon, tg$deploy_lon)
  lat <- c(rx$deploy_lat, tg$deploy_lat)
  local <- sprintf("+proj=aeqd +lat_0=%.8f +lon_0=%.8f +datum=WGS84 +units=m",
                   mean(range(lat)), mean(range(lon)))
  pts <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  xy <- sf::st_coordinates(sf::st_transform(pts, local))
  n_rx <- nrow(rx)
  d <- sqrt(outer(xy[seq_len(n_rx), 1], xy[-seq_len(n_rx), 1], "-")^2 +
              outer(xy[seq_len(n_rx), 2], xy[-seq_len(n_rx), 2], "-")^2)
  round(d, 1)
}

# Every whole bin of each pair.
rt_bins <- function(pairs, bin_s) {
  first <- ceiling(pairs$from / bin_s) * bin_s
  last <- floor(pairs$to / bin_s) * bin_s - bin_s
  n <- pmax(0, (last - first) / bin_s + 1)
  idx <- rep(seq_len(nrow(pairs)), n)
  rows <- pairs[idx, setdiff(names(pairs), c("from", "to")), drop = FALSE]
  offset <- sequence(n) - 1
  rows$datetime <- as.POSIXct(first[idx] + offset * bin_s, origin = "1970-01-01", tz = "UTC")
  attr(rows, "extra") <- attr(pairs, "extra")
  rows
}

rt_finish <- function(rows) {
  over <- rows$heard > rows$sent
  if (any(over)) {
    message(wrap_text(sprintf(
      paste0("%s of %s rows had more detections than the delay allows ",
             "transmissions; `sent` was raised to the count in those rows."),
      fmt_n(sum(over)), fmt_n(nrow(rows))
    )))
    rows$sent[over] <- rows$heard[over]
  }
  cols <- c("datetime", "receiver_serial", "station_id", intersect("depth_m", names(rows)),
            "transmitter", attr(rows, "extra"), "distance_m", "sent", "heard")
  out <- rows[order(rows$datetime, rows$receiver_serial, rows$transmitter), cols]
  rownames(out) <- NULL
  out
}
