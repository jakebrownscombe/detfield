#' Count detections per animal, time bin and receiver
#'
#' @description
#' Turns a detection table into the counts that [location_surface()] reads.
#' For each animal and each time bin in which the animal was detected, the
#' result has one row per receiver listening in the bin: the detections that
#' receiver logged (`y`, zero for a silent receiver) and the transmissions the
#' animal sent in the bin (`n`).
#'
#' @details
#' **Transmissions.** `n` is the bin length divided by the tag's mean delay
#' between transmissions, rounded to a whole number. The delay is the mean
#' time from the start of one transmission to the start of the next. Where a
#' tag's nominal delay leaves out the transmission itself, a few seconds,
#' add it: otherwise `n` is a few per cent too high (4% for a 120 s delay and
#' a 5 s transmission), and silent receivers count as slightly stronger
#' evidence than they are. Check how the manufacturer defines the delay. When an animal carries
#' more than one transmitter, `codes` says how they share the tag's
#' transmissions:
#' * `"shared"`: the tag sends one transmission per delay and the codes take
#'   turns, so `n` is the bin length over the mean of the codes' delays. Give
#'   each code the tag's delay. In detections of V13AP and V9AP tags most
#'   consecutive detections switch code, so their codes take turns. For V16AT
#'   tags the evidence is mixed and the rule is unresolved; check them with
#'   the rate check below.
#' * `"independent"`: each code keeps its own schedule, so `n` sums the
#'   codes' rates.
#'
#' There is no default for animals with more than one transmitter, because
#' the wrong rule doubles or halves `n`. Give `codes`, or a `codes` column in
#' `tag_delay` to set it per transmitter.
#'
#' **The floor.** A receiver cannot log more transmissions than were sent, so
#' where a receiver logged more detections in a bin than `n`, `n` is raised to
#' that count. The number of bins raised is given in a message and in the
#' attribute `"n_raised"`.
#'
#' **Rate check.** For each animal with at least 100 detected bins, the
#' attribute `"rate_check"` gives the 99th percentile of the busiest
#' receiver's count per bin over `n`. With the right delay and rule it lies
#' near 1; near 0.5 means `n` is doubled, as when the codes of a tag that
#' takes turns are summed. A warning names animals below 0.6 or above 1.5.
#'
#' **Listening.** `listening` takes one of three forms.
#' * Receiver-days: one row per receiver and UTC date on which it listened,
#'   with the columns `receiver_serial` and `date` (class Date), and optionally
#'   `station_id` and a logical column `listening` (rows with `FALSE` are not
#'   used). A receiver listens in every bin of a listed date, so with this form
#'   a bin must divide a day. The form cannot see the hours of a deployment or
#'   recovery day when the receiver was out of the water.
#' * Deployments: one row per receiver deployment, with
#'   `receiver_serial`, `station_id`, `deploy_datetime` and `recover_datetime`
#'   (`NA` while deployed). A receiver listens in a bin when its deployment
#'   covers at least `min_cover` of the bin.
#' * Both, as `list(deployments = , days = )`: a receiver listens in a bin
#'   when its deployment covers at least `min_cover` of the bin and the bin's
#'   UTC date is a listed day. Use this form when receiver-days mark days that
#'   failed a check, so that neither the failed days nor the hours out of the
#'   water count as silence.
#'
#' Receiver-bins left out for covering less than `min_cover` are counted in
#' the attribute `"cover_dropped"`. At a station, the cover of all its
#' receivers in the bin is summed (at most 1), so a swap within a bin does not
#' drop the station.
#'
#' **Receiver swaps.** When the listening table carries `station_id` and two
#' receivers at one station both listen in a bin, as on the day one replaces
#' the other, the bin keeps one row for the station: the receiver deployed
#' later, with the detections of both. Without this the station would count
#' twice as a silent receiver. The number of merges is in the attribute
#' `"swaps_merged"`.
#'
#' **Detections outside the listening table.** A detection at a receiver not
#' listening in its bin keeps a row with `listening = FALSE`, and a warning
#' gives the count. [location_surface()] leaves such rows out; add the
#' receiver to the listening table, or lower `min_cover`, to use them.
#'
#' **Deployments with no serial.** Some deployment sheets, such as OTN's
#' simplified ones, give the station but not the receiver's serial. Rows of
#' `listening` with no `receiver_serial` but a `station_id` are then matched
#' to detections by station, when `detections` has a `station_id` column:
#' both take the key `"station:<station_id>"` in `receiver_serial`, and a
#' message gives the count. Without that column a warning says that those
#' rows cannot be matched. Use a field keyed by `station_id` with such
#' counts.
#'
#' Bins are aligned to UTC: hour bins start on the hour and day bins at
#' midnight UTC. Only bins with a detection are returned, so the surfaces
#' describe detected bins.
#'
#' @param detections Detections, one row each, with the columns `datetime`
#'   (POSIXct, UTC), `receiver_serial` (text), `transmitter` (text, the
#'   tag's code) and `animal_id` (text, the animal that carried the tag), and
#'   optionally `station_id`. Rows with no `animal_id`,
#'   such as reference tags, are not counted.
#' @param listening Which receivers listened when: receiver-days,
#'   deployments, or a list of both (see Details).
#' @param tag_delay Mean delay between transmissions, in seconds (see
#'   Details): one number for every transmitter, two numbers giving the
#'   shortest and longest delay (their mean is used), or a data frame with the
#'   column
#'   `transmitter` and either `delay` or both `min_delay` and `max_delay`
#'   (the mean of the two is used). An optional `animal_id` column lists
#'   every code an animal carries, including codes never detected, and an
#'   optional `codes` column sets the rule per transmitter.
#' @param bin Bin length: `"hour"`, `"day"` or a number of seconds.
#' @param codes How the transmitters of one tag share its transmissions:
#'   `"shared"` or `"independent"` (see Details). Needed when an animal
#'   carries more than one transmitter and `tag_delay` has no `codes` column.
#' @param min_cover Share of a bin, from 0 to 1, that a deployment must cover
#'   for its receiver to listen in the bin. Not used with receiver-days alone.
#' @return A data frame with one row per animal, bin and listening receiver:
#'   `animal_id`, `bin_start` and `bin_end` (POSIXct, UTC), `receiver_serial`,
#'   `station_id` when the listening table has it, `y` (integer), `n`
#'   (integer) and `listening` (logical). Attributes: `"bin"` (seconds),
#'   `"n_raised"`, `"swaps_merged"`, `"cover_dropped"`, `"not_listening"`
#'   (detections at receivers outside the listening table), `"codes"` (per
#'   animal: transmitters, and the rule when there are two or more) and
#'   `"rate_check"` (per animal: `bins`, median `n`, `q99_best` and `ratio`,
#'   which is `NA` under 100 bins).
#' @seealso [location_surface()], which reads the result.
#' @export
#' @examples
#' # The Stoney example: five animals, 70 receivers, 48 hours from 1 June 2024.
#' counts <- bin_detections(stoney_track$detections, stoney_track$receivers,
#'                          stoney_track$tags)
#' # One row per animal, hour and listening receiver; y is the count, n the
#' # transmissions sent in the hour.
#' head(counts[counts$y > 0, c("animal_id", "bin_start", "station_id", "y", "n")])
#' table(unique(counts[c("animal_id", "bin_start")])$animal_id)  # detected hours
#'
#' # The same detections by day, with the delay given as the shortest and
#' # longest delay in seconds.
#' daily <- bin_detections(stoney_track$detections, stoney_track$receivers,
#'                         tag_delay = c(120, 180), bin = "day")
#' daily[daily$y > 0, c("animal_id", "bin_start", "station_id", "y", "n")][1:5, ]
bin_detections <- function(detections, listening, tag_delay, bin = "hour", codes = NULL,
                           min_cover = 0.5) {
  bin_s <- bd_bin_seconds(bin)
  share_ok <- is.numeric(min_cover) && length(min_cover) == 1 && isTRUE(min_cover >= 0)
  if (!share_ok || min_cover > 1) stop("`min_cover` must be one number from 0 to 1.", call. = FALSE)
  keyed <- bd_station_keys(detections, listening)
  listening <- keyed$listening
  det <- bd_detections(keyed$detections, bin_s)
  delays <- bd_delays(tag_delay, det, codes)
  bins <- unique(det[, c("animal_id", "bin_start")])
  covered <- bd_cover(bd_listening_bins(bins, listening, bin_s), min_cover)
  swapped <- bd_merge_swaps(covered$heard, det)
  rows <- bd_join_counts(swapped$listening, swapped$detections)
  rows <- bd_transmissions(rows, delays, bin_s)
  out <- bd_finish(rows, bin_s, swapped$merged, covered$dropped)
  bd_rate_check(out, delays)
}

# ---- Deployments with no serial ----
# Listening rows that name a station but no receiver are matched to the
# detections at that station through one key, "station:<station_id>".
bd_station_keys <- function(detections, listening) {
  both <- is.list(listening) && !is.data.frame(listening)
  tabs <- if (both) listening else list(listening)
  bare <- unique(unlist(lapply(tabs, bd_bare_stations)))
  if (!length(bare)) return(list(detections = detections, listening = listening))
  if (!is.data.frame(detections) || !"station_id" %in% names(detections)) {
    bd_warn_bare(bare)
    return(list(detections = detections, listening = listening))
  }
  known <- unique(unlist(lapply(tabs, function(t) if (is.data.frame(t)) t$receiver_serial)))
  known <- known[!is.na(known)]
  at <- as.character(detections$station_id) %in% bare & !detections$receiver_serial %in% known
  detections$receiver_serial[at] <- paste0("station:", detections$station_id[at])
  message(wrap_text(sprintf(
    paste0("%s %s listening rows with no receiver_serial; %s there %s matched by ",
           "station (receiver_serial \"station:<station_id>\")."),
    plural(length(bare), "station"), agree(length(bare), "has", "have"),
    plural(sum(at), "detection"), agree(sum(at), "is", "are")
  )))
  list(detections = detections,
       listening = if (both) lapply(listening, bd_key_bare) else bd_key_bare(listening))
}

bd_key_bare <- function(t) {
  if (!length(bd_bare_stations(t))) return(t)
  na <- is.na(t$receiver_serial) & !is.na(t$station_id)
  t$receiver_serial[na] <- paste0("station:", t$station_id[na])
  t
}

bd_warn_bare <- function(bare) {
  warning(sprintf(paste0("%s in `listening` (%s) %s rows with no receiver_serial, ",
                         "and `detections` has no station_id column, so those rows cannot ",
                         "be matched to detections and count as no receiver. Give ",
                         "`detections` a station_id column, or the deployments their serials."),
                  plural(length(bare), "station"), bd_first(bare),
                  agree(length(bare), "has", "have")),
          call. = FALSE)
}

bd_bare_stations <- function(t) {
  if (!is.data.frame(t) || !all(c("receiver_serial", "station_id") %in% names(t))) return(NULL)
  unique(as.character(t$station_id[is.na(t$receiver_serial) & !is.na(t$station_id)]))
}

# ---- Bin length ----
bd_bin_seconds <- function(bin) {
  named <- c(hour = 3600, day = 86400)
  if (is.character(bin) && length(bin) == 1 && bin %in% names(named)) return(named[[bin]])
  if (bd_positive(bin)) return(as.numeric(bin))
  stop("`bin` must be \"hour\", \"day\" or one positive number of seconds.", call. = FALSE)
}

bd_positive <- function(x) {
  is.numeric(x) && length(x) == 1 && isTRUE(is.finite(x) && x > 0)
}

# ---- Column checks ----
# The contract asks a function that reads a shared table to stop at the
# first column that is missing or of the wrong type, and to name it.
bd_check_columns <- function(x, spec, what) {
  if (!is.data.frame(x)) stop(sprintf("`%s` must be a data frame.", what), call. = FALSE)
  for (col in names(spec)) {
    if (!col %in% names(x)) {
      stop(sprintf("`%s` has no column `%s`.", what, col), call. = FALSE)
    }
    if (!inherits(x[[col]], spec[[col]]) && !(spec[[col]] == "numeric" && is.numeric(x[[col]]))) {
      stop(wrap_text(sprintf("Column `%s` of `%s` must be %s, not %s.", col, what,
                             spec[[col]], class(x[[col]])[1])), call. = FALSE)
    }
  }
  invisible(x)
}

bd_check_utc <- function(x, col, what) {
  tz <- attr(x[[col]], "tzone")
  if (is.null(tz) || !identical(tz, "UTC")) {
    stop(wrap_text(sprintf("Column `%s` of `%s` must be in UTC, not \"%s\".", col, what,
                           if (is.null(tz)) "" else tz[1])), call. = FALSE)
  }
}

# ---- Detections ----
bd_detections <- function(detections, bin_s) {
  bd_check_columns(detections, c(datetime = "POSIXct", receiver_serial = "character",
                                 transmitter = "character", animal_id = "character"),
                   "detections")
  bd_check_utc(detections, "datetime", "detections")
  det <- data.table::data.table(
    animal_id = detections$animal_id, receiver_serial = detections$receiver_serial,
    transmitter = detections$transmitter, datetime = detections$datetime
  )
  no_animal <- is.na(det$animal_id)
  if (any(no_animal)) {
    message(wrap_text(sprintf("%s of %s detections have no animal_id and were not counted.",
                              fmt_n(sum(no_animal)),
                              fmt_n(nrow(det)))))
    det <- det[!no_animal]
  }
  if (!nrow(det)) stop("`detections` holds no detection with an animal_id.", call. = FALSE)
  if (anyNA(det$datetime) || anyNA(det$receiver_serial)) {
    stop("`datetime` and `receiver_serial` of `detections` must have no NA.", call. = FALSE)
  }
  start <- floor(as.numeric(det$datetime) / bin_s) * bin_s
  bin_start <- NULL
  det[, bin_start := as.POSIXct(start, origin = "1970-01-01", tz = "UTC")]
  det
}

# ---- Tag delays ----
# Returns one row per animal and transmitter code with its mean delay and the
# animal's rule for sharing transmissions between codes.
bd_delays <- function(tag_delay, det, codes) {
  if (!is.null(codes) && !identical(codes, "shared") && !identical(codes, "independent")) {
    stop("`codes` must be NULL, \"shared\" or \"independent\".", call. = FALSE)
  }
  seen <- unique(det[, c("animal_id", "transmitter")])
  tag_delay <- bd_delay_range(tag_delay)
  if (is.numeric(tag_delay) && length(tag_delay) == 1) {
    if (!bd_positive(tag_delay)) {
      stop("`tag_delay` must be a positive number of seconds.", call. = FALSE)
    }
    out <- seen[, c("delay", "rule") := list(as.numeric(tag_delay), NA_character_)]
    return(bd_code_rule(out, codes))
  }
  tab <- bd_delay_table(tag_delay)
  if ("animal_id" %in% names(tab)) {
    listed <- unique(tab[!is.na(animal_id) & animal_id %in% seen$animal_id,
                         c("animal_id", "transmitter")])
    seen <- unique(rbind(seen, listed))
  }
  out <- merge(seen, unique(tab[, c("transmitter", "delay", "rule")]), by = "transmitter",
               all.x = TRUE)
  missing <- unique(out$transmitter[is.na(out$delay)])
  if (length(missing)) {
    stop(sprintf("`tag_delay` has no delay for %s: %s.", plural(length(missing), "transmitter"),
                 bd_first(missing)), call. = FALSE)
  }
  bd_code_rule(out, codes)
}

# Two numbers are the shortest and longest delay, as a simulator or a tag
# sheet gives them; their mean is the delay.
bd_delay_range <- function(tag_delay) {
  if (!is.numeric(tag_delay) || length(tag_delay) != 2) return(tag_delay)
  ok <- all(is.finite(tag_delay)) && all(tag_delay > 0) && tag_delay[1] <= tag_delay[2]
  if (!ok) {
    stop("Two numbers in `tag_delay` must be the shortest and longest delay, both positive.",
         call. = FALSE)
  }
  mean(tag_delay)
}

# One row per transmitter: its delay and, from an optional `codes` column,
# its rule.
bd_delay_table <- function(tag_delay) {
  if (!is.data.frame(tag_delay) || !"transmitter" %in% names(tag_delay)) {
    stop(wrap_text(paste0("`tag_delay` must be one number, two numbers (shortest and ",
                          "longest delay) or a data frame with a `transmitter` column.")),
         call. = FALSE)
  }
  tab <- data.table::as.data.table(tag_delay)
  delay <- k <- rule <- NULL
  if ("delay" %in% names(tab)) {
    tab[, delay := suppressWarnings(as.numeric(delay))]
  } else if (all(c("min_delay", "max_delay") %in% names(tab))) {
    lo <- suppressWarnings(as.numeric(tab$min_delay))
    hi <- suppressWarnings(as.numeric(tab$max_delay))
    tab[, delay := (lo + hi) / 2]
  } else {
    stop("`tag_delay` needs a `delay` column, or `min_delay` and `max_delay`.", call. = FALSE)
  }
  bad <- !is.na(tab$delay) & tab$delay <= 0
  if (any(bad)) stop("Every delay in `tag_delay` must be positive.", call. = FALSE)
  tab[, transmitter := as.character(transmitter)]
  tab[, rule := if ("codes" %in% names(tab)) as.character(tab$codes) else NA_character_]
  odd <- setdiff(stats::na.omit(tab$rule), c("shared", "independent"))
  if (length(odd)) {
    stop(sprintf("The `codes` column of `tag_delay` holds \"%s\"; use \"shared\" or ",
                 odd[1]), "\"independent\".", call. = FALSE)
  }
  tab <- tab[!is.na(delay)]
  dup <- tab[, list(k = max(data.table::uniqueN(delay), data.table::uniqueN(rule, na.rm = TRUE))),
             by = "transmitter"][k > 1]
  if (nrow(dup)) {
    stop(sprintf("`tag_delay` gives two delays or two rules for %s.", bd_first(dup$transmitter)),
         call. = FALSE)
  }
  tab[, rule := rule[!is.na(rule)][1], by = "transmitter"]
  tab
}

# Each animal gets one rule: from the `codes` column where given, else from
# the argument. An animal with two or more codes and no rule stops the run.
bd_code_rule <- function(delays, codes) {
  rule <- k <- rules <- NULL
  if (!is.null(codes)) delays[is.na(rule), rule := codes]
  per <- delays[, list(k = .N, rules = data.table::uniqueN(rule, na.rm = TRUE)),
                by = "animal_id"]
  if (any(per$rules > 1)) {
    stop(sprintf("Animal %s: its transmitters are given both rules in `tag_delay$codes`.",
                 per$animal_id[per$rules > 1][1]), call. = FALSE)
  }
  open <- per[k > 1 & rules == 0]$animal_id
  if (length(open)) {
    many <- length(open) > 1
    msg <- paste0("%s carr%s more than one transmitter (%s). Say how a tag's codes ",
                  "share its transmissions: codes = \"shared\" when the codes take turns on ",
                  "one schedule, or codes = \"independent\" when each keeps its own; or give ",
                  "`tag_delay` a `codes` column. See ?bin_detections.")
    stop(sprintf(msg, plural(length(open), "animal"), if (many) "y" else "ies",
                 bd_first(open)), call. = FALSE)
  }
  delays[, rule := rule[!is.na(rule)][1], by = "animal_id"]
  delays
}

bd_first <- function(x, k = 5) {
  x <- unique(x)
  paste0(paste(utils::head(x, k), collapse = ", "), if (length(x) > k) ", ..." else "")
}

# ---- Listening receivers per bin ----
# Returns one row per animal, bin and listening receiver, with the station, a
# rank that orders receivers at one station by when they went in, and the
# share of the bin the receiver was in the water.
bd_listening_bins <- function(bins, listening, bin_s) {
  both <- is.list(listening) && !is.data.frame(listening)
  if (both && !all(c("deployments", "days") %in% names(listening))) {
    stop("A list given as `listening` must hold `deployments` and `days`.", call. = FALSE)
  }
  bd_check_rows(if (both) listening[c("deployments", "days")] else list(listening = listening),
                both)
  heard <- if (both) {
    bd_from_both(bins, listening, bin_s)
  } else if (is.data.frame(listening) && "deploy_datetime" %in% names(listening)) {
    bd_from_deployments(bins, listening, bin_s)
  } else {
    bd_from_days(bins, listening, bin_s)
  }
  # A receiver moved between stations within one bin keeps its later
  # deployment.
  data.table::setorderv(heard, c("animal_id", "bin_start", "receiver_serial", "rank"))
  last <- !duplicated(heard, by = c("animal_id", "bin_start", "receiver_serial"), fromLast = TRUE)
  heard[last, c("animal_id", "bin_start", "receiver_serial", "station_id", "rank", "cover")]
}

# An empty listening table would make every receiver-bin silent or absent.
bd_check_rows <- function(tabs, both) {
  for (nm in names(tabs)) {
    if (is.data.frame(tabs[[nm]]) && !nrow(tabs[[nm]])) {
      stop(sprintf("`%s` has no rows, so no receiver listened in any bin.",
                   if (both) paste0("listening$", nm) else "listening"), call. = FALSE)
    }
  }
  invisible(NULL)
}

bd_from_deployments <- function(bins, listening, bin_s) {
  bd_check_columns(listening, c(receiver_serial = "character", station_id = "character",
                                deploy_datetime = "POSIXct", recover_datetime = "POSIXct"),
                   "listening")
  bd_check_utc(listening, "deploy_datetime", "listening")
  bd_check_utc(listening, "recover_datetime", "listening")
  end <- as.numeric(listening$recover_datetime)
  start <- bs <- be <- d_start <- d_end <- cover <- NULL
  dep <- data.table::data.table(
    receiver_serial = listening$receiver_serial, station_id = listening$station_id,
    start = as.numeric(listening$deploy_datetime), end = data.table::fifelse(is.na(end), Inf, end)
  )
  dep[, c("rank", "d_start", "d_end") := list(start, start, end)]
  b <- data.table::data.table(animal_id = bins$animal_id, bin_start = bins$bin_start,
                              bs = as.numeric(bins$bin_start))
  b[, be := bs + bin_s]
  # In a non-equi join the join columns of the result hold the values of b,
  # so the deployment's own times travel as rank, d_start and d_end.
  out <- dep[b, on = list(start < be, end > bs), nomatch = NULL, allow.cartesian = TRUE]
  bs0 <- as.numeric(out$bin_start)
  out[, cover := (pmin(d_end, bs0 + bin_s) - pmax(d_start, bs0)) / bin_s]
  out[, c("animal_id", "bin_start", "receiver_serial", "station_id", "rank", "cover")]
}

bd_from_days <- function(bins, listening, bin_s) {
  days <- bd_days(listening, bin_s)
  bin_start <- date <- NULL
  b <- data.table::data.table(animal_id = bins$animal_id, bin_start = bins$bin_start)
  b[, date := as.Date(bin_start, tz = "UTC")]
  out <- merge(b, days, by = "date", allow.cartesian = TRUE)
  out[, c("date", "cover") := list(NULL, 1)]
  out
}

# Deployments give the hours in the water; the listed days remove days that
# failed a check.
bd_from_both <- function(bins, listening, bin_s) {
  days <- bd_days(listening$days, bin_s)
  heard <- bd_from_deployments(bins, listening$deployments, bin_s)
  date <- bin_start <- NULL
  heard[, date := as.Date(bin_start, tz = "UTC")]
  ok <- unique(days[, c("receiver_serial", "date")])
  out <- merge(heard, ok, by = c("receiver_serial", "date"))
  out[, date := NULL]
  out
}

# Listed receiver-days, with a rank: a receiver's first listed day at a
# station stands in for its deployment time, so the receiver that arrived
# later at a station ranks higher.
bd_days <- function(listening, bin_s) {
  bd_check_columns(listening, c(receiver_serial = "character", date = "Date"), "listening")
  if (86400 %% bin_s != 0) {
    stop("With a receiver-day listening table the bin must divide a day; give deployments ",
         "alone for other bins.", call. = FALSE)
  }
  date <- station_id <- rank <- NULL
  days <- data.table::data.table(receiver_serial = listening$receiver_serial,
                                 date = listening$date)
  days[, station_id := if ("station_id" %in% names(listening)) {
    as.character(listening$station_id)
  } else {
    NA_character_
  }]
  if ("listening" %in% names(listening)) days <- days[listening$listening %in% TRUE]
  days <- unique(days[!is.na(date) & !is.na(receiver_serial)])
  if (nrow(days)) days[, rank := as.numeric(min(date)), by = c("receiver_serial", "station_id")]
  days
}

# A receiver listens in a bin it covers for at least min_cover; at a station
# the covers of its receivers add up, so a swap within the bin keeps it.
bd_cover <- function(heard, min_cover) {
  share <- cover <- NULL
  heard[, share := cover]
  at_station <- !is.na(heard$station_id)
  heard[at_station, share := pmin(1, sum(cover)), by = c("animal_id", "bin_start", "station_id")]
  keep <- heard$share > 0 & heard$share >= min_cover
  out <- heard[keep]
  out[, c("cover", "share") := NULL]
  list(heard = out, dropped = sum(!keep))
}

# ---- Receiver swaps ----
bd_merge_swaps <- function(heard, det) {
  y <- station_id <- rank <- keep <- NULL
  counts <- det[, list(y = .N), by = c("animal_id", "bin_start", "receiver_serial")]
  if (all(is.na(heard$station_id))) {
    return(list(listening = heard, detections = counts, merged = 0L))
  }
  by_station <- c("animal_id", "bin_start", "station_id")
  heard[, keep := receiver_serial[order(rank, receiver_serial)][.N], by = by_station]
  heard[is.na(station_id), keep := receiver_serial]
  moved <- heard[receiver_serial != keep, c("animal_id", "bin_start", "receiver_serial", "keep")]
  if (nrow(moved)) {
    counts <- merge(counts, moved, by = c("animal_id", "bin_start", "receiver_serial"),
                    all.x = TRUE)
    counts[!is.na(keep), receiver_serial := keep]
    counts <- counts[, list(y = sum(y)), by = c("animal_id", "bin_start", "receiver_serial")]
  }
  heard <- heard[receiver_serial == keep]
  heard[, keep := NULL]
  list(listening = heard, detections = counts, merged = nrow(moved))
}

# ---- Join listening receivers and counts ----
bd_join_counts <- function(heard, counts) {
  y <- rank <- listening <- NULL
  rows <- merge(heard, counts, by = c("animal_id", "bin_start", "receiver_serial"), all = TRUE)
  rows[, listening := !is.na(rank)]
  rows[is.na(y), y := 0L]
  rows[, rank := NULL]
  rows
}

# ---- Transmissions per bin and the floor ----
bd_transmissions <- function(rows, delays, bin_s) {
  n <- y <- delay <- n_nominal <- y_max <- rule <- NULL
  # Codes that take turns share the tag's transmissions; independent codes
  # add theirs.
  per_animal <- delays[, list(n_nominal = round(if (identical(rule[1], "shared")) {
    bin_s / mean(delay)
  } else {
    sum(bin_s / delay)
  })), by = "animal_id"]
  rows <- merge(rows, per_animal, by = "animal_id", all.x = TRUE)
  rows[, y_max := max(y), by = c("animal_id", "bin_start")]
  rows[, n := as.integer(pmax(n_nominal, y_max))]
  raised <- unique(rows[y_max > n_nominal, c("animal_id", "bin_start")])
  attr(rows, "n_raised") <- nrow(raised)
  rows[, c("n_nominal", "y_max") := NULL]
  rows
}

# ---- Assemble the result ----
bd_finish <- function(rows, bin_s, merged, dropped) {
  listening <- NULL
  n_raised <- attr(rows, "n_raised")
  n_bins <- nrow(unique(rows[, c("animal_id", "bin_start")]))
  data.table::setorderv(rows, c("animal_id", "bin_start", "receiver_serial"))
  out <- data.frame(
    animal_id = rows$animal_id, bin_start = rows$bin_start,
    bin_end = rows$bin_start + bin_s, receiver_serial = rows$receiver_serial,
    stringsAsFactors = FALSE
  )
  if (!all(is.na(rows$station_id))) out$station_id <- rows$station_id
  out$y <- as.integer(rows$y)
  out$n <- rows$n
  out$listening <- rows$listening
  off <- rows[listening == FALSE]
  if (nrow(off)) {
    warning(sprintf(paste0("%s in %s %s logged by receivers not listening ",
                           "in the bin (not listed, or in the water for less than `min_cover` ",
                           "of it); their rows have listening = FALSE and location_surface() ",
                           "leaves them out."), plural(sum(off$y), "detection"),
                    plural(nrow(unique(off[, 1:2])), "bin"), agree(sum(off$y), "was", "were")),
            call. = FALSE)
  }
  if (n_raised > 0) {
    message(wrap_text(sprintf("n was raised to the largest receiver count in %s of %s bins.",
                              fmt_n(n_raised), fmt_n(n_bins))))
  }
  attr(out, "bin") <- bin_s
  attr(out, "n_raised") <- n_raised
  attr(out, "swaps_merged") <- merged
  attr(out, "cover_dropped") <- dropped
  attr(out, "not_listening") <- sum(off$y)
  out
}

# ---- Rate check ----
# With the right delay and rule, the busiest receiver of a well-detected
# animal logs close to every transmission in its best bins, so the 99th
# percentile of its count sits near n. The ratio is then near 1 for a
# well-detected animal, and near one half for a two-code animal whose codes
# share the tag's transmissions but whose rates are summed. The cut-offs sit
# between.
bd_rate_check <- function(out, delays) {
  rule <- y <- n <- bins <- q99_best <- ratio <- NULL
  info <- delays[, list(transmitters = .N, rule = if (.N > 1) rule[1] else NA_character_),
                 by = "animal_id"]
  d <- data.table::data.table(animal_id = out$animal_id, bin_start = out$bin_start, y = out$y,
                              n = out$n)[out$listening]
  # .N guards max() when no receiver listened: data.table evaluates j once
  # on empty columns to learn their types.
  best <- d[, list(y = if (.N) max(y) else NA_integer_, n = n[1]),
            by = c("animal_id", "bin_start")]
  # as.numeric: the median of an odd number of integers is an integer, and
  # data.table needs one type across groups.
  chk <- best[, list(bins = .N, n = as.numeric(stats::median(n)),
                     q99_best = as.numeric(stats::quantile(y, 0.99, names = FALSE))),
              by = "animal_id"]
  chk[, ratio := ifelse(bins >= 100, q99_best / n, NA_real_)]
  off <- chk$animal_id[!is.na(chk$ratio) & (chk$ratio < 0.6 | chk$ratio > 1.5)]
  if (length(off)) {
    warning(sprintf(paste0("For %s the busiest receiver's 99th-percentile count per ",
                           "bin is below 0.6 or above 1.5 times n, so n may be wrong: %s. ",
                           "Check `tag_delay` and `codes`; see attr(, \"rate_check\")."),
                    plural(length(off), "animal"), bd_first(off)), call. = FALSE)
  }
  attr(out, "codes") <- as.data.frame(info)
  attr(out, "rate_check") <- as.data.frame(chk)
  out
}
