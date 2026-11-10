#' Location surfaces from detection counts
#'
#' @description
#' For each animal and time bin, computes the posterior probability that the
#' animal was in each water cell of the field's grid, given which receivers
#' were listening and how many detections each logged. The surface is the
#' likelihood of the bin's counts in each water cell under a flat prior over
#' the water cells, normalised to sum to one over them. It is therefore a
#' posterior by Bayes' rule. The default likelihood is binomial and uses the
#' detections and the misses (binomial grid localization). With
#' `likelihood = "presence_only"` it is the product of the detection
#' probabilities of the detections alone (presence-only grid localization).
#' Each bin is computed on its own; no movement model links bins.
#'
#' @details
#' **The model.** In bin \eqn{t}, an animal in cell \eqn{c} sends \eqn{n_t}
#' transmissions, the column `n` of `counts`, taken as known. Each listening
#' receiver's count is binomial in \eqn{n_t} with the field's detection
#' probability at that cell, the receivers independent given the cell. The
#' surface is the product of these terms over the listening receivers,
#' normalised over every water cell of the grid, including cells no receiver
#' can hear: a flat prior over water. Silent receivers count: a receiver that
#' heard nothing makes the cells near it less likely. The article Hourly
#' positions against the true track writes the model out and explains how
#' \eqn{n_t} is set, in its section The model:
#' <https://jakebrownscombe.github.io/detfield/articles/positioning.html#the-model>.
#'
#' The model makes six assumptions. The animal stays in one cell for the bin.
#' The tag sits at the cell centre, and its depth is ignored. Each receiver's
#' detection probability is constant within the bin. The field is known.
#' Transmissions are detected independently given the cell. Every detection is
#' a true transmission of this tag, with no false detections or code
#' collisions. A receiver whose deployment covers at least
#' `min_cover` of the bin (0.5 by default in [bin_detections()]) is taken to
#' listen for all of it, so it counts every transmission of the bin.
#' With `blocked = 0` in [de_model()], the clear-line curve is read across
#' land; [de_field()] gives the count of such pairs in a message.
#'
#' **Two readings of the curve.** `likelihood = "binomial"`, the default,
#' is binomial grid localization: each listening receiver's count is
#' binomial in the transmissions sent, so its misses and the silence of other
#' receivers are evidence too. It is the binomial count model of Winton et
#' al. (2018) with the detection function fixed in advance rather than
#' estimated, on a grid of water cells rather than in continuous space.
#' Winton et al. take the number of trials as the maximum expected number of
#' transmissions in the interval, and note that a Poisson count may suit
#' tags whose transmission intervals are programmed to vary, as random
#' delays are. In the terms of Hostetter and Royle (2020) it is independent
#' localization of one bin, with all the bin's transmissions sent from one
#' position. They describe fixing the detection function from a range test as
#' the common practice, and recommend estimating the detection parameters
#' from the telemetry data instead.
#'
#' Use the binomial reading by default. It places an animal that one
#' receiver heard rarely on the ring of cells where that receiver's
#' probability matches its count. Its credible regions are calibrated when
#' the curve is right for tagged animals at the bin scale and the animal
#' stays in one cell for the bin. In the example study the curve is the one
#' that drew the detections, but the animals swim within the hour: the 95%
#' regions held 35.3% of the true minutes and the hour's mean true position
#' in 156 of 212 hours (74%). The article Hourly positions against the true
#' track computes both, in its section Credible regions and their coverage:
#' <https://jakebrownscombe.github.io/detfield/articles/positioning.html>.
#'
#' `likelihood = "presence_only"` is presence-only grid localization: the
#' log-likelihood of cell \eqn{c} is \eqn{\sum_s y_{st} \log p_s(c)}, the
#' product of \eqn{p_s(c)} over the detections, with misses and silent
#' receivers dropped. It never places an animal away from the receivers that
#' heard it; use it as a check when the curve's level is not trusted. Its
#' credible regions are not calibrated: the surface narrows as the count
#' grows, whatever the curve's level. On the example data (`stoney_track`,
#' 212 detected hours of five animals, 12,720 minutes), its 95% regions held
#' the animal's true cell in 23.8% of the minutes, against 35.3% for the
#' binomial reading; neither reaches 95% there, because the animals move within the
#' hour. Compare its region sizes between bins only. Under it `overdispersion` does not apply, and
#' `"station_hour"` is ignored with a warning.
#' `misfit_p` and the near diagnostics are computed against the binomial
#' prediction under either reading.
#'
#' **Overdispersion.** Detection probability varies between receivers and
#' hours more than a range curve says. With `overdispersion =
#' "station_hour"`, each receiver-bin carries its own logit-scale effect,
#' normal with standard deviation `field$model$sigma$station_hour`, which is
#' integrated out by Gauss-Hermite quadrature. This standard deviation is the
#' whole spread of one receiver's logit about the curve within one bin, with
#' all sources of variation combined.
#'
#' `"station_hour"` is the default. It applies the spread when the field's
#' model gives a standard deviation above 0, and runs as the plain binomial
#' when it gives 0. When the model gives none, as [de_model()] and
#' [de_logistic()] do by default, the call stops. This is why every example
#' sets `sigma = list(station_hour = 0)` in [de_model()].
#'
#' The standard deviation is defined only together with the curve it was
#' estimated around. Paired with a curve that is too flat or too low at the
#' bin scale, a large standard deviation lets a strong count be explained
#' from far away, and the surface moves away from the receivers that heard
#' the animal. Estimate the curve and the standard deviation together, from
#' data at the bin scale (for example range-test windows of the bin length),
#' and pass both; a value estimated around one curve does not transfer to
#' another. A value that is too small for its curve makes regions too small.
#' The form also assumes the effect is normal on the logit scale. When the
#' truth differs, for example when a receiver either hears the tag in an hour
#' or does not, regions can hold the truth less often than their level.
#'
#' The standard deviation is a property of the detection model and is never
#' estimated here. It describes hour bins: under `"station_hour"`, the
#' default, with a standard deviation above 0, a warning says so when a bin
#' is not an hour. With `"none"` the receivers' counts are independent
#' binomials, which is overconfident whenever detection varies between hours.
#'
#' **Computation.** All sums are in log space, so a bin with thousands of
#' transmissions and no detection stays finite. Under `"station_hour"` each
#' receiver's term is read from a table over the logit, built once per count
#' and number of transmissions by adaptive Gauss-Hermite quadrature, which
#' transforms the variable of integration so the nodes sample the integrand
#' where its mass lies (Liu and Pierce 1994), and
#' interpolated by cubic Hermite polynomials with exact derivatives. The
#' number of nodes grows with the standard deviation (25 per logit, from 20
#' to 200). The script
#' <https://github.com/jakebrownscombe/detfield/blob/main/validation/02_quadrature_error.R>
#' compares each term with numerical integration over counts 0 to n, for n of
#' 30, 60 and 720. It is in the repository and not in the installed package.
#' The largest error in log probability is 1.2e-6 at a standard deviation of
#' 2.55, 6.0e-6 at 3.4 and 1.1e-4 at 6. A smaller `control$nodes` gives a
#' warning, because the error grows fast as nodes fall: 20 nodes give 3.2e-4,
#' 1.8e-3 and 1.8e-2 at the same standard deviations. On the example data
#' (`stoney_track`, 212 detected animal-hours, 4,818 water cells) at a
#' standard deviation of 2.55, surfaces were compared with surfaces from 200
#' nodes. Surfaces from 1 node were off by up to 0.0339 in total variation,
#' from 20 nodes by up to 1.15e-4, and from the default of 64 nodes by up to
#' 2.43e-7. The script
#' <https://github.com/jakebrownscombe/detfield/blob/main/validation/01_node_error.R>
#' computes these; it is in the repository and not in the installed package.
#' Above a standard deviation of 8 the 200-node cap applies, also with a
#' warning. Below a logit at which the detection probability is
#' negligible, the term is computed in closed form, so a field with very
#' negative logits costs no more than one without.
#'
#' **Cost.** Each call has a fixed cost per animal, for the tables and the
#' silent receivers' terms, before any bin. On the example data at a standard
#' deviation of 2.55, one hour of `stoney-01` took 0.17 s and all 43 of its
#' hours took 2.59 s. These are medians of five runs. The script
#' <https://github.com/jakebrownscombe/detfield/blob/main/validation/03_location_surface_timing.R>
#' prints them; it is in the repository and not in the installed package. Give
#' all of an animal's bins in one call: a call per day pays the fixed cost
#' once a day.
#'
#' **Large runs.** The kept cells of every bin come back as one table, and a
#' table holds at most 2,147,483,647 rows. A run that would exceed it stops
#' with a message. With `control$write_dir`, each animal's surfaces are
#' written to their own file as they are computed, so the session never
#' holds them all; the call then returns a table of the files.
#'
#' **What is kept.** Each surface keeps only the smallest set of cells that
#' holds `control$kept_mass` (default 0.999) of its probability, because
#' every water cell of every bin would not fit in memory for a large array.
#' Probabilities are not rescaled, so they sum to the kept mass.
#'
#' **Misfit.** For each bin, `misfit_p` in the bin table is the smallest,
#' over the receivers that detected the animal, of the posterior predictive
#' probability of logging at least the count that receiver logged, and
#' `misfit_receiver` names that receiver. A value near zero means no cell of
#' the surface explains that receiver's count under the field: the field may
#' be wrong near that receiver, or the animal moved within the bin.
#'
#' **Distance from the receivers that heard the animal.** A surface can fit
#' every count and still sit away from the receivers that logged them. This
#' happens with a curve that is too high near a receiver. A low count is then
#' best explained by a ring of cells at the distance where the curve gives
#' that rate, and `misfit_p` stays high. Three columns of the bin table show this. `mode_dist_m`
#' is the straight-line distance from the centre of the most probable cell
#' to the nearest receiver that detected the animal in the bin. `near_mass`
#' is the probability within `near_m` metres (straight line) of
#' `near_receiver`, the receiver that logged the most detections in the bin
#' (the first in row order on a tie). `near_m` is `control$near_m` when given;
#' otherwise it is the longest clear-line distance (`in_water`) from that
#' receiver to a water cell centre at which the field's detection probability
#' is still 0.5 or more, read from the field's distance table. Because it is
#' read off cell centres, it falls short of the curve's d50 by up to one cell
#' diagonal. On the example data it is 224 m for receiver S069, whose curve
#' gives d50 253 m at its depth of 1.2 m. Where the field holds
#' no distance table, or the receiver has no such pair, `near_m` and
#' `near_mass` are `NA`; set `control$near_m`. `near_mass` is computed over
#' every water cell, before cells are dropped by `kept_mass`. Distances are in
#' the units of the grid's CRS, metres for any grid [water_distances()]
#' accepts. [print()] counts the bins with `near_mass` below 0.001.
#'
#' **Checks.** A surface that is not a distribution is degenerate: a
#' log-likelihood that is NaN, zero likelihood in every cell, or the same
#' likelihood in every cell. By default it stops the run with the animal and
#' bin named; with `control$on_degenerate = "skip"`, the bin is left out,
#' recorded in the attribute `skipped`, and a warning gives the count.
#'
#' **Matching receivers.** Rows of `counts` match columns of `field` by
#' `receiver_serial`; when the field's columns are stations, by
#' `station_id`. A listening receiver with no column in the field stops the
#' run, because leaving it out would drop its evidence without notice. Every
#' cell of the field must be a water cell (not `NA`) of its grid.
#'
#' **Detections at receivers not listening.** Rows with `listening = FALSE`
#' are left out. A bin whose detections all sit on such rows is left out as
#' well, because without them it would be computed as a bin in which every
#' listening receiver was silent, which puts the animal away from every
#' receiver. Such bins are counted in a warning and listed in the attribute
#' `not_listening_bins`.
#'
#' @param counts Counts from [bin_detections()], or a data frame with the same
#'   columns: `animal_id`, `bin_start`, `bin_end`, `receiver_serial`, `y`, `n`
#'   and `listening`. Rows with `listening = FALSE` are left out. When the
#'   field is keyed by `station_id` and `counts` has no `station_id`, each
#'   `receiver_serial` is mapped to its station on the bin's date through the
#'   field's deployment lookup, when the field holds one.
#' @param field A detection field (class `"de_field"`): a list holding `eta`,
#'   the logit of the per-transmission detection probability with one row per
#'   water cell and one column per receiver; `cells`, the grid cell numbers of
#'   the rows; `receivers`, the column names; `grid`, the raster or the path
#'   of a raster file; and `model`, whose `sigma$station_hour` sets the
#'   overdispersion.
#' @param overdispersion `"station_hour"` (the default) or `"none"`. Not used
#'   when `likelihood = "presence_only"`.
#' @param likelihood `"binomial"` (the default) or `"presence_only"`; see Two
#'   readings of the curve.
#' @param control A list of settings: `cores`, the number of processes for
#'   animals (default 1; forked processes, so 1 on Windows); `nodes`, the
#'   Gauss-Hermite nodes (default `NULL`, set from the standard deviation as
#'   above); `kept_mass`, the share of each surface kept (default 0.999);
#'   `table_step`, the logit step of the term tables (default 0.01);
#'   `on_degenerate`, `"stop"` (default) or `"skip"`; `write_dir`, a
#'   directory to which each animal's surfaces are written (default `NULL`,
#'   return them); `near_m`, the radius in metres of `near_mass` (default
#'   `NULL`, per receiver the farthest cell centre at p >= 0.5; see Details).
#' @return A data frame of class `"loc_surface"` with the columns
#'   `animal_id`, `bin_start`, `cell` and `prob`, one row per kept cell. Its
#'   attributes hold `bins` (one row per animal and bin: `bin_end`,
#'   `n_transmissions`, `n_listening`, `n_detecting`, `y_total`, `log_lik`,
#'   the log of the mean likelihood over water cells, `kept_mass`, `n_cells`,
#'   `misfit_p`, `misfit_receiver`, `mode_cell`, `mode_dist_m`,
#'   `near_receiver`, `near_m` and `near_mass`), `likelihood` (the reading
#'   used), `overdispersion` (form,
#'   sigma, nodes),
#'   `kept` (the rule), `grid` (packed), `n_water_cells`, `receivers` (the
#'   listening rows used), `receiver_xy` (receiver coordinates, when the field
#'   holds them), `skipped` (degenerate bins left out: animal, bin and
#'   reason) and `not_listening_bins` (bins left out because all their
#'   detections were at receivers not listening: animal, bin and detections).
#'
#'   With `control$write_dir`, each animal with at least one surface is saved
#'   by [saveRDS()] as a `"loc_surface"` of its own, and the call returns a
#'   data frame with one row per animal: `animal_id`, `file`, `bins` and
#'   `rows`, with the attributes `bins`, `overdispersion`, `skipped` and
#'   `not_listening_bins` as above.
#' @references
#' Hostetter, N. J. and Royle, J. A. (2020). Movement-assisted localization
#' from acoustic telemetry data. *Movement Ecology*, 8, 15.
#' \doi{10.1186/s40462-020-00199-6}
#'
#' Liu, Q. and Pierce, D. A. (1994). A note on Gauss-Hermite quadrature.
#' *Biometrika*, 81, 624-629. \doi{10.1093/biomet/81.3.624}
#'
#' Winton, M. V., Kneebone, J., Zemeckis, D. R. and Fay, G. (2018). A spatial
#' point process model to estimate individual centres of activity from
#' passive acoustic telemetry data. *Methods in Ecology and Evolution*, 9,
#' 2262-2272. \doi{10.1111/2041-210X.13080}
#' @seealso [bin_detections()] for the counts; [credible_region()] and
#'   [sample_locations()] for what to do with the surfaces.
#' @export
#' @examples
#' # Stoney Lake: a detection field from distances through water, the
#' # detections counted per hour, and one surface per detected hour.
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
#' surface
#' bins <- attr(surface, "bins")
#' head(bins[, c("animal_id", "bin_start", "n_cells", "misfit_p", "mode_dist_m")])
location_surface <- function(counts, field, overdispersion = c("station_hour", "none"),
                             likelihood = c("binomial", "presence_only"), control = list()) {
  likelihood <- match.arg(likelihood)
  asked <- !missing(overdispersion)
  overdispersion <- match.arg(overdispersion)
  overdispersion <- ls_reading(likelihood, overdispersion, asked)
  ctrl <- ls_control(control)
  rows <- ls_counts(counts)
  lost <- attr(rows, "not_listening_bins")
  rows <- ls_station_from_serial(rows, field)
  fld <- ls_field(field, rows)
  if (!is.null(fld$near)) fld$near$near_m <- ctrl$near_m
  ls_check_duplicates(rows, fld$key)
  col <- NULL
  rows[, col := match(rows[[fld$key]], fld$receivers)]
  od <- ls_overdispersion(overdispersion, fld$sigma, ctrl, rows)
  od$likelihood <- likelihood
  ls_warn_size(rows, fld, ctrl)
  pieces <- split(rows, by = "animal_id", keep.by = TRUE)
  res <- parallel_apply(pieces, function(a) ls_animal_out(a, fld, od, ctrl), ctrl$cores)
  if (!is.null(ctrl$write_dir)) return(ls_files(res, od, lost))
  ls_assemble(res, rows, fld, od, ctrl, lost)
}

# ---- Settings ----
ls_control <- function(control) {
  defaults <- list(cores = 1L, nodes = NULL, kept_mass = 0.999, table_step = 0.01,
                   on_degenerate = "stop", write_dir = NULL, near_m = NULL)
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop(sprintf("Unknown setting in `control`: %s. Known: %s.", paste(unknown, collapse = ", "),
                 paste(names(defaults), collapse = ", ")), call. = FALSE)
  }
  ctrl <- utils::modifyList(defaults, control, keep.null = TRUE)
  ok <- vapply(names(ls_settings), function(nm) isTRUE(ls_settings[[nm]](ctrl[[nm]])), logical(1))
  if (!all(ok)) {
    stop(sprintf("`control$%s` is out of range; see ?location_surface.", names(ok)[!ok][1]),
         call. = FALSE)
  }
  ctrl$cores <- as.integer(ctrl$cores)
  if (!is.null(ctrl$nodes)) ctrl$nodes <- as.integer(ctrl$nodes)
  if (!is.null(ctrl$write_dir)) {
    dir.create(ctrl$write_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(ctrl$write_dir)) {
      stop(sprintf("`control$write_dir` %s could not be created.", ctrl$write_dir), call. = FALSE)
    }
    ctrl$write_dir <- normalizePath(ctrl$write_dir)
  }
  ctrl
}

# The allowed range of each setting.
ls_settings <- list(
  cores = function(v) ls_whole(v) && v >= 1,
  nodes = function(v) is.null(v) || (ls_whole(v) && v >= 1 && v <= 200),
  kept_mass = function(v) bd_positive(v) && v <= 1,
  table_step = function(v) bd_positive(v) && v <= 0.1,
  on_degenerate = function(v) identical(v, "stop") || identical(v, "skip"),
  write_dir = function(v) is.null(v) || is_one_string(v),
  near_m = function(v) is.null(v) || (bd_positive(v) && is.finite(v))
)

ls_whole <- function(x) {
  is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x)
}

# ---- Counts ----
ls_counts <- function(counts) {
  bd_check_columns(counts, c(animal_id = "character", bin_start = "POSIXct",
                             bin_end = "POSIXct", receiver_serial = "character",
                             y = "numeric", n = "numeric", listening = "logical"), "counts")
  bd_check_utc(counts, "bin_start", "counts")
  keep <- intersect(c("animal_id", "bin_start", "bin_end", "receiver_serial", "station_id",
                      "y", "n", "listening"), names(counts))
  # Columns are taken by name with [[, because `[` on a keyed data.table
  # with a character index joins on the key instead of selecting columns.
  rows <- data.table::as.data.table(lapply(stats::setNames(keep, keep), function(k) counts[[k]]))
  whole <- function(v) !anyNA(v) && all(v >= 0 & v == round(v))
  if (!whole(rows$y) || !whole(rows$n)) {
    stop("`y` and `n` of `counts` must be whole numbers of at least 0, not NA.", call. = FALSE)
  }
  off <- !(rows$listening %in% TRUE)
  lost <- NULL
  if (any(off)) {
    bin_of <- paste(rows$animal_id, as.numeric(rows$bin_start))
    gone <- length(setdiff(bin_of[off], bin_of[!off]))
    message(wrap_text(sprintf("%s (%s) with listening = FALSE %s left out%s",
                              plural(sum(off), "row"), plural(sum(rows$y[off]), "detection"),
                              agree(sum(off), "was", "were"),
                              if (gone > 0) {
                                sprintf(", and with them %s with no listening receiver.",
                                        plural(gone, "bin"))
                              } else {
                                "."
                              })))
    lost <- ls_unheard_bins(rows, off)
    rows <- rows[!off & !bin_of %in% lost$key]
  }
  if (!nrow(rows)) {
    stop(sprintf(paste0("`counts` has no listening receiver: all %s %s listening = ",
                        "FALSE. Check that the listening table names the receivers that ",
                        "logged the detections (bin_detections() warns when it does not)."),
                 plural(length(off), "row"), agree(length(off), "has", "have")), call. = FALSE)
  }
  rows <- ls_floor(rows)
  data.table::setattr(rows, "not_listening_bins", lost$table)
  rows
}

# Bins whose detections all sit on rows that are not listening. Computed
# from their listening rows alone, they would be bins in which every
# receiver was silent, so they are left out and listed.
ls_unheard_bins <- function(rows, off) {
  y_off <- y_on <- live <- NULL
  d <- data.table::data.table(animal_id = rows$animal_id, bin_start = rows$bin_start,
                              y_off = ifelse(off, rows$y, 0), y_on = ifelse(off, 0, rows$y),
                              live = !off)
  per <- d[, list(detections = sum(y_off), heard = sum(y_on), live = any(live)),
           by = c("animal_id", "bin_start")]
  per <- per[per$detections > 0 & per$heard == 0 & per$live]
  table <- data.frame(animal_id = per$animal_id, bin_start = per$bin_start,
                      detections = as.integer(per$detections), stringsAsFactors = FALSE)
  if (nrow(table)) {
    warning(sprintf(paste0("%s had detections only at receivers not listening (%s); ",
                           "without them %s would read as unheard, so %s left out. ",
                           "See attr(, \"not_listening_bins\")."),
                    plural(nrow(table), "bin"), plural(sum(table$detections), "detection"),
                    agree(nrow(table), "the bin", "the bins"),
                    agree(nrow(table), "it is", "they are")),
            call. = FALSE)
  }
  list(table = table, key = paste(per$animal_id, as.numeric(per$bin_start)))
}

# The floor n >= max y also holds for counts built by hand.
ls_floor <- function(rows) {
  n <- y <- n_bins <- bin_end <- NULL
  span <- as.numeric(rows$bin_end) - as.numeric(rows$bin_start)
  if (anyNA(span) || any(span <= 0)) stop("Every `bin_end` must follow its `bin_start`.",
                                          call. = FALSE)
  by <- c("animal_id", "bin_start")
  n_bins <- rows[, list(k = data.table::uniqueN(n), e = data.table::uniqueN(bin_end),
                        m = max(y), n = n[1]), by = by]
  for (what in c("k", "e")) {
    if (any(n_bins[[what]] > 1)) {
      bad <- n_bins[n_bins[[what]] > 1][1]
      stop(sprintf("Animal %s, bin starting %s UTC: rows give different `%s`.", bad$animal_id,
                   ls_time(bad$bin_start), if (what == "k") "n" else "bin_end"), call. = FALSE)
    }
  }
  raised <- n_bins$m > n_bins$n
  if (any(raised)) {
    message(wrap_text(sprintf("n was raised to the largest receiver count in %s of %s bins.",
                              fmt_n(sum(raised)),
                              fmt_n(nrow(n_bins)))))
    rows[, n := pmax(n, max(y)), by = by]
  }
  if (any(n_bins$n == 0 & n_bins$m == 0)) {
    bad <- n_bins[n_bins$n == 0 & n_bins$m == 0][1]
    stop(sprintf(paste0("Animal %s, bin starting %s UTC: n is 0, so the counts say nothing ",
                        "about where the animal was. Check the delays given to ",
                        "bin_detections()."), bad$animal_id, ls_time(bad$bin_start)),
         call. = FALSE)
  }
  rows
}

# Two rows of one bin on one field column would count that receiver twice.
ls_check_duplicates <- function(rows, key) {
  dup <- duplicated(rows[, c("animal_id", "bin_start", key), with = FALSE])
  if (any(dup)) {
    stop(sprintf("Animal %s, bin starting %s UTC: %s %s appears twice.",
                 rows$animal_id[dup][1], ls_time(rows$bin_start[dup][1]), key,
                 rows[[key]][dup][1]), call. = FALSE)
  }
  invisible(rows)
}

ls_time <- function(t) format(t, "%Y-%m-%d %H:%M:%S", tz = "UTC")

# ---- Field ----
ls_field <- function(field, rows) {
  if (!is.list(field)) stop("`field` must be a detection field (a list).", call. = FALSE)
  # Reads go through unclass(), because a de_field's own `$` method computes
  # elements such as `p` on the fly.
  f <- unclass(field)
  eta <- ls_field_eta(f)
  receivers <- if (is.null(f$receivers)) colnames(eta) else f$receivers
  receivers <- if (is.data.frame(receivers)) receivers[[1]] else as.character(receivers)
  cells <- f$cells
  if (length(receivers) != ncol(eta) || length(cells) != nrow(eta) || anyNA(eta)) {
    stop("`field$eta` must have one row per cell in `field$cells`, one column per receiver ",
         "in `field$receivers`, and no NA.", call. = FALSE)
  }
  grid <- ls_read_grid(f$grid)
  ls_check_cells(grid, cells)
  key <- ls_key(f, rows, receivers)
  model <- if (is.null(f$model)) list() else unclass(f$model)
  check_station_day(model$sigma$station_day)
  xy <- ls_receiver_xy(f, key)
  near <- ls_near_setup(grid, cells, xy, f, eta, receivers, key)
  # terra warns when it packs a template grid that holds no values.
  if (!is.null(grid)) grid <- suppressWarnings(terra::wrap(grid))
  list(eta = eta, receivers = receivers, cells = as.integer(cells), grid = grid,
       sigma = model$sigma$station_hour, key = key, xy = xy, near = near)
}

# Cell centres, receiver positions and each receiver's d50 in the field, for
# the distance columns of the bin table. NULL when the field holds
# no grid or no receiver positions, or the grid is in degrees.
ls_near_setup <- function(grid, cells, xy, f, eta, receivers, key) {
  if (is.null(grid) || is.null(xy) || terra::is.lonlat(grid)) return(NULL)
  cxy <- terra::xyFromCell(grid, cells)
  list(receiver = as.character(xy$receiver), x = xy$x, y = xy$y, cx = cxy[, 1],
       cy = cxy[, 2], d50 = ls_d50(f$distances, eta, cells, receivers, key,
                                   as.character(xy$receiver)))
}

# The longest clear-line distance at which the field gives a receiver a
# detection probability of 0.5 or more (logit 0), per receiver; NA where the
# field has no distance table or the receiver no such pair.
ls_d50 <- function(dist, eta, cells, receivers, key, ids) {
  out <- rep(NA_real_, length(ids))
  need <- c("cell", key, "distance_m", "in_water")
  if (!is.data.frame(dist) || !all(need %in% names(dist))) return(out)
  row <- match(dist$cell, cells)
  col <- match(as.character(dist[[key]]), receivers)
  ok <- !is.na(row) & !is.na(col) & dist$in_water %in% TRUE
  e <- eta[cbind(row[ok], col[ok])]
  half <- e >= 0
  if (!any(half)) return(out)
  d50 <- tapply(dist$distance_m[ok][half], receivers[col[ok][half]], max)
  out[ids %in% names(d50)] <- as.numeric(d50[ids[ids %in% names(d50)]])
  out
}

# The field stores logits; a field that holds probabilities is converted.
ls_field_eta <- function(f) {
  eta <- f$eta
  if (is.null(eta) && !is.null(f$p)) eta <- stats::qlogis(f$p)
  if (!is.matrix(eta) || !is.numeric(eta)) {
    stop("`field` must hold `eta`, a numeric matrix of cells by receivers.", call. = FALSE)
  }
  eta
}

# A field cell on land would take probability that no animal can hold. A
# template grid without values cannot be checked; after terra::wrap() such a
# grid reads back as all NA, which no water grid is, so that counts as none.
ls_check_cells <- function(grid, cells) {
  if (is.null(grid) || !terra::hasValues(grid)) return(invisible(NULL))
  v <- suppressWarnings(terra::values(grid[[1]], mat = FALSE))
  if (all(is.na(v))) return(invisible(NULL))
  inside <- !is.na(cells) & cells >= 1 & cells <= length(v)
  bad <- !inside
  bad[inside] <- is.na(v[cells[inside]])
  if (any(bad)) {
    stop(sprintf("%s of the field's %s cells are not water cells of its grid (NA or outside ",
                 fmt_n(sum(bad)), fmt_n(length(cells))), sprintf("it): %s.", bd_first(cells[bad])),
         call. = FALSE)
  }
  invisible(NULL)
}

# Counts binned against a listening table keyed by receiver_serial carry no
# station_id. A field keyed by station_id that holds the deployment lookup
# (a field built from deployments or from an arraydesign layout) maps each
# serial to its station on the bin's date, as system_de() does.
ls_station_from_serial <- function(rows, field) {
  if (!is.list(field)) return(rows)
  f <- unclass(field)
  keyed_by_serial <- identical(f$key, "station_id") && !"station_id" %in% names(rows) &&
    "receiver_serial" %in% names(rows)
  if (!keyed_by_serial || is.null(f$lookup)) return(rows)
  station <- serial_to_station(f$lookup, rows$receiver_serial,
                               as.Date(rows$bin_start, tz = "UTC"))
  if (anyNA(station)) {
    lost <- unique(rows$receiver_serial[is.na(station)])
    stop(sprintf(paste0("The field is keyed by station_id and `counts` has none; %s ",
                        "%s no deployment in the field's lookup on the bin's date: %s."),
                 plural(length(lost), "receiver_serial value"), agree(length(lost), "has", "have"),
                 bd_first(lost)), call. = FALSE)
  }
  rows <- data.table::copy(rows)
  rows[, station_id := station]
  rows
}

ls_key <- function(f, rows, receivers) {
  keys <- intersect(c("receiver_serial", "station_id"), names(rows))
  if (is.character(f$key) && length(f$key) == 1) {
    if (!f$key %in% names(rows)) {
      stop(sprintf("The field's columns are %s, a column `counts` lacks.", f$key),
           call. = FALSE)
    }
    keys <- f$key
  }
  missing <- lapply(keys, function(k) unique(rows[[k]][!rows[[k]] %in% receivers]))
  names(missing) <- keys
  ok <- keys[lengths(missing) == 0]
  if (length(ok)) return(ok[1])
  best <- keys[which.min(lengths(missing))]
  n_miss <- length(missing[[best]])
  stop(sprintf("%s in `counts` %s no column in `field` (matched by %s): %s.",
               plural(n_miss, "listening receiver"), agree(n_miss, "has", "have"), best,
               bd_first(missing[[best]])), call. = FALSE)
}

# Receiver coordinates, from any table in the field with the key, x and y.
ls_receiver_xy <- function(f, key) {
  for (el in f) {
    if (is.data.frame(el) && all(c(key, "x", "y") %in% names(el)) && nrow(el) < 1e5) {
      xy <- unique(data.frame(receiver = as.character(el[[key]]), x = el$x, y = el$y))
      return(xy[!duplicated(xy$receiver), ])
    }
  }
  NULL
}

# A surface carries its grid. A path would depend on the working directory,
# so a path is read here and the raster is stored packed (a
# SpatRaster holds a pointer that does not survive saveRDS() or a fork).
ls_read_grid <- function(grid) {
  if (is.null(grid)) return(NULL)
  if (inherits(grid, "PackedSpatRaster")) return(terra::rast(grid))
  if (is.character(grid) && length(grid) == 1) {
    if (!file.exists(grid)) {
      stop(sprintf("The field's grid file %s does not exist (working directory %s).", grid,
                   getwd()), call. = FALSE)
    }
    grid <- terra::rast(normalizePath(grid))
  }
  if (!inherits(grid, "SpatRaster")) {
    stop("`field$grid` must be a SpatRaster, a packed SpatRaster or a raster file.",
         call. = FALSE)
  }
  grid
}

# ---- Reading of the curve ----
# The presence-only reading has no term for misses, so no spread to integrate:
# it runs as "none", and says so only when "station_hour" was asked for.
ls_reading <- function(likelihood, overdispersion, asked) {
  if (likelihood == "binomial") return(overdispersion)
  if (asked && overdispersion == "station_hour") {
    warning(paste0("overdispersion = \"station_hour\" is ignored under likelihood = ",
                   "\"presence_only\", which does not model misses."), call. = FALSE)
  }
  "none"
}

# ---- Overdispersion ----
ls_overdispersion <- function(form, sigma, ctrl, rows) {
  none <- list(form = "none", sigma = 0, nodes = NA_integer_)
  if (form == "none") return(none)
  # A spread of 0 is the binomial, and is labelled so.
  if (identical(sigma, 0) || identical(sigma, 0L)) return(none)
  if (!bd_positive(sigma)) {
    # The code line is kept whole, outside the wrapped text, so it can be
    # copied as it stands.
    stop(paste0(
      wrap_text(paste0("The field's detection model gives no `sigma$station_hour`, which ",
                       "overdispersion = \"station_hour\", the default, needs. For plain ",
                       "binomial terms, build the field with a spread of 0:")),
      "\n  de_field(distances, de_model(model, sigma = list(station_hour = 0)))\n",
      wrap_text(paste0("Give a spread estimated around the same curve in place of 0, or ",
                       "use overdispersion = \"none\" for plain binomial terms."))
    ), call. = FALSE)
  }
  span <- unique(as.numeric(rows$bin_end) - as.numeric(rows$bin_start))
  if (any(span != 3600)) {
    msg <- paste0("sigma$station_hour describes hour bins; these bins last %s s. ",
                  "The overdispersion may be wrong for them.")
    warning(sprintf(msg, paste(utils::head(span, 3), collapse = ", ")), call. = FALSE)
  }
  nodes <- ls_check_nodes(sigma, ctrl$nodes)
  gh <- gauss_hermite(nodes)
  list(form = "station_hour", sigma = sigma, nodes = nodes, z = gh$z, w = gh$w,
       step = ctrl$table_step)
}

# Nodes for a standard deviation. Against integrate(), one receiver's term
# is off by at most 1.2e-6 in log probability at sigma 2.55 (64 nodes),
# 6.0e-6 at 3.4 (85) and 1.1e-4 at 6 (150), over counts 0 to n for n of 30,
# 60 and 720; 20 nodes are off by 3.2e-4, 1.8e-3 and 1.8e-2
# (<https://github.com/jakebrownscombe/detfield/blob/main/validation/02_quadrature_error.R>).
ls_nodes <- function(sigma) {
  as.integer(min(200, max(20, ceiling(25 * sigma))))
}

# Fewer nodes than the rule, or a spread whose rule passes the cap, gives a
# surface off by more than the rule allows, so either warns.
ls_check_nodes <- function(sigma, nodes) {
  rule <- ls_nodes(sigma)
  if (25 * sigma > 200) {
    warning(sprintf(paste0("sigma %g needs about %s Gauss-Hermite nodes and 200 is the most ",
                           "allowed; the surfaces may be off by more than 1e-4."), sigma,
                    fmt_n(as.integer(ceiling(25 * sigma)))), call. = FALSE)
  }
  if (is.null(nodes)) return(rule)
  if (nodes < rule) {
    warning(sprintf(paste0("control$nodes = %s is below the %s nodes the default gives at ",
                           "sigma %g; the surfaces may be off (see ?location_surface, ",
                           "Computation)."), fmt_n(nodes), fmt_n(rule), sigma), call. = FALSE)
  }
  nodes
}

# Gauss-Hermite nodes and weights by the Golub-Welsch eigenvalue method
# (Golub and Welsch 1969, Mathematics of Computation 23, 221-230,
# doi:10.1090/S0025-5718-69-99647-1): the nodes are the eigenvalues of the
# symmetric tridiagonal Jacobi matrix and the weights the squared first
# components of its eigenvectors, scaled so that sum(w * f(sqrt(2) * z)) is
# the mean of f over a standard normal.
gauss_hermite <- function(nodes) {
  if (nodes == 1) return(list(z = 0, w = 1))
  k <- seq_len(nodes - 1)
  jac <- matrix(0, nodes, nodes)
  jac[cbind(k, k + 1)] <- sqrt(k / 2)
  jac[cbind(k + 1, k)] <- sqrt(k / 2)
  e <- eigen(jac, symmetric = TRUE)
  w <- e$vectors[1, ]^2
  list(z = e$values, w = w / sum(w))
}

# ---- One animal ----
ls_animal <- function(a, fld, od, ctrl) {
  cols <- sort(unique(a$col))
  a <- data.table::copy(a)
  j <- NULL
  a[, j := match(col, cols)]
  starts <- sort(unique(a$bin_start))
  idx <- split(seq_len(nrow(a)), match(a$bin_start, starts))
  n_bin <- vapply(idx, function(r) a$n[r[1]], numeric(1))
  prep <- ls_prepare(fld$eta[, cols, drop = FALSE], od, n_bin)
  out <- vector("list", length(idx))
  for (i in seq_along(idx)) {
    out[[i]] <- ls_bin(a[idx[[i]]], n_bin[i], prep, fld, ctrl$kept_mass)
    if (!is.null(out[[i]]$why) && ctrl$on_degenerate == "stop") {
      ls_stop_degenerate(out[[i]]$why, a$animal_id[1], starts[i])
    }
  }
  ls_animal_tables(a, starts, idx, n_bin, out)
}

# One bin: its surface, what is kept, and the numbers a reader checks it by;
# or, for a degenerate surface, the reason.
ls_bin <- function(b, n, prep, fld, kept_mass) {
  ll <- ls_loglik(b$j, b$y, n, prep)
  why <- ls_degenerate(ll)
  if (!is.null(why)) return(list(why = why))
  top <- max(ll)
  w <- exp(ll - top)
  total <- sum(w)
  prob <- w / total
  keep <- ls_keep(prob, kept_mass)
  fit <- ls_misfit(b, n, prob, prep, fld$key)
  info <- list(n_listening = nrow(b), n_detecting = sum(b$y > 0), y_total = sum(b$y),
               log_lik = top + log(total) - log(length(ll)), kept_mass = sum(prob[keep]),
               n_cells = length(keep), misfit_p = fit$p, misfit_receiver = fit$receiver)
  list(cell = fld$cells[keep], prob = prob[keep], info = c(info, ls_near(b, prob, fld)))
}

ls_animal_tables <- function(a, starts, idx, n_bin, out) {
  ok <- vapply(out, function(o) is.null(o$why), logical(1))
  skipped <- data.frame(animal_id = rep(a$animal_id[1], sum(!ok)), bin_start = starts[!ok],
                        reason = vapply(out[!ok], `[[`, character(1), "why"),
                        stringsAsFactors = FALSE)
  if (!any(ok)) return(list(surface = NULL, bins = NULL, skipped = skipped))
  out <- out[ok]
  k <- vapply(out, function(x) length(x$cell), integer(1))
  surface <- data.table::data.table(
    animal_id = a$animal_id[1], bin_start = rep(starts[ok], k),
    cell = unlist(lapply(out, `[[`, "cell"), use.names = FALSE),
    prob = unlist(lapply(out, `[[`, "prob"), use.names = FALSE)
  )
  first <- vapply(idx[ok], `[`, integer(1), 1)
  bins <- cbind(data.table::data.table(animal_id = a$animal_id[1], bin_start = starts[ok],
                                       bin_end = a$bin_end[first], n_transmissions = n_bin[ok]),
                data.table::rbindlist(lapply(out, `[[`, "info")))
  list(surface = surface, bins = bins, skipped = skipped)
}

# A surface that is not a distribution is never written as a flat or empty
# surface: it stops the run or, when asked, is left out and recorded.
ls_degenerate <- function(ll) {
  if (anyNA(ll)) {
    "the log-likelihood is NaN in some cells"
  } else if (all(ll == -Inf)) {
    "every cell has zero likelihood"
  } else if (length(ll) > 1 && max(ll) == min(ll)) {
    "every cell has the same likelihood"
  }
}

ls_stop_degenerate <- function(why, animal, start) {
  msg <- paste0("Animal %s, bin starting %s UTC: the surface is degenerate, because %s. ",
                "Check the listening receivers and the field for this bin, or set ",
                "control$on_degenerate = \"skip\".")
  detfield_abort(sprintf(msg, animal, ls_time(start), why), "degenerate_bin")
}

# The smallest set of cells holding `mass`; ties go to the lower row. The
# probabilities sum to 1 within rounding, far inside the 1e-12 allowance, so
# a mass of 1 keeps every cell.
ls_keep <- function(prob, mass) {
  o <- order(prob, decreasing = TRUE, method = "radix")
  o[seq_len(which(cumsum(prob[o]) >= mass * (1 - 1e-12))[1])]
}

# ---- Distance from the detecting receivers ----
# The mode's straight-line distance to the nearest detecting receiver, and
# the mass within near_m (or the field's d50) of the receiver with the most
# detections. A ring of cells away from the receivers fits low counts, so
# misfit_p alone does not show it.
ls_near <- function(b, prob, fld) {
  mode <- which.max(prob)
  out <- list(mode_cell = fld$cells[mode], mode_dist_m = NA_real_,
              near_receiver = NA_character_, near_m = NA_real_, near_mass = NA_real_)
  det <- which(b$y > 0)
  if (!length(det)) return(out)
  top <- det[which.max(b$y[det])]
  out$near_receiver <- as.character(b[[fld$key]][top])
  nr <- fld$near
  if (is.null(nr)) return(out)
  k <- match(as.character(b[[fld$key]][det]), nr$receiver)
  k <- k[!is.na(k)]
  if (length(k)) {
    out$mode_dist_m <- min(sqrt((nr$x[k] - nr$cx[mode])^2 + (nr$y[k] - nr$cy[mode])^2))
  }
  kt <- match(out$near_receiver, nr$receiver)
  if (is.na(kt)) return(out)
  r <- if (is.null(nr$near_m)) nr$d50[kt] else nr$near_m
  if (is.na(r)) return(out)
  out$near_m <- r
  out$near_mass <- sum(prob[(nr$cx - nr$x[kt])^2 + (nr$cy - nr$y[kt])^2 <= r^2])
  out
}

# ---- Misfit ----
# Posterior predictive probability that each detecting receiver logs at
# least its count; the smallest, and the receiver that gives it.
ls_misfit <- function(b, n, prob, prep, key) {
  det <- which(b$y > 0)
  if (!length(det)) return(list(p = NA_real_, receiver = NA_character_))
  tails <- vapply(det, function(k) sum(prob * ls_upper(prep$eta[, b$j[k]], b$y[k], n, prep)),
                  numeric(1))
  i <- which.min(tails)
  list(p = min(1, tails[i]), receiver = as.character(b[[key]][det[i]]))
}

# P(Y >= y) for each logit, with Y binomial given the station-hour effect.
# Under the station-hour form it is read from a table over the logit by
# linear interpolation; below the table it is 0 to far better than 1e-12.
ls_upper <- function(eta, y, n, prep) {
  if (prep$form == "binomial") {
    return(stats::pbinom(y - 1, n, stats::plogis(eta), lower.tail = FALSE))
  }
  key <- paste(y, n)
  tab <- get0(key, envir = prep$tails, inherits = FALSE)
  if (is.null(tab)) {
    x <- seq(ls_lo(prep, y), prep$hi, by = prep$step)
    u <- outer(x, prep$od$sigma * sqrt(2) * prep$od$z, "+")
    up <- matrix(stats::pbinom(y - 1, n, stats::plogis(u), lower.tail = FALSE), nrow(u))
    tab <- list(x = x, v = drop(up %*% prep$od$w))
    assign(key, tab, envir = prep$tails)
  }
  out <- as.numeric(eta > 0)
  fin <- is.finite(eta)
  out[fin] <- stats::approx(tab$x, tab$v, eta[fin], rule = 2)$y
  out
}

# ---- Log-likelihood of one bin ----
ls_prepare <- function(eta, od, n_bin) {
  if (od$sigma == 0) {
    return(list(form = "binomial", eta = eta, l1 = stats::plogis(eta, log.p = TRUE),
                l0 = stats::plogis(-eta, log.p = TRUE),
                presence_only = identical(od$likelihood, "presence_only")))
  }
  fin <- eta[is.finite(eta)]
  if (!length(fin)) fin <- 0
  step <- od$step
  lo <- floor(min(fin) / step) * step - step
  hi <- ceiling(max(fin) / step) * step + step
  prep <- list(form = "quadrature", eta = eta, od = od, lo = lo, hi = hi, step = step,
               tables = new.env(), silent = new.env(), tails = new.env())
  # Silent receivers dominate the cost, so their terms are computed once for
  # every number of transmissions that two or more bins share.
  shared <- as.numeric(names(which(table(n_bin) >= 2)))
  for (n in shared) {
    assign(as.character(n), ls_term(eta, 0, n, prep), envir = prep$silent)
  }
  prep
}

ls_loglik <- function(j, y, n, prep) {
  if (isTRUE(prep$presence_only)) {
    # Each detection weighs its cell by p; misses and silent receivers add nothing.
    det <- y > 0
    return(drop(prep$l1[, j[det], drop = FALSE] %*% y[det]))
  }
  if (prep$form == "binomial") {
    det <- y > 0
    rest <- n - y > 0
    ll <- prep$l1[, j[det], drop = FALSE] %*% y[det] +
      prep$l0[, j[rest], drop = FALSE] %*% (n - y)[rest]
    return(drop(ll) + sum(lchoose(n, y)))
  }
  silent <- y == 0
  pre <- get0(as.character(n), envir = prep$silent, inherits = FALSE)
  ll <- if (!is.null(pre)) {
    rowSums(pre[, j[silent], drop = FALSE])
  } else {
    rowSums(ls_term(prep$eta[, j[silent], drop = FALSE], 0, n, prep))
  }
  for (k in which(!silent)) ll <- ll + ls_term(prep$eta[, j[k]], y[k], n, prep)
  ll
}

# ---- Receiver terms under the station-hour effect ----
# log E[ Binomial(y; n, plogis(eta + sigma * Z)) ] for each value of eta, with Z
# standard normal: from the table for (y, n), or in closed form below it.
ls_term <- function(eta, y, n, prep) {
  key <- paste(y, n)
  tab <- get0(key, envir = prep$tables, inherits = FALSE)
  if (is.null(tab)) {
    lo <- ls_lo(prep, y)
    tab <- ls_table(seq(lo, prep$hi, by = prep$step), y, n, prep$od)
    tab$lo <- lo
    assign(key, tab, envir = prep$tables)
  }
  inside <- eta >= tab$lo & eta < Inf
  if (all(inside)) return(ls_hermite(eta, tab, tab$lo, prep$step))
  out <- eta
  out[inside] <- ls_hermite(eta[inside], tab, tab$lo, prep$step)
  low <- eta < tab$lo
  out[low] <- ls_low(eta[low], y, n, prep$od$sigma)
  out[eta == Inf] <- if (y == n) 0 else -Inf
  out
}

# Lowest logit of the table for count y. Below it, under the normal tilted by
# the binomial, eta + sigma * Z stays under -30 except with probability below
# 1e-15, so p is below 1e-13 and the closed form of ls_low() holds; without
# this bound the table spans every logit of the field, which for one cell at
# -10,000 takes seconds and gigabytes.
ls_lo <- function(prep, y) {
  s <- prep$od$sigma
  cut <- floor((-30 - 8 * s - (1 + y) * s^2) / prep$step) * prep$step
  min(max(prep$lo, cut), prep$hi - prep$step)
}

# The term where p = plogis(eta + s * Z) is tiny: log E[(1 - p)^n] is
# -n E[p] for y = 0, and log choose(n, y) E[p^y] for y > 0, with E[p^y] =
# exp(y * eta + y^2 * s^2 / 2) for p = exp(eta + s * Z).
ls_low <- function(eta, y, n, s) {
  if (y == 0) return(-n * exp(eta + s^2 / 2))
  lchoose(n, y) + y * eta + y^2 * s^2 / 2
}

# The term and its derivative at each grid value, by adaptive Gauss-Hermite
# quadrature in log space: the nodes sit on the mode of each integrand and
# scale with its curvature, because with many transmissions the binomial is
# too narrow in the logit for fixed nodes to find. The derivative of
# log Binomial(y; n, p) in the logit is y - n * p.
ls_table <- function(grid, y, n, od) {
  s <- od$sigma
  mode <- ls_mode(grid, y, n, s)
  p <- stats::plogis(grid + s * mode)
  tau <- 1 / sqrt(s^2 * n * p * (1 - p) + 1)
  z <- mode + sqrt(2) * outer(tau, od$z)
  u <- grid + s * z
  lw <- lchoose(n, y) + y * stats::plogis(u, log.p = TRUE) +
    (n - y) * stats::plogis(-u, log.p = TRUE) - z^2 / 2 +
    rep(od$z^2 + log(od$w), each = length(grid))
  top <- lw[cbind(seq_along(grid), max.col(lw, ties.method = "first"))]
  e <- exp(lw - top)
  total <- rowSums(e)
  ls_cubic(log(tau) + top + log(total), rowSums(e * (y - n * stats::plogis(u))) / total,
           grid[2] - grid[1])
}

# Coefficients of the cubic Hermite polynomial on each grid interval, in the
# interval's own coordinate t from 0 to 1, so evaluation is one Horner step.
ls_cubic <- function(g, dg, step) {
  k <- length(g)
  g0 <- g[-k]
  g1 <- g[-1]
  d0 <- step * dg[-k]
  d1 <- step * dg[-1]
  list(c0 = g0, c1 = d0, c2 = 3 * (g1 - g0) - 2 * d0 - d1, c3 = 2 * (g0 - g1) + d0 + d1)
}

# Mode in z of log Binomial(y; n, plogis(eta + s * z)) - z^2 / 2. Its slope
# falls with z and has one root, between s * (y - n) and s * y. Plain Newton
# can cycle there, so the bracket is first halved 64 times, which brings it
# below 1e-15 of its width, and three Newton steps kept inside it polish the
# root.
ls_mode <- function(eta, y, n, s) {
  slope <- function(z) s * (y - n * stats::plogis(eta + s * z)) - z
  lo <- rep(s * (y - n), length(eta))
  hi <- rep(s * y, length(eta))
  for (it in seq_len(64)) {
    mid <- (lo + hi) / 2
    up <- slope(mid) > 0
    lo <- ifelse(up, mid, lo)
    hi <- ifelse(up, hi, mid)
  }
  z <- (lo + hi) / 2
  for (it in 1:3) {
    p <- stats::plogis(eta + s * z)
    z <- pmin(pmax(z + slope(z) / (s^2 * n * p * (1 - p) + 1), lo), hi)
  }
  z
}

# Cubic Hermite interpolation on a uniform grid that starts at lo. The grid
# extends one step beyond the field's logits, so every finite logit at or
# above lo falls inside it.
ls_hermite <- function(eta, tab, lo, step) {
  pos <- (eta - lo) / step
  i <- as.integer(pos)
  t <- pos - i
  i <- i + 1L
  tab$c0[i] + t * (tab$c1[i] + t * (tab$c2[i] + t * tab$c3[i]))
}

# ---- Assemble the result ----
# The most rows one data.table holds.
ls_max_rows <- function() .Machine$integer.max

# Before any bin is computed: when every water cell of every bin could not
# fit in one table, say how to avoid a stop at the end of the run.
ls_warn_size <- function(rows, fld, ctrl) {
  bins <- nrow(unique(rows[, c("animal_id", "bin_start")]))
  most <- as.numeric(bins) * nrow(fld$eta)
  if (is.null(ctrl$write_dir) && most > ls_max_rows()) {
    message(wrap_text(sprintf(
      paste0("%s bins x %s water cells could keep more rows than one table ",
             "holds (%s). If the run stops at the end for this reason, set ",
             "control$write_dir or give fewer animals per call."),
      fmt_n(bins), fmt_n(nrow(fld$eta)), fmt_n(ls_max_rows())
    )))
  }
  invisible(most)
}

# One table holds at most ls_max_rows() rows; a run past it stops here, with
# the remedy, rather than inside data.table.
ls_check_rows <- function(total, limit = ls_max_rows()) {
  if (total > limit) {
    stop(sprintf(paste0("The kept surfaces hold %s rows, more than one table can hold (%s). ",
                        "Set control$write_dir to write each animal's surfaces to its own ",
                        "file, or give fewer animals per call."),
                 fmt_n(total), fmt_n(limit)), call. = FALSE)
  }
  invisible(total)
}

# One animal, in a worker. With write_dir the surfaces are saved there and
# only the bin table goes back to the session.
ls_animal_out <- function(a, fld, od, ctrl) {
  r <- ls_animal(a, fld, od, ctrl)
  if (is.null(ctrl$write_dir)) return(r)
  file <- NA_character_
  if (!is.null(r$surface)) {
    file <- file.path(ctrl$write_dir, ls_file_name(a$animal_id[1]))
    saveRDS(ls_structure(r$surface, r$bins, r$skipped, a, fld, od, ctrl, NULL), file)
  }
  list(surface = NULL, bins = r$bins, skipped = r$skipped,
       file = data.frame(animal_id = a$animal_id[1], file = file, bins = NROW(r$bins),
                         rows = NROW(r$surface), stringsAsFactors = FALSE))
}

# A file name from an animal_id: characters a file system may refuse become
# "_", and a short hash of the id keeps two ids that differ only there apart.
ls_file_name <- function(id) {
  hash <- sum(utf8ToInt(id) * seq_along(utf8ToInt(id))) %% 1e6
  sprintf("%s_%06d.rds", gsub("[^A-Za-z0-9._-]", "_", id), as.integer(hash))
}

ls_files <- function(res, od, lost) {
  bins <- data.table::rbindlist(lapply(res, `[[`, "bins"))
  skipped <- ls_skipped(res, nrow(bins))
  out <- do.call(rbind, lapply(res, `[[`, "file"))
  rownames(out) <- NULL
  structure(out, bins = as.data.frame(bins), overdispersion = od[c("form", "sigma", "nodes")],
            likelihood = od$likelihood, skipped = skipped, not_listening_bins = lost)
}

# Degenerate bins left out: a warning with the count, or a stop when no bin
# is left.
ls_skipped <- function(res, n_bins) {
  skipped <- do.call(rbind, lapply(res, `[[`, "skipped"))
  if (nrow(skipped)) {
    if (!n_bins) {
      msg <- "All %s bins are degenerate; the first: animal %s, bin starting %s UTC, %s."
      detfield_abort(sprintf(msg, fmt_n(nrow(skipped)), skipped$animal_id[1],
                             ls_time(skipped$bin_start[1]), skipped$reason[1]), "degenerate")
    }
    warning(sprintf("%s of %s bins were degenerate and left out; see attr(, \"skipped\").",
                    fmt_n(nrow(skipped)),
                    fmt_n(nrow(skipped) + n_bins)), call. = FALSE)
  }
  skipped
}

ls_assemble <- function(res, rows, fld, od, ctrl, lost = NULL) {
  ls_check_rows(sum(vapply(res, function(r) as.numeric(NROW(r$surface)), numeric(1))))
  surface <- data.table::rbindlist(lapply(res, `[[`, "surface"))
  bins <- data.table::rbindlist(lapply(res, `[[`, "bins"))
  skipped <- ls_skipped(res, nrow(bins))
  ls_structure(surface, bins, skipped, rows, fld, od, ctrl, lost)
}

ls_structure <- function(surface, bins, skipped, rows, fld, od, ctrl, lost) {
  out <- data.frame(animal_id = surface$animal_id, bin_start = surface$bin_start,
                    cell = surface$cell, prob = surface$prob, stringsAsFactors = FALSE)
  used <- data.frame(animal_id = rows$animal_id, bin_start = rows$bin_start,
                     receiver = rows[[fld$key]], y = as.integer(rows$y),
                     stringsAsFactors = FALSE)
  structure(
    out,
    class = c("loc_surface", "data.frame"),
    bins = as.data.frame(bins),
    overdispersion = od[c("form", "sigma", "nodes")],
    likelihood = od$likelihood,
    kept = list(mass = ctrl$kept_mass,
                rule = "smallest set of cells holding the mass; probabilities not rescaled"),
    grid = fld$grid,
    n_water_cells = nrow(fld$eta),
    key = fld$key,
    receivers = used,
    receiver_xy = fld$xy,
    near_rule = if (is.null(ctrl$near_m)) "d50" else "control$near_m",
    skipped = skipped,
    not_listening_bins = lost
  )
}
