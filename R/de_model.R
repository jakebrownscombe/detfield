#' A detection model for positioning
#'
#' Wraps a range model so that the rest of the package can ask one question
#' of it: the probability, on the logit scale, that one transmission at a
#' given distance from a receiver is detected. The model is used as given:
#' `de_model()` never fits anything, and the between-station overdispersion
#' and the time-varying modulators are arguments the user gives.
#'
#' `x` can be:
#' * a fitted `glm` (including an mgcv `gam`) or `glmmTMB` model, called
#'   through [stats::predict()] on the link scale. The data frame it receives
#'   holds `distance_m`, `depth_m`, `tag_model` and one column per modulator;
#'   give `covariates` when the model was fitted on other names or units. For
#'   `glmmTMB` the random effects are set to zero (`re.form = NA`) unless
#'   `...` says otherwise. glmmTMB fits the wiggly part of an `s()` smooth
#'   as a random effect, so `re.form = NA` drops it and a glmmTMB smooth over
#'   time reaches the field as a straight line. `de_model()` warns when a
#'   glmmTMB formula holds an `s()`, `te()`, `t2()` or `ti()` term and `...`
#'   gives no `re.form`. Fit a smooth over time in mgcv (`gam()` or `bam()`)
#'   and pass its random-effect smooths to `exclude` through `...`; or pass
#'   `re.form = NULL`, which keeps the smooth and also adds the fitted random
#'   effect of any grouping column the curve is given.
#' * a data frame of curve values with `distance_m` and either `eta` (logit)
#'   or `p` (probability), and optionally `tag_model` and `depth_m`. Values
#'   between rows are interpolated linearly on the logit scale, in distance
#'   and then in depth. Beyond the last distance the logit continues along
#'   the slope of the last two rows when that slope is negative and is held
#'   otherwise; depth is held at the nearest tabulated depth.
#' * a function of `distance_m`, `depth_m`, `tag_model` and `...` that
#'   returns the logit. Modulator values arrive through `...` by name.
#' * a `de_model`, to replace its `sigma` or `modulators`.
#'
#' @param x A fitted model, a data frame of curve values, a function or a
#'   `de_model`; see Details.
#' @param sigma A list of standard deviations, in logits. `station_hour` is
#'   the whole spread of one receiver's logit about this curve within one
#'   hour bin. It is defined only together with the curve it was estimated
#'   around: estimate both from the same bin-scale data, and do not pair a
#'   value estimated around one curve with another (see
#'   [location_surface()], Overdispersion). `station_day` is reserved and
#'   must be `NA`; the likelihood has no day-scale term. `NA` means not given.
#'   The values are stored for the likelihood and do not change the curve.
#' @param modulators A named list of covariates that change with date. Each
#'   element is a constant, a vector named by receiver key, or a function of
#'   one date that returns either. [de_field()] evaluates them for its date and
#'   passes them to the model by name.
#' @param blocked A logit offset added to the curve for receiver and cell
#'   pairs whose straight line meets land (`in_water` `FALSE` in
#'   [water_distances()]). `0`, the default, reads the curve at the distance
#'   through water, as positionR 1.x did. `-Inf` makes a receiver deaf to the
#'   cells behind land. Range models fitted on clear lines of sight say
#'   nothing about these pairs, so the value is the user's assumption or an
#'   estimate made outside the package; it is never fitted here. Reference
#'   tags behind land can log a small share of the detections the clear-line
#'   curve predicts, so 0 can be far from right; [de_field()]
#'   counts the pairs it leaves at the clear-line value.
#' @param covariates For a fitted model, an optional function that takes the
#'   data frame described above and returns the `newdata` the model's
#'   `predict()` method needs.
#' @param ... For a fitted model, further arguments to `predict()`.
#'
#' @return An object of class `de_model`: a list with `predict`, a function
#'   of `distance_m`, `depth_m`, `tag_model` and `...` that returns the logit,
#'   one value per element of the longest of the three (an argument of length
#'   one is recycled; other lengths that differ stop);
#'   `sigma`; `modulators`; `blocked`; `tag_models`, the tag models a
#'   table or a fitted model knows, or `NULL`; `uses`, the arguments the curve is known to read,
#'   which [de_field()] checks before it calls the curve; and `source`, one
#'   line saying where the model came from. A function that takes `...`, or a
#'   fitted model given with `covariates`, declares only `distance_m`, because
#'   which other arguments it reads cannot be known in advance.
#'
#' @seealso [de_logistic()] for a curve from two distances; [de_field()].
#' @export
#' @examples
#' # The Stoney range test: hourly counts at receivers in 1 m and 25 m of water.
#' fit <- glm(cbind(heard, sent - heard) ~ distance_m * depth_m,
#'            family = binomial, data = stoney_range_test)
#' # The counts show no spread beyond binomial chance about this curve.
#' m <- de_model(fit, sigma = list(station_hour = 0))
#' m
#' plogis(m$predict(distance_m = c(0, 250, 500), depth_m = 1))
#'
#' # The same curve as a table of values, interpolated on the logit scale.
#' tab <- data.frame(distance_m = seq(0, 1600, by = 100), depth_m = 1)
#' tab$eta <- m$predict(tab$distance_m, depth_m = 1)
#' plogis(de_model(tab)$predict(c(0, 250, 500), depth_m = 1))
de_model <- function(x, sigma = list(station_hour = NA, station_day = NA), modulators = list(),
                     blocked = 0, covariates = NULL, ...) {
  sigma <- check_sigma(sigma)
  modulators <- check_modulators(modulators)
  if (!is_one_number(blocked) || blocked > 0) {
    stop("`blocked` must be one number of 0 or less, in logits; -Inf is allowed.", call. = FALSE)
  }
  if (inherits(x, "de_model")) {
    x$sigma <- sigma
    x$modulators <- modulators
    x$blocked <- blocked
    return(x)
  }
  core <- if (inherits(x, "glmmTMB") || inherits(x, "glm")) {
    model_curve(x, covariates, list(...))
  } else if (is.data.frame(x)) {
    table_curve(x)
  } else if (is.function(x)) {
    function_curve(x)
  } else {
    stop("`x` must be a fitted glm or glmmTMB model, a data frame of curve values, ",
         "a function or a de_model.", call. = FALSE)
  }
  new_de_model(core, sigma, modulators, blocked)
}

new_de_model <- function(core, sigma, modulators, blocked = 0) {
  structure(list(predict = recycled_curve(core$predict), sigma = sigma, modulators = modulators,
                 blocked = blocked, tag_models = core$tag_models, uses = core$uses,
                 source = core$source), class = "de_model")
}

# The curves take their length from `distance_m`, so one distance with four
# depths would give one value and no warning. Arguments of length one are
# recycled to the longest; any other length mismatch stops.
recycled_curve <- function(curve) {
  force(curve)
  function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
    given <- list(distance_m = distance_m, depth_m = depth_m, tag_model = tag_model)
    given <- given[lengths(given) > 0 | names(given) == "distance_m"]
    lens <- lengths(given)
    n <- max(lens)
    if (length(distance_m) && any(!lens %in% c(1L, n))) {
      got <- paste(sprintf("%s %s", names(lens), fmt_n(lens)), collapse = ", ")
      stop(sprintf(paste0("`distance_m`, `depth_m` and `tag_model` must have length 1 or a ",
                          "common length; got %s."), got), call. = FALSE)
    }
    if (length(distance_m) == 1L && n > 1L) distance_m <- rep(distance_m, n)
    if (length(depth_m) == 1L && n > 1L) depth_m <- rep(depth_m, n)
    if (length(tag_model) == 1L && n > 1L) tag_model <- rep(tag_model, n)
    curve(distance_m, depth_m, tag_model, ...)
  }
}

check_sigma <- function(sigma) {
  if (!is.list(sigma)) stop("`sigma` must be a list.", call. = FALSE)
  unknown <- setdiff(names(sigma), c("station_hour", "station_day"))
  if (length(unknown) || (length(sigma) && is.null(names(sigma)))) {
    stop("`sigma` takes only `station_hour` and `station_day`.", call. = FALSE)
  }
  out <- list(station_hour = NA_real_, station_day = NA_real_)
  for (nm in names(sigma)) out[[nm]] <- check_one_sigma(sigma[[nm]], nm)
  check_station_day(out$station_day)
  out
}

check_one_sigma <- function(v, nm) {
  ok <- length(v) == 1 && (is.na(v) || (is.numeric(v) && v >= 0))
  if (!ok) stop(sprintf("`sigma$%s` must be one number of at least 0, or NA.", nm), call. = FALSE)
  as.numeric(v)
}

check_modulators <- function(modulators) {
  if (!is.list(modulators)) stop("`modulators` must be a list.", call. = FALSE)
  if (length(modulators) && (is.null(names(modulators)) || any(!nzchar(names(modulators))))) {
    stop("Every element of `modulators` needs a name: the covariate it feeds.", call. = FALSE)
  }
  clash <- intersect(names(modulators), c("distance_m", "depth_m", "tag_model"))
  if (length(clash)) {
    stop(sprintf("Modulator %s %s %s reserved for the curve's own arguments.",
                 agree(length(clash), "name", "names"), paste(clash, collapse = ", "),
                 agree(length(clash), "is", "are")), call. = FALSE)
  }
  modulators
}

#' @export
print.de_model <- function(x, ...) {
  cat_wrapped(paste("<de_model>", x$source))
  cat("  reads:", paste(x$uses, collapse = ", "), "\n")
  s <- x$sigma$station_hour
  given <- if (is.na(s)) "not given" else plural(s, "logit")
  cat_wrapped(sprintf("  sigma: station_hour %s\n", given))
  cat_wrapped(sprintf("  blocked line of sight: %s\n", blocked_label(x$blocked)))
  if (length(x$modulators)) cat("  modulators:", paste(names(x$modulators), collapse = ", "), "\n")
  invisible(x)
}

blocked_label <- function(blocked) {
  if (is.null(blocked) || blocked == 0) return("curve read at the distance through water")
  if (is.infinite(blocked)) return("receiver deaf behind land")
  sprintf("curve plus %g logits", blocked)
}

# The likelihood takes one per-bin spread. A station-day value would be
# ignored without notice, so it stops instead.
check_station_day <- function(v) {
  if (!is.null(v) && length(v) == 1 && !is.na(v)) {
    stop(paste0("`sigma$station_day` is not used: location_surface() takes one spread per ",
                "receiver and bin, `sigma$station_hour`. Leave station_day NA."), call. = FALSE)
  }
  invisible(NULL)
}
