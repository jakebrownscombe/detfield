# ---- Curve builders behind de_model() ----
curve_args <- c("distance_m", "depth_m", "tag_model")

curve_frame <- function(distance_m, depth_m, tag_model, mods) {
  if (!is.numeric(distance_m)) stop("`distance_m` must be numeric.", call. = FALSE)
  n <- length(distance_m)
  d <- data.frame(distance_m = distance_m, depth_m = rep_len(as.numeric(depth_m), n),
                  tag_model = rep_len(as.character(tag_model), n), stringsAsFactors = FALSE)
  for (nm in names(mods)) d[[nm]] <- rep_len(mods[[nm]], n)
  d
}

check_eta <- function(eta, n) {
  eta <- as.numeric(eta)
  if (length(eta) != n) {
    stop(sprintf("The detection model returned %s for %s.", plural(length(eta), "value"),
                 plural(n, "distance")),
         call. = FALSE)
  }
  if (anyNA(eta)) {
    stop(sprintf("The detection model returned NA for %s of %s distances.",
                 fmt_n(sum(is.na(eta))), fmt_n(n)),
         call. = FALSE)
  }
  eta
}

model_curve <- function(fit, covariates, dots) {
  dots <- model_curve_dots(fit, covariates, dots)
  # glmmTMB needs the grouping columns in newdata even when random effects are
  # set to zero, so they are added as NA.
  groups <- if (inherits(fit, "glmmTMB")) unlist(fit$modelInfo$grpVar) else character()
  predict <- function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
    d <- curve_frame(distance_m, depth_m, tag_model, list(...))
    nd <- if (is.null(covariates)) d else covariates(d)
    for (g in setdiff(groups, names(nd))) nd[[g]] <- NA
    eta <- do.call(stats::predict, c(list(fit, newdata = nd, type = "link"), dots))
    if (is.list(eta)) eta <- eta$fit
    check_eta(eta, length(distance_m))
  }
  # With a covariates function the variables the model reads are unknown here,
  # so only distance is declared and missing inputs surface as NA predictions.
  vars <- all.vars(stats::formula(fit)[[3]])
  uses <- if (is.null(covariates)) intersect(curve_args, vars) else "distance_m"
  list(predict = predict, uses = uses, tag_models = fit_tag_models(fit),
       source = sprintf("%s, %s", class(fit)[1], deparse1(stats::formula(fit))))
}

# The checks on `covariates` and the arguments passed on to predict(), and
# the defaults model_curve() sets among them.
model_curve_dots <- function(fit, covariates, dots) {
  if (!is.null(covariates) && !is.function(covariates)) {
    stop("`covariates` must be a function or NULL.", call. = FALSE)
  }
  if ("type" %in% names(dots)) stop("`type` is fixed at \"link\".", call. = FALSE)
  if (inherits(fit, "glmmTMB") && !"re.form" %in% names(dots)) {
    warn_dropped_smooth(fit)
    dots$re.form <- NA
  }
  # predict.bam() on a discrete fit discretises the new data together, so a
  # pair's prediction moves with the other pairs predicted beside it.
  # discrete = FALSE predicts each row on its own; predict.bam()
  # ignores it for a fit that was not discrete.
  if (inherits(fit, "bam") && !"discrete" %in% names(dots)) {
    dots$discrete <- FALSE
  }
  dots
}

# glmmTMB fits the wiggly part of an s() smooth as a random effect, so
# re.form = NA predicts the smooth as its straight part with no message.
# The smooth is read from the model's formulas, so the check runs whatever
# glmmTMB version, or none, is installed.
warn_dropped_smooth <- function(fit) {
  smooths <- smooth_terms(fit)
  if (!length(smooths)) return(invisible(character()))
  warning(wrap_text(sprintf(paste0(
    "This glmmTMB model has the smooth term%s %s. glmmTMB fits the wiggly part of a smooth ",
    "as a random effect, and de_model() predicts with re.form = NA, so the smooth reaches ",
    "the field as a straight line. Fit the smooth in mgcv (gam() or bam()), or pass ",
    "re.form = NULL to de_model() to keep it; re.form = NULL also adds the fitted random ",
    "effect of any grouping column the curve is given."
  ), if (length(smooths) > 1) "s" else "", paste(smooths, collapse = ", "))), call. = FALSE)
  invisible(smooths)
}

# The smooth terms (s(), te(), t2(), ti(), bare or with mgcv::) in the
# conditional, zero-inflation and dispersion formulas of a glmmTMB fit.
smooth_terms <- function(fit) {
  forms <- tryCatch(fit$modelInfo$allForm[c("formula", "ziformula", "dispformula")],
                    error = function(e) NULL)
  forms <- Filter(function(f) inherits(f, "formula"), forms)
  if (!length(forms)) {
    forms <- list(tryCatch(stats::formula(fit), error = function(e) NULL))
  }
  found <- character()
  for (f in forms) if (inherits(f, "formula")) found <- c(found, smooth_calls(f[[length(f)]]))
  unique(found)
}

# The smooth calls in one expression, depth first, in the order they are written.
smooth_calls <- function(e) {
  if (!is.call(e)) return(character())
  head <- e[[1]]
  if (is.call(head) && identical(head[[1]], as.name("::")) && length(head) == 3) head <- head[[3]]
  if (is.name(head) && as.character(head) %in% c("s", "te", "t2", "ti")) return(deparse1(e))
  found <- character()
  for (a in as.list(e)[-1]) found <- c(found, smooth_calls(a))
  found
}

# The tag models a fitted model knows, when it has a tag_model term.
fit_tag_models <- function(fit) {
  if (!inherits(fit, "glmmTMB")) return(fit$xlevels$tag_model)
  v <- fit$frame$tag_model
  if (is.factor(v)) levels(v) else sort(unique(v))
}

function_curve <- function(f) {
  fmls <- names(formals(f))
  if (!any(c("distance_m", "...") %in% fmls)) {
    stop("The function must take an argument `distance_m`.", call. = FALSE)
  }
  open <- "..." %in% fmls
  predict <- function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
    args <- c(list(distance_m = distance_m, depth_m = depth_m, tag_model = tag_model), list(...))
    if (!open) args <- args[names(args) %in% fmls]
    check_eta(do.call(f, args), length(distance_m))
  }
  list(predict = predict, uses = if (open) "distance_m" else intersect(curve_args, fmls),
       source = "function")
}
