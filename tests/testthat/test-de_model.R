range_test <- function() {
  data.frame(distance_m = rep(c(50, 200, 400, 600, 800), 2), depth_m = rep(c(4, 10), each = 5),
             heard = c(92, 71, 40, 14, 4, 95, 80, 52, 22, 7))
}

test_that("a glm is read on the link scale", {
  fit <- stats::glm(cbind(heard, 100 - heard) ~ distance_m + depth_m, family = stats::binomial,
                    data = range_test())
  m <- de_model(fit, sigma = list(station_hour = 1.36))
  nd <- data.frame(distance_m = c(0, 250, 900), depth_m = c(3, 6, 12))
  expect_equal(m$predict(nd$distance_m, depth_m = nd$depth_m),
               unname(stats::predict(fit, nd, type = "link")))
  expect_identical(m$uses, c("distance_m", "depth_m"))
  expect_equal(m$sigma, list(station_hour = 1.36, station_day = NA_real_))
  expect_match(m$source, "^glm, cbind")
  expect_error(de_model(fit, type = "response"), "fixed at")
  with_se <- de_model(fit, se.fit = TRUE)
  expect_equal(with_se$predict(250, depth_m = 6), m$predict(250, depth_m = 6))
  expect_error(m$predict("far", depth_m = 6), "numeric")
})

test_that("modulators reach a fitted model as columns", {
  d <- range_test()
  d$ice <- rep(c(0, 1), 5)
  fit <- stats::glm(cbind(heard, 100 - heard) ~ distance_m + ice, family = stats::binomial,
                    data = d)
  m <- de_model(fit, modulators = list(ice = 0))
  expect_equal(m$predict(100, ice = 1) - m$predict(100, ice = 0), unname(stats::coef(fit)["ice"]))
})

test_that("covariates map the standard columns onto the model's own", {
  d <- range_test()
  d$dist_km <- d$distance_m / 1000
  fit <- stats::glm(cbind(heard, 100 - heard) ~ dist_km, family = stats::binomial, data = d)
  m <- de_model(fit, covariates = function(x) data.frame(dist_km = x$distance_m / 1000))
  expect_equal(m$predict(c(100, 500)),
               unname(stats::predict(fit, data.frame(dist_km = c(0.1, 0.5)))))
  expect_identical(m$uses, "distance_m")
  expect_error(de_model(fit, covariates = "dist_km"), "function or NULL")
})

test_that("a glmmTMB fit is read with its random effects at zero", {
  # glmmTMB warns on load when its TMB build differs from the installed TMB;
  # that is a property of the library, not of this package.
  suppressWarnings(skip_if_not_installed("glmmTMB"))
  d <- range_test()
  d <- rbind(d, d)
  d$station <- factor(rep(c("a", "b"), each = 10))
  d$heard <- pmin(100, d$heard + rep(c(-3, 3), each = 10))
  fit <- suppressWarnings(glmmTMB::glmmTMB(cbind(heard, 100 - heard) ~ distance_m + (1 | station),
                                           family = stats::binomial, data = d))
  m <- de_model(fit)
  fixed <- glmmTMB::fixef(fit)$cond
  expect_equal(m$predict(c(0, 400)), unname(fixed[1] + fixed[2] * c(0, 400)), tolerance = 1e-6)
  expect_match(m$source, "^glmmTMB")
})

# A range test over six days whose logit wanders by day: a smooth over day
# has a wiggly part to drop.
smooth_test <- function() {
  d <- expand.grid(distance_m = c(50, 200, 400, 600, 800), day = 0:5, station = c("a", "b"))
  wiggle <- c(0.8, -0.6, 0.4, -0.9, 0.7, -0.2)[d$day + 1]
  d$heard <- round(100 * stats::plogis(2.5 - 0.006 * d$distance_m + wiggle +
                                         ifelse(d$station == "a", -0.1, 0.1)))
  d
}

# A glmmTMB fit with a smooth, or a skip when glmmTMB or its smooth support
# is missing (glmmTMB before 1.1.10 refuses s()).
fit_tmb_smooth <- function(d) {
  suppressWarnings(skip_if_not_installed("glmmTMB"))
  fit <- tryCatch(
    suppressWarnings(glmmTMB::glmmTMB(
      cbind(heard, 100 - heard) ~ distance_m + s(day, k = 4) + (1 | station),
      family = stats::binomial, data = d
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) skip("this glmmTMB does not fit s() smooths")
  fit
}

test_that("a glmmTMB smooth that re.form = NA would drop gives one warning", {
  d <- smooth_test()
  fit <- fit_tmb_smooth(d)
  w <- character()
  m <- withCallingHandlers(de_model(fit, modulators = list(day = 0)), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  expect_length(w, 1)
  w <- gsub("\\s+", " ", w)  # the message is wrapped to the console width
  expect_match(w, "smooth term s(day, k = 4)", fixed = TRUE)
  expect_match(w, "straight line")
  expect_match(w, "mgcv")
  expect_match(w, "re.form = NULL", fixed = TRUE)
  # The warning describes what happens: through re.form = NA the logit is
  # linear in day; with re.form = NULL it is not.
  by_day <- vapply(0:5, function(k) m$predict(300, day = k), numeric(1))
  expect_equal(diff(diff(by_day)), rep(0, 4), tolerance = 1e-8)
  kept <- expect_warning(de_model(fit, modulators = list(day = 0), re.form = NULL), NA)
  by_day_kept <- vapply(0:5, function(k) kept$predict(300, day = k), numeric(1))
  expect_gt(max(abs(diff(diff(by_day_kept)))), 0.1)
})

test_that("no smooth warning for glmmTMB without smooths, for mgcv, or with re.form given", {
  d <- smooth_test()
  suppressWarnings(skip_if_not_installed("glmmTMB"))
  plain <- suppressWarnings(glmmTMB::glmmTMB(cbind(heard, 100 - heard) ~ distance_m + day +
                                               (1 | station), family = stats::binomial, data = d))
  expect_warning(de_model(plain), NA)
  expect_length(smooth_terms(plain), 0)
  gam_fit <- mgcv::gam(cbind(heard, 100 - heard) ~ distance_m + s(day, k = 4),
                       family = stats::binomial, data = d)
  expect_warning(de_model(gam_fit, modulators = list(day = 0)), NA)
  fit <- fit_tmb_smooth(d)
  expect_warning(de_model(fit, modulators = list(day = 0), re.form = NA), NA)
  expect_warning(de_model(fit, modulators = list(day = 0), re.form = ~ (1 | station)), NA)
})

test_that("smooth terms are found in the formula alone, written bare or with mgcv::", {
  fake <- function(f) {
    list(modelInfo = list(allForm = list(formula = f, ziformula = ~0, dispformula = ~1)))
  }
  expect_identical(smooth_terms(fake(y ~ x + s(day) + te(x, z) + (1 | g))),
                   c("s(day)", "te(x, z)"))
  expect_identical(smooth_terms(fake(y ~ mgcv::t2(x, z) + log(s))), "mgcv::t2(x, z)")
  expect_identical(smooth_terms(fake(y ~ x + ti(x, z))), "ti(x, z)")
  expect_length(smooth_terms(fake(y ~ s + x)), 0)
  # A fit without allForm (older glmmTMB) falls back to formula().
  expect_length(smooth_terms(structure(list(), class = "not_a_model")), 0)
})

test_that("a table of curve values is interpolated on the logit scale", {
  tab <- data.frame(distance_m = c(0, 500, 1000), eta = c(1, -1, -3))
  m <- de_model(tab)
  expect_equal(m$predict(c(0, 250, 750, 1000)), c(1, 0, -2, -3))
  # Beyond the table the falling logit continues; before it, it is held.
  expect_equal(m$predict(1500), -5)
  expect_equal(m$predict(-10), 1)
  flat <- de_model(data.frame(distance_m = c(0, 500, 1000), eta = c(1, -1, -1)))
  expect_equal(flat$predict(2000), -1)
  pt <- de_model(data.frame(distance_m = c(0, 1000), p = c(0.5, stats::plogis(-2))))
  expect_equal(pt$predict(500), -1)
  expect_identical(m$uses, "distance_m")
})

test_that("tables by tag model and depth interpolate in depth and refuse unknown models", {
  tab <- expand.grid(distance_m = c(0, 1000), depth_m = c(2, 10), tag_model = c("V9", "V13"),
                     stringsAsFactors = FALSE)
  tab$eta <- 1 - tab$distance_m / 1000 * ifelse(tab$tag_model == "V9", 4, 2) - tab$depth_m / 10
  m <- de_model(tab)
  expect_identical(m$uses, c("distance_m", "depth_m", "tag_model"))
  truth <- function(d, z, tm) 1 - d / 1000 * ifelse(tm == "V9", 4, 2) - z / 10
  expect_equal(m$predict(c(250, 600), depth_m = c(6, 3), tag_model = c("V9", "V13")),
               truth(c(250, 600), c(6, 3), c("V9", "V13")))
  expect_equal(m$predict(500, depth_m = 30, tag_model = "V13"), truth(500, 10, "V13"))
  expect_error(m$predict(500, depth_m = 5, tag_model = "V16"), "Unknown tag_model V16")
  expect_error(m$predict(500, tag_model = "V13"), "needs `depth_m`")
  expect_error(m$predict("a", depth_m = 5, tag_model = "V13"), "numeric")
})

test_that("bad curve tables stop", {
  expect_error(de_model(data.frame(eta = 1:2)), "needs `distance_m`")
  expect_error(de_model(data.frame(distance_m = 1:2)), "`eta` or `p`")
  expect_error(de_model(data.frame(distance_m = 1:2, p = c(1, 0.5))), "strictly between")
  expect_error(de_model(data.frame(distance_m = 1, eta = 1)), "at least two rows")
  expect_error(de_model(data.frame(distance_m = c(1, 1), eta = 1:2)), "distinct distances")
  expect_error(de_model(data.frame(distance_m = 1:2, eta = c(1, NA))), "Column `eta`")
})

test_that("a function is called with the arguments it takes", {
  m <- de_model(function(distance_m) 2 - distance_m / 100)
  expect_equal(m$predict(c(0, 100), depth_m = 4, tag_model = "V7"), c(2, 1))
  expect_identical(m$uses, "distance_m")
  open <- de_model(function(distance_m, ...) {
    a <- list(...)
    -distance_m / 100 + a$depth_m + a$ice
  }, modulators = list(ice = -1))
  expect_equal(open$predict(100, depth_m = 3, ice = -1), 1)
  expect_identical(open$uses, "distance_m")
  expect_error(de_model(function(d) d), "argument `distance_m`")
  expect_error(de_model(function(distance_m) NA)$predict(1), "returned NA for 1 of 1")
  expect_error(de_model(function(distance_m) 1:3)$predict(1:2), "returned 3 values for 2 distances")
})

test_that("sigma and modulators are checked and carried, never fitted", {
  m <- de_logistic(300, 800)
  m2 <- de_model(m, sigma = list(station_hour = 1.362), modulators = list(season = "summer"))
  expect_equal(m2$sigma$station_hour, 1.362)
  expect_true(is.na(m2$sigma$station_day))
  expect_identical(m2$modulators, list(season = "summer"))
  expect_identical(m2$predict, m$predict)
  expect_error(de_model(m, sigma = list(station_week = 1)), "takes only")
  expect_error(de_model(m, sigma = list(1)), "takes only")
  expect_error(de_model(m, sigma = list(station_hour = -1)), "at least 0")
  expect_error(de_model(m, sigma = 1), "must be a list")
  expect_error(de_model(m, modulators = "ice"), "must be a list")
  expect_error(de_model(m, modulators = list(1)), "needs a name")
  expect_error(de_model(m, modulators = list(depth_m = 1)), "reserved")
  expect_error(de_model(1:3), "must be a fitted glm")
})

test_that("printing names the source, the arguments and sigma", {
  m <- de_model(de_logistic(300, 800), sigma = list(station_hour = 1.362),
                modulators = list(ice = 0))
  expect_output(print(m), "de_logistic, d50 300 m, d05 800 m")
  expect_output(print(m), "sigma: station_hour 1.362 logits")
  expect_output(print(de_logistic(300, 800)), "sigma: station_hour not given")
  expect_output(print(m), "modulators: ice")
})

test_that("blocked is checked and carried, and printed", {
  m <- de_model(de_logistic(300, 800), blocked = -4.2)
  expect_identical(m$blocked, -4.2)
  expect_identical(de_logistic(300, 800)$blocked, 0)
  expect_identical(de_model(m)$blocked, 0)
  expect_output(print(m), "blocked line of sight: curve plus -4.2 logits")
  expect_error(de_model(m, blocked = 1), "0 or less")
  expect_error(de_model(m, blocked = c(-1, -2)), "one number")
  expect_error(de_model(m, blocked = NA_real_), "one number")
})

test_that("tables and fitted models report the tag models they know", {
  tab <- data.frame(distance_m = c(0, 1000, 0, 1000), eta = c(1, -3, 1, -5),
                    tag_model = c("V13", "V13", "V7", "V7"))
  expect_identical(de_model(tab)$tag_models, c("V13", "V7"))
  d <- range_test()
  d$tag_model <- rep(c("V13", "V9"), 5)
  fit <- stats::glm(cbind(heard, 100 - heard) ~ distance_m + tag_model, family = stats::binomial,
                    data = d)
  expect_identical(de_model(fit)$tag_models, c("V13", "V9"))
  expect_null(de_logistic(300, 800)$tag_models)
  suppressWarnings(skip_if_not_installed("glmmTMB"))
  d$station <- factor(rep(c("a", "b"), each = 5))
  tmb <- suppressWarnings(glmmTMB::glmmTMB(cbind(heard, 100 - heard) ~ distance_m + tag_model +
                                             (1 | station), family = stats::binomial, data = d))
  expect_identical(de_model(tmb)$tag_models, c("V13", "V9"))
  d$tag_model <- factor(d$tag_model, levels = c("V9", "V13"))
  tmb_f <- suppressWarnings(glmmTMB::glmmTMB(cbind(heard, 100 - heard) ~ distance_m + tag_model +
                                               (1 | station), family = stats::binomial, data = d))
  expect_identical(de_model(tmb_f)$tag_models, c("V9", "V13"))
})

test_that("a curve recycles arguments of length one and stops on other mismatches", {
  l <- de_logistic(d50 = c(200, 400), d05 = c(600, 1000), depth = c(2, 20))
  expect_equal(l$predict(distance_m = 300, depth_m = c(2, 20)),
               l$predict(distance_m = c(300, 300), depth_m = c(2, 20)))
  expect_length(l$predict(distance_m = 300, depth_m = c(2, 20)), 2)
  tab <- data.frame(distance_m = rep(c(0, 500, 1000), 2), eta = c(2, 0, -2, 3, 1, -1),
                    depth_m = rep(c(5, 15), each = 3))
  m <- de_model(tab)
  expect_equal(m$predict(0, depth_m = c(5, 15)), c(2, 3))
  fit <- stats::glm(cbind(h, 10 - h) ~ distance_m + depth_m, binomial,
                    data.frame(distance_m = rep(c(0, 300, 600), 2),
                               depth_m = rep(c(5, 15), each = 3), h = c(9, 5, 1, 10, 7, 3)))
  g <- de_model(fit)
  expect_equal(g$predict(0, depth_m = c(5, 15)), g$predict(c(0, 0), depth_m = c(5, 15)))
  expect_error(l$predict(c(100, 200, 300), depth_m = c(2, 20)), "common length")
  expect_length(l$predict(numeric(), depth_m = 2), 0)
})
