test_that("each pair carries the model's logit, and missing pairs are -Inf", {
  grid <- peninsula_lake()
  d <- water_distances(grid, peninsula_receivers(), max_distance = 500)
  m <- de_logistic(250, 600)
  f <- de_field(d, m)
  expect_s3_class(f, "de_field")
  expect_identical(dim(f$eta), c(length(attr(d, "cells")), 2L))
  expect_identical(colnames(f$eta), c("west", "east"))
  expect_identical(f$receivers, c("west", "east"))
  expect_identical(f$key, "station_id")
  expect_equal(f$receiver_xy, data.frame(station_id = c("west", "east"), x = c(285, 885),
                                         y = c(435, 435)))
  rows <- match(d$cell, f$cells)
  cols <- match(d$station_id, f$receivers)
  expect_equal(f$eta[cbind(rows, cols)], m$predict(d$distance_m))
  expect_equal(sum(is.infinite(f$eta)), 2 * length(f$cells) - nrow(d))
  expect_true(all(f$eta[is.infinite(f$eta)] < 0))
  expect_equal(f$p, stats::plogis(f$eta))
  expect_true(all(f$p[is.infinite(f$eta)] == 0))
  expect_null(f$not_there)
})

test_that("as.data.frame() gives one row per pair", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  f <- de_field(d, de_logistic(250, 600))
  df <- as.data.frame(f)
  expect_identical(names(df), c("cell", "station_id", "distance_m", "in_water", "eta", "p"))
  expect_identical(df$in_water, d$in_water)
  expect_equal(nrow(df), nrow(d))
  expect_equal(df$p, stats::plogis(de_logistic(250, 600)$predict(d$distance_m)))
})

test_that("depth comes from deploy_z, one number or a named vector", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  m <- de_logistic(d50 = c(200, 400), d05 = c(600, 1000), depth = c(2, 20))
  f <- de_field(d, m)
  expect_equal(f$inputs$depth, c(west = 4, east = 6))
  west <- d$station_id == "west"
  expect_equal(f$eta[cbind(match(d$cell[west], f$cells), 1)],
               m$predict(d$distance_m[west], depth_m = 4))
  expect_equal(de_field(d, m, depth = 10)$inputs$depth, c(west = 10, east = 10))
  expect_equal(de_field(d, m, depth = c(east = 3, west = 9))$inputs$depth, c(west = 9, east = 3))
  expect_error(de_field(d, m, depth = c(1, 2)), "named by receiver key")
  expect_error(de_field(d, m, depth = c(west = 1)), "every receiver")
  expect_error(de_field(d, m, depth = "deep"), "numeric")
  no_z <- water_distances(peninsula_lake(), peninsula_receivers()[, 1:3])
  expect_error(de_field(no_z, m), "2 receivers have no depth")
  expect_true(all(is.na(de_field(no_z, de_logistic(250, 600))$inputs$depth)))
})

test_that("a model that reads tag_model needs one", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  tab <- data.frame(distance_m = c(0, 1000, 0, 1000), eta = c(1, -3, 1, -5),
                    tag_model = c("V13", "V13", "V7", "V7"))
  m <- de_model(tab)
  expect_error(de_field(d, m), "reads `tag_model`")
  expect_error(de_field(d, m, tag_model = c("V13", "V7")), "one tag model")
  f13 <- de_field(d, m, tag_model = "V13")
  f7 <- de_field(d, m, tag_model = "V7")
  expect_true(all(f7$eta <= f13$eta))
  expect_identical(f7$inputs$tag_model, "V7")
})

test_that("modulators are evaluated at the date and passed by name", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  curve <- function(distance_m, ...) {
    a <- list(...)
    -distance_m / 200 + a$ice + a$weed
  }
  m <- de_model(curve, modulators = list(
    ice = function(date) if (format(date, "%m") %in% c("01", "02")) -1 else 0,
    weed = c(east = 0.5, west = -0.5)
  ))
  winter <- de_field(d, m, date = as.Date("2024-01-15"))
  summer <- de_field(d, m, date = as.POSIXct("2024-07-01 12:00", tz = "UTC"))
  expect_equal(winter$eta - summer$eta, matrix(-1, nrow(winter$eta), 2,
                                               dimnames = dimnames(winter$eta)))
  expect_equal(summer$inputs$modulators, list(ice = c(0, 0), weed = c(-0.5, 0.5)))
  expect_error(de_field(d, m), "function of date; give `date`")
  expect_error(de_field(d, m, date = as.POSIXct("2024-07-01", tz = "America/Toronto")), "UTC")
  expect_error(de_field(d, m, date = "2024-07-01"), "Date or a POSIXct")
  expect_error(de_field(d, m, date = as.Date(c("2024-07-01", "2024-07-02"))), "one date")
  short <- de_model(curve, modulators = list(ice = 0, weed = c(east = 1, north = 1)))
  expect_error(de_field(d, short), "no value for receiver west")
  wrong <- de_model(curve, modulators = list(ice = c(1, 2), weed = 0))
  expect_error(de_field(d, wrong), "one value or one named value")
})

test_that("bad inputs stop", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  expect_error(de_field(as.data.frame(unclass(d)), de_logistic(250, 600)), "water_distances")
  expect_error(de_field(d, function(x) x), "must be a de_model")
})

test_that("printing summarises the field", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  f <- de_field(d, de_logistic(250, 600), date = as.Date("2024-06-01"))
  expect_output(print(f), "786 water cells x 2 receivers \\(by station_id\\)")
  expect_output(print(f), "date: 2024-06-01")
  tab <- data.frame(distance_m = c(0, 1000), eta = c(1, -3), tag_model = "V13")
  expect_output(print(de_field(d, de_model(tab), tag_model = "V13")), "tag_model: V13")
})

test_that("the blocked offset moves only pairs whose line of sight meets land", {
  # With no offset, the curve is read on blocked pairs as if sound travelled
  # round the land in open water.
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  base <- as.data.frame(de_field(d, de_logistic(250, 600)))
  shift <- as.data.frame(de_field(d, de_model(de_logistic(250, 600), blocked = -4.2)))
  blocked <- !d$in_water
  expect_equal(sum(blocked), 750)
  expect_equal(shift$eta[blocked] - base$eta[blocked], rep(-4.2, 750))
  expect_identical(shift$eta[!blocked], base$eta[!blocked])
  deaf <- de_field(d, de_model(de_logistic(250, 600), blocked = -Inf))
  expect_true(all(as.data.frame(deaf)$p[blocked] == 0))
  expect_identical(deaf$inputs$blocked, -Inf)
  # Printed lines wrap at the console width; width = 200 keeps each on a line.
  expect_output(print(deaf), "blocked line of sight: 750 of 1,572 pairs, receiver deaf behind land",
                width = 200)
  expect_output(print(de_field(d, de_logistic(250, 600))), "read at the distance through water",
                width = 200)
  old <- de_logistic(250, 600)
  old$blocked <- NULL
  expect_identical(de_field(d, old)$eta, de_field(d, de_logistic(250, 600))$eta)
})

test_that("field[[\"p\"]] gives the probability matrix", {
  # [[ and $ give the same matrix.
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  expect_identical(dim(f[["p"]]), c(786L, 2L))
  expect_identical(f[["p"]], f$p)
  expect_identical(f[["key"]], "station_id")
})

test_that("a field that stores p instead of eta returns the stored p", {
  # A field built with `p` and no `eta` returns `p`.
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  p <- f$p
  stored <- unclass(f)
  stored$eta <- NULL
  stored$p <- p
  class(stored) <- "de_field"
  expect_identical(stored$p, p)
  expect_identical(stored[["p"]], p)
  expect_null(stored$eta)
})

# ---- The blocked offset ----
test_that("blocked = 0 says how many heard pairs have land between them", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  base <- suppressMessages(de_field(d, de_logistic(250, 600)))
  heard <- base$eta[cbind(match(d$cell, base$cells), match(d$station_id, base$receivers))] >=
    stats::qlogis(0.05)
  n_shadowed <- sum(!d$in_water & heard)
  expect_gt(n_shadowed, 0)
  expect_message(de_field(d, de_logistic(250, 600)),
                 sprintf("%d of %d receiver-cell pairs with p >= 0.05 have land", n_shadowed,
                         sum(heard)))
  expect_silent(de_field(d, de_model(de_logistic(250, 600), blocked = -4.2)))
})
