# nolint start: object_usage_linter.
# The fixtures come from helper-fixtures.R, which the linter cannot see.

surface_two_bins <- function() {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "a", 0),
                  counts_for(field, c(0, 0, 2, 5), 12, "a", 1),
                  counts_for(field, c(0, 4, 0, 0), 12, "b", 30))
  location_surface(counts, field)
}

test_that("draws are with replacement and follow the probabilities", {
  # Draws without replacement would spread to improbable cells once the
  # probable ones are taken.
  s <- surface_two_bins()
  one <- s[s$animal_id == "b", ]
  d <- sample_locations(one, 1e5, seed = 3)
  expect_equal(nrow(d), 1e5)
  expect_gt(anyDuplicated(d$cell), 0)
  freq <- table(factor(d$cell, levels = one$cell)) / 1e5
  expect_lt(max(abs(as.numeric(freq) - one$prob / sum(one$prob))), 0.005)
  expect_true(all(d$cell %in% one$cell))
})

test_that("each animal-bin gets n draws with cell centres", {
  s <- surface_two_bins()
  d <- sample_locations(s, 50, seed = 1)
  expect_named(d, c("animal_id", "bin_start", "draw", "cell", "x", "y"))
  expect_equal(as.vector(table(paste(d$animal_id, d$bin_start))), rep(50L, 3))
  xy <- terra::xyFromCell(sl_grid(s), d$cell)
  expect_equal(cbind(d$x, d$y), unname(xy))
  expect_equal(attr(d, "epsg"), 32617L)
  expect_match(attr(d, "estimand"), "each detected bin")
})

test_that("a seed gives the same draws and leaves the caller's stream alone", {
  s <- surface_two_bins()
  set.seed(99)
  before <- .Random.seed
  a <- sample_locations(s, 10, seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(sample_locations(s, 10, seed = 5), a)
  rm(".Random.seed", envir = globalenv())
  sample_locations(s, 1, seed = 5)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(sample_locations(s, 1), "data.frame")
})

test_that("bad arguments stop", {
  s <- surface_two_bins()
  expect_error(sample_locations(s, 0), "`n` must be")
  expect_error(sample_locations(s, 1.5), "`n` must be")
  expect_error(sample_locations(s[0, ], 1), "holds no cells")
  attr(s, "grid") <- NULL
  expect_true(all(is.na(sample_locations(s, 2)$x)))
  expect_true(is.na(attr(sample_locations(s, 2), "epsg")))
})

# nolint end
