# nolint start: object_usage_linter.
# The fixtures come from helper-fixtures.R, which the linter cannot see.

surface_toy <- function(...) {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "a", 0),
                  counts_for(field, c(0, 0, 2, 5), 12, "a", 1),
                  counts_for(field, c(0, 4, 0, 0), 12, "b", 0))
  location_surface(counts, field, ...)
}

test_that("the region is the smallest set of cells holding the level", {
  s <- surface_toy()
  r <- credible_region(s, 0.9)
  expect_s3_class(r, "credible_region")
  summ <- attr(r, "summary")
  expect_equal(nrow(summ), 3)
  for (i in seq_len(nrow(summ))) {
    one <- s[s$animal_id == summ$animal_id[i] & s$bin_start == summ$bin_start[i], ]
    p <- sort(one$prob, decreasing = TRUE)
    expect_equal(summ$n_cells[i], which(cumsum(p) >= 0.9)[1])
    expect_gte(summ$mass[i], 0.9)
  }
  expect_equal(summ$area_m2, summ$n_cells * 900)
  expect_output(print(r), "level 0.9: 3 animal-bins")
})

test_that("a region matches brute-force ranking on the square lake's inside case", {
  field <- fixture_square_lake()
  s <- location_surface(counts_for(field, c(9, 7, 3, 12, 3), 20), field, "none")
  expect_equal(attr(credible_region(s), "summary")$n_cells, 183)
  expect_equal(attr(credible_region(s), "summary")$area_m2, 183 * 900)
})

test_that("the level cannot exceed the kept mass", {
  s <- surface_toy(control = list(kept_mass = 0.9))
  expect_error(credible_region(s, 0.95), "at most the kept mass of the surface, 0.9")
  expect_error(credible_region(s, 0), "above 0")
  expect_error(credible_region(data.frame(a = 1)), "must be a location surface")
})

test_that("a surface cut below the level warns that its regions fall short", {
  # After subsetting, the cells left can hold less than the level (0.507
  # against 0.95 here); the region says so.
  s <- surface_toy()
  top <- s[s$prob > 0.02, ]
  expect_warning(r <- credible_region(top, 0.95), "3 of 3 bins store less than the level 0.95")
  expect_equal(attr(r, "summary")$n_cells, attr(top, "bins")$n_cells)
  expect_silent(credible_region(s, 0.95))
})

test_that("the raster holds one layer per animal-bin with NA outside the region", {
  s <- surface_toy()
  r <- credible_region(s, 0.5, raster = TRUE)
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::nlyr(r), 3)
  tab <- credible_region(s, 0.5)
  v <- terra::values(r[[1]], mat = FALSE)
  first <- tab[tab$animal_id == "a" & tab$bin_start == min(tab$bin_start), ]
  expect_equal(sort(which(!is.na(v))), sort(first$cell))
  expect_equal(sum(v, na.rm = TRUE), sum(first$prob))
  expect_match(names(r)[1], "^a 2024-06-01 00:00:00$")
})

test_that("areas use each cell's own size on a longitude-latitude grid", {
  g <- terra::rast(nrows = 10, ncols = 10, xmin = -78, xmax = -77.9, ymin = 44.5, ymax = 44.6,
                   crs = "EPSG:4326", vals = 1)
  expect_equal(cr_cell_area(g, 1:3), terra::cellSize(g, unit = "m")[1:3][[1]])
  expect_true(all(is.na(cr_cell_area(NULL, 1:3))))
  s <- surface_toy()
  attr(s, "grid") <- NULL
  expect_true(all(is.na(attr(credible_region(s), "summary")$area_m2)))
  expect_error(credible_region(s, raster = TRUE), "no grid")
})

# nolint end
