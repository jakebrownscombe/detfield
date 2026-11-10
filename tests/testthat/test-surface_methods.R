# nolint start: object_usage_linter.
# The fixtures come from helper-fixtures.R, which the linter cannot see.

surface_for_methods <- function() {
  field <- fixture_field(fixture_grid(10, land = 1:5), fixture_receivers(),
                         function(d) 0.5 - d / 80)
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "a", 0),
                  counts_for(field, c(0, 0, 2, 5), 12, "b", 1))
  location_surface(counts, field)
}

test_that("subsetting keeps the attributes and trims the bin table", {
  s <- surface_for_methods()
  b <- s[s$animal_id == "b", ]
  expect_s3_class(b, "loc_surface")
  expect_equal(attr(b, "bins")$animal_id, "b")
  expect_equal(unique(attr(b, "receivers")$animal_id), "b")
  expect_equal(attr(b, "overdispersion"), attr(s, "overdispersion"))
  expect_false(inherits(s[, c("cell", "prob")], "loc_surface"))
  # The bin table follows the rows left.
  top <- s[s$prob > 0.02, ]
  b <- attr(top, "bins")
  expect_equal(b$n_cells, as.vector(table(paste(top$animal_id, top$bin_start))))
  expect_equal(b$kept_mass, as.vector(tapply(top$prob, paste(top$animal_id, top$bin_start),
                                             sum)))
  expect_equal(s[2, "prob"], s$prob[2])
})

test_that("plot maps a bin and stops on a missing one", {
  s <- surface_for_methods()
  mapped <- function(q) sub("^~", "", deparse(q))
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "", USE.NAMES = FALSE)
  draws <- function(p) {
    withr::local_pdf(NULL)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
  p <- plot(s)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), c("GeomRaster", "GeomRaster", "GeomPath", "GeomPath", "GeomPoint",
                               "GeomPoint", "GeomPoint", "GeomPoint", "ScaleBarRect",
                               "ScaleBarText"))
  # The title is one line above the panel and the keys, so it cannot run
  # into them.
  expect_identical(p$labels$title, "a, 2024-06-01 00:00 UTC")
  expect_identical(p$theme$plot.title.position, "plot")
  # The axis titles name the grid's CRS.
  expect_identical(c(p$labels$x, p$labels$y),
                   c("Easting (m, UTM 17N)", "Northing (m, UTM 17N)"))
  # The stored cells carry their probability; other water cells are NA.
  bin <- s[s$animal_id == "a", ]
  expect_equal(sum(p$data$prob, na.rm = TRUE), sum(bin$prob))
  expect_equal(sum(!is.na(p$data$prob)), nrow(bin))
  expect_identical(mapped(p$layers[[2]]$mapping$fill), "prob")
  # The fill runs from 0 to the bin's largest cell probability, which is the
  # top label of the key.
  fill <- p$scales$get_scales("fill")
  expect_identical(fill$name, "Probability\nper cell")
  expect_equal(fill$limits, c(0, max(bin$prob)))
  expect_equal(max(fill$get_breaks()), max(bin$prob))
  # The receivers are drawn by layer_receivers(), the one receiver mark of
  # the suite: silent receivers are white rings over a dark ring (layers 5
  # and 6), keyed "No detections"; receivers that heard the animal are
  # yellow rings over a dark ring (layers 7 and 8), sized by count.
  expect_identical(mapped(p$layers[[8]]$mapping$size), ".count")
  expect_true(all(p$layers[[8]]$data$y_count > 0))
  expect_true(all(p$layers[[5]]$data$y_count == 0))
  expect_identical(p$layers[[5]]$aes_params$colour, "#1A1A1A")
  expect_identical(p$layers[[6]]$aes_params$colour, "#FFFFFF")
  expect_true(is.na(p$layers[[6]]$aes_params$fill))
  expect_gt(p$layers[[5]]$aes_params$stroke, p$layers[[6]]$aes_params$stroke)
  expect_identical(mapped(p$layers[[6]]$mapping$shape), ".key")
  expect_identical(unique(p$layers[[6]]$data$.key), "No detections")
  expect_identical(p$scales$get_scales("shape")$name, NULL)
  # Detecting receivers are yellow over a wider dark ring,
  # and the 50 and 95 percent regions are drawn in two keyed colours that
  # differ from each other and from the bubbles.
  expect_identical(p$layers[[8]]$aes_params$colour, "#FDE725")
  expect_identical(p$layers[[7]]$aes_params$colour, "#1A1A1A")
  expect_gt(p$layers[[7]]$aes_params$stroke, p$layers[[8]]$aes_params$stroke)
  expect_identical(p$layers[[3]]$aes_params$colour, "grey10")
  expect_identical(mapped(p$layers[[4]]$mapping$colour), "region")
  expect_identical(levels(p$layers[[4]]$data$region), c("50%", "95%"))
  region_scale <- p$scales$get_scales("colour")
  expect_identical(region_scale$name, "Credible region")
  region_cols <- region_scale$palette(2)
  expect_identical(unname(region_cols), c("#56B4E9", "#009E73"))
  expect_false(any(toupper(region_cols) %in% toupper(p$layers[[8]]$aes_params$colour)))
  # The 50 percent region lies inside the 95 percent region.
  r <- p$layers[[4]]$data
  expect_true(all(range(r$x[r$region == "50%"]) >= min(r$x[r$region == "95%"])))
  expect_true(all(range(r$x[r$region == "50%"]) <= max(r$x[r$region == "95%"])))
  one <- plot(s, level = 0.95)
  expect_identical(levels(one$layers[[4]]$data$region), "95%")
  expect_identical(unname(one$scales$get_scales("colour")$palette(1)), "#009E73")
  expect_error(plot(s, level = c(0.5, 0.8, 0.95)), "one or two numbers")
  expect_error(plot(s, level = 1.5), "one or two numbers")
  # The map takes in every receiver that detected the animal.
  rx <- p$layers[[8]]$data
  lim <- p$coordinates$limits
  expect_true(all(rx$x >= lim$x[1] & rx$x <= lim$x[2] & rx$y >= lim$y[1] & rx$y <= lim$y[2]))
  expect_equal(p$coordinates$ratio, 1)
  draws(p)
  b <- plot(s, animal_id = "b", bin_start = "2024-06-01 01:00:00", level = 0.5)
  expect_identical(b$labels$title, "b, 2024-06-01 01:00 UTC")
  draws(b)
  expect_false("GeomPath" %in% geoms(plot(s, level = NULL)))
  wide <- plot(s, ext = c(0, 300, 0, 300))
  expect_equal(unname(wide$coordinates$limits$x), c(0, 300))
  # Receivers outside the map are not drawn or keyed: the silent receivers
  # at y = 225 fall outside y 0 to 150, so no white ring and no
  # "No detections" key.
  south <- plot(s, ext = c(0, 300, 0, 150))
  expect_identical(sum(geoms(south) == "GeomPoint"), 2L)
  expect_true(all(south$layers[[5]]$data$y_count > 0))
  expect_null(south$scales$get_scales("shape"))
  draws(south)
  expect_warning(plot(s, xlab = "x"), "Ignored `xlab`")
  expect_error(plot(s, animal_id = "z"), "Animal z is not in the surface")
  expect_error(plot(s, animal_id = "a", bin_start = "2025-01-01"), "has no bin starting")
  attr(s, "receiver_xy") <- NULL
  expect_false("GeomPoint" %in% geoms(plot(s)))
  draws(plot(s))
  attr(s, "grid") <- NULL
  expect_error(plot(s), "no grid")
})

test_that("the size of a detecting bubble is proportional in area to its count", {
  s <- surface_for_methods()
  p <- plot(s)
  sz <- ggplot2::ggplot_build(p)$plot$scales$get_scales("size")
  # Area from zero: the largest count draws at 6 mm, a smaller one at 6 mm
  # times the square root of its share of the largest.
  # The scale starts at zero, so a key count below the bin's smallest count
  # still draws.
  lim <- sz$get_limits()
  expect_equal(lim[1], 0)
  expect_equal(sz$map(lim[2]), 6)
  expect_equal(sz$map(lim[2] / 4), 3, tolerance = 1e-6)
  expect_false(anyNA(sz$map(sz$get_breaks())))
  big <- ggplot2::scale_size_area(max_size = 6)
  big$train(c(2, 10))
  drawn <- big$map(c(2, 10))
  expect_equal(drawn[2] / drawn[1], sqrt(5), tolerance = 1e-6)
  # The key shows at most three whole counts, inside the limits.
  expect_identical(count_breaks(c(0, 10)), c(5, 10))
  expect_identical(count_breaks(c(0, 17)), c(5, 10, 15))
  expect_identical(count_breaks(c(0, 1)), 1)
  expect_identical(count_breaks(c(0, 3)), c(1, 2, 3))
  expect_lte(length(count_breaks(c(0, 140))), 3)
  expect_identical(count_breaks(c(0, 0)), numeric())
})

test_that("rings of a region small on the map are drawn at half width", {
  rings <- data.frame(x = c(0, 30, 30, 0), y = c(0, 0, 60, 60))
  # A 60 m region in a 3,000 m view: half width; in a 300 m view: full width.
  expect_identical(ring_width(rings, c(xmin = 0, xmax = 3000, ymin = 0, ymax = 2000)), 0.3)
  expect_identical(ring_width(rings, c(xmin = 0, xmax = 300, ymin = 0, ymax = 200)), 0.6)
  expect_identical(ring_width(rings, NULL), 0.6)
  s <- surface_for_methods()
  p <- plot(s)
  expect_equal(p$layers[[3]]$aes_params$linewidth, p$layers[[4]]$aes_params$linewidth + 0.5)
  wide <- plot(s, ext = c(-20000, 20000, -20000, 20000))
  expect_identical(wide$layers[[4]]$aes_params$linewidth, 0.3)
  expect_identical(wide$layers[[3]]$aes_params$linewidth, 0.8)
})

test_that("axis titles name the CRS when the grid carries one", {
  utm <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                     crs = "EPSG:32617")
  expect_identical(axis_titles(utm), c(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)"))
  south <- terra::rast(utm)
  terra::crs(south) <- "EPSG:32733"
  expect_identical(crs_label(south), "UTM 33S")
  mtm <- terra::rast(utm)
  terra::crs(mtm) <- "EPSG:3161"
  expect_identical(crs_label(mtm), "EPSG:3161")
  none <- terra::rast(utm)
  terra::crs(none) <- ""
  expect_null(crs_label(none))
  lonlat <- terra::rast(nrows = 2, ncols = 2, crs = "EPSG:4326")
  expect_identical(axis_titles(lonlat), c(x = "x", y = "y"))
})

test_that("the zoom leaves out silent receivers far from the bin", {
  grid <- terra::rast(nrows = 10, ncols = 100, xmin = 0, xmax = 1000, ymin = 0, ymax = 100,
                      crs = "EPSG:32617", vals = 1)
  rx <- data.frame(x = c(55, 105, 950), y = c(50, 50, 50), y_count = c(4, 0, 0))
  e <- as.vector(pl_zoom(grid, terra::cellFromXY(grid, cbind(c(45, 55), c(45, 55))), rx))
  # The near silent receiver widens the map; the one 900 m off does not.
  expect_gte(e[["xmax"]], 105)
  expect_lt(e[["xmax"]], 950)
})

test_that("the region outline closes around the region's cells", {
  grid <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
                      crs = "EPSG:32617", vals = 1)
  rings <- cell_outline(grid, c(6, 7))
  # The outline lies 1% of a cell outside the cells' edges.
  expect_equal(range(rings$x), c(1, 3) + c(-0.01, 0.01), tolerance = 1e-9)
  expect_equal(range(rings$y), c(2, 3) + c(-0.01, 0.01), tolerance = 1e-9)
  expect_length(unique(rings$ring), 1)
  expect_equal(nrow(cell_outline(grid, integer())), 0)
})

test_that("a grid stored as a file path is read back", {
  s <- surface_for_methods()
  path <- tempfile(fileext = ".tif")
  terra::writeRaster(sl_grid(s), path)
  attr(s, "grid") <- path
  expect_s4_class(sl_grid(s), "SpatRaster")
  attr(s, "grid") <- terra::rast(path)
  expect_s4_class(sl_grid(s), "SpatRaster")
  attr(s, "grid") <- 5
  expect_warning(expect_null(sl_grid(s)), "grid \\(numeric\\) cannot be read")
  expect_error(sl_grid(s, strict = TRUE), "cannot be read")
})

# nolint end

test_that("a region is outlined as one shape, cells touching at a corner included", {
  grid <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 400, ymin = 0, ymax = 400,
                      crs = "EPSG:32617")
  # Cells 6 and 11 touch only at a corner; cell 1 stands apart.
  joined <- cell_outline(grid, c(6, 11))
  expect_length(unique(joined$ring), 1)
  apart <- cell_outline(grid, c(1, 11))
  expect_length(unique(apart$ring), 2)
  # The outline moves by 1% of a cell at most.
  expect_lte(max(abs(range(joined$x) - c(100, 300))), 1 + 1e-9)
})
