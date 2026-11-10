# ---- Helpers ----
# The geom of each layer, by class name, so a test reads what a plot draws.
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "", USE.NAMES = FALSE)

# The column an aesthetic maps, as text.
mapped <- function(q) sub("^~", "", deparse(q))

# Every plot must also draw: building the gtable runs each layer's draw code.
expect_draws <- function(p) {
  withr::local_pdf(NULL)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
}

# ---- plot.de_model ----
test_that("a detection model plots one line per depth on the positionR theme", {
  m <- de_logistic(d50 = c(200, 400), d05 = c(600, 1000), depth = c(2, 20))
  p <- plot(m, depth = c(2, 10, 20))
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), "GeomLine")
  expect_identical(levels(p$data$curve), c("2 m depth", "10 m depth", "20 m depth"))
  expect_identical(mapped(p$mapping$x), "distance")
  expect_identical(mapped(p$mapping$y), "p")
  expect_identical(mapped(p$mapping$colour), "curve")
  expect_identical(p$labels$x, "Distance through water (m)")
  expect_identical(p$labels$y, "Detection probability per transmission")
  expect_true(all(p$data$p >= 0 & p$data$p <= 1))
  expect_identical(p$scales$get_scales("colour")$guide, "legend")
  expect_equal(p$scales$get_scales("y")$limits, c(0, 1))
  expect_equal(p$theme, theme_positionr())
  expect_draws(p)
  one <- plot(de_logistic(300, 800), distance = seq(0, 1000, by = 50))
  expect_equal(range(one$data$distance), c(0, 1000))
  expect_identical(one$scales$get_scales("colour")$guide, "none")
  expect_draws(one)
  expect_draws(plot(de_model(function(distance_m) rep(2, length(distance_m)))))
})

test_that("a tag-model curve plots every model it knows by default", {
  # A tag-model table with no tag_model chosen draws every model.
  tab <- data.frame(distance_m = c(0, 1000, 0, 1000), eta = c(1, -3, 1, -5),
                    tag_model = c("V13", "V13", "V7", "V7"))
  p <- plot(de_model(tab))
  expect_identical(levels(p$data$curve), c("V13", "V7"))
  expect_identical(levels(plot(de_model(tab), tag_model = c("V7", "V13"))$data$curve),
                   c("V7", "V13"))
  expect_draws(p)
  by_tag <- de_model(function(distance_m, tag_model) -distance_m / 300)
  expect_error(plot(by_tag), "name the tag models to draw")
  expect_identical(levels(plot(by_tag, tag_model = "V9")$data$curve), "V9")
})

test_that("arguments the plot methods do not use are ignored with a warning", {
  m <- de_logistic(300, 800)
  expect_warning(p <- plot(m, NULL, NA, NULL, "blue", main = "x"),
                 "an unnamed argument, `main`")
  expect_s3_class(p, "ggplot")
  expect_silent(plot(m))
})

# ---- plot.de_field ----
test_that("a field maps the whole array from 0 to 1 with receivers and a scale bar", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(300, 800))
  p <- plot(f)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), c("GeomRaster", "GeomPoint", "GeomPoint", "GeomPoint",
                               "GeomPoint", "ScaleBarRect", "ScaleBarText"))
  expect_identical(p$labels$title, "Detection probability, any receiver")
  expect_identical(c(p$labels$x, p$labels$y),
                   c("Easting (m, UTM 17N)", "Northing (m, UTM 17N)"))
  # One raster cell per water cell; land is not drawn.
  expect_equal(nrow(p$data), length(f$cells))
  fill <- p$scales$get_scales("fill")
  expect_equal(fill$limits, c(0, 1))
  expect_identical(fill$name, "Probability")
  expect_equal(p$coordinates$ratio, 1)
  # Receivers are open circles: shape 21 with no fill, yellow over a wider
  # dark ring, as the detecting receivers of plot.loc_surface().
  pts <- p$layers[[3]]
  expect_identical(pts$aes_params$shape, 21)
  expect_true(is.na(pts$aes_params$fill))
  expect_identical(pts$aes_params$colour, "#FDE725")
  expect_identical(p$layers[[2]]$aes_params$colour, "grey10")
  expect_gt(p$layers[[2]]$aes_params$stroke, pts$aes_params$stroke)
  marker_colours <- vapply(p$layers[2:5], function(l) l$aes_params$colour, "")
  expect_false(any(c("green", "#00FF00") %in% marker_colours))
  expect_equal(nrow(pts$data) + nrow(p$layers[[5]]$data), 2)
  expect_draws(p)
})

test_that("a single receiver is mapped on its own scale, its marker heavier", {
  # A receiver peaking near 0.5 is drawn on its own range of probability.
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(300, 800))
  p <- plot(f, receiver = "east")
  expect_identical(p$labels$title, "Detection probability, station east")
  expect_equal(p$scales$get_scales("fill")$limits, c(0, max(f$p[, "east"])))
  expect_lt(max(f$p[, "east"]), 1)
  expect_identical(p$layers[[5]]$data$station_id, "east")
  expect_gt(p$layers[[5]]$aes_params$stroke, p$layers[[3]]$aes_params$stroke)
  expect_draws(p)
  expect_error(plot(f, receiver = "north"), "one of the field's receivers")
  silent <- de_model(function(distance_m) rep(-Inf, length(distance_m)))
  g <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), silent)
  expect_equal(plot(g, receiver = "west")$scales$get_scales("fill")$limits, c(0, 1))
  serial <- data.frame(receiver_serial = "R1", x = 285, y = 435)
  h <- de_field(water_distances(peninsula_lake(), serial), de_logistic(250, 600))
  expect_identical(plot(h, receiver = "R1")$labels$title, "Detection probability, receiver R1")
})

test_that("a field map spans its water and receivers, not the whole grid", {
  # A grid with 2 km of land on its west side.
  grid <- terra::rast(nrows = 20, ncols = 100, xmin = 0, xmax = 3000, ymin = 0, ymax = 600,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowColCombine(grid, 1:20, 1:67)] <- NA
  receivers <- data.frame(station_id = c("a", "b"), x = c(2200, 2800), y = c(300, 300))
  f <- de_field(water_distances(grid, receivers), de_logistic(250, 600))
  lim <- plot(f)$coordinates$limits
  expect_equal(unname(lim$x), c(2010, 3000))
  # A lake of water to the grid's edge keeps the grid's extent.
  full <- plot(de_field(water_distances(peninsula_lake(), peninsula_receivers()),
                        de_logistic(300, 800)))
  expect_equal(unname(unlist(full$coordinates$limits)), c(0, 1200, 0, 600))
  # A receiver outside the water widens the map to take it in.
  e <- field_extent(data.frame(x = c(15, 45), y = c(15, 15)), data.frame(x = -100, y = 15),
                    c(30, 30))
  expect_equal(unname(e[c("xmin", "xmax")]), c(-130, 60))
})

test_that("the map's scale bar is a round length that fits the lake", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(300, 800))
  built <- ggplot2::ggplot_build(plot(f))
  bar <- scale_bar_parts(built$layout$panel_params[[1]], "bl", 0.15, FALSE)
  expect_identical(bar$text$label, "100 m")
})
