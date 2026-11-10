# ---- Helpers ----
# Every layer of a list, as the colour, stroke and data rows it draws.
ring_table <- function(layers) {
  points <- Filter(function(l) inherits(l, "LayerInstance"), unlist(list(layers)))
  data.frame(colour = vapply(points, function(l) l$aes_params$colour, ""),
             stroke = vapply(points, function(l) l$aes_params$stroke, 0),
             rows = vapply(points, function(l) nrow(l$data), 0))
}

rx <- data.frame(station_id = c("A", "B", "C"), x = c(0, 100, 200), y = c(0, 50, 0),
                 detections = c(0, 4, 12))

# ---- layer_receivers ----
test_that("without counts every receiver is a yellow ring over a dark ring", {
  layers <- layer_receivers(rx)
  rings <- ring_table(layers)
  expect_identical(rings$colour, c("#1A1A1A", "#FDE725"))
  expect_equal(rings$stroke, c(2.2, 1.1))
  expect_equal(rings$rows, c(3, 3))
  expect_true(all(vapply(layers, function(l) is.na(l$aes_params$fill), TRUE)))
  expect_false(any(vapply(layers, function(l) l$inherit.aes, TRUE)))
})

test_that("with counts the heard receivers are yellow and sized, the silent white", {
  layers <- layer_receivers(rx, count = "detections", max_size = 5)
  rings <- ring_table(layers)
  expect_identical(rings$colour, c("#1A1A1A", "#FFFFFF", "#1A1A1A", "#FDE725"))
  expect_equal(rings$rows, c(1, 1, 2, 2))
  scales <- Filter(function(l) inherits(l, "Scale"), layers)
  expect_setequal(vapply(scales, function(s) s$aesthetics[1], ""), c("size", "shape"))
  p <- ggplot2::ggplot(rx, ggplot2::aes(x, y)) + layer_receivers(rx, count = "detections") +
    theme_positionr()
  withr::local_pdf(NULL)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("a count column with no silent receiver draws no white ring", {
  heard_only <- rx[rx$detections > 0, ]
  rings <- ring_table(layer_receivers(heard_only, count = "detections"))
  expect_false("#FFFFFF" %in% rings$colour)
})

test_that("bad input stops with a plain message", {
  expect_error(layer_receivers(rx[, c("station_id", "x")]), "columns `x` and `y`")
  expect_error(layer_receivers(rx, count = "heard"), "name of a column")
  expect_error(layer_receivers(rx[0, ]), "no rows")
  rx$detections[2] <- NA
  expect_error(layer_receivers(rx, count = "detections"), "counts of 0 or more")
})
