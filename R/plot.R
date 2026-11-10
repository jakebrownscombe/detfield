#' Plot a detection model or a detection field
#'
#' `plot.de_model()` draws detection probability against distance, one line
#' per tag model and receiver depth asked for. `plot.de_field()` maps
#' one receiver's detection probability over the water cells, or the
#' probability that at least one receiver detects a transmission from each
#' cell when every receiver listens. The whole array is drawn on the fixed
#' scale 0 to 1. A single receiver is drawn from 0 to its own highest
#' probability, so that a receiver peaking near 0.5 is readable. Both use the
#' magma palette.
#'
#' Both methods build the figure with ggplot2 and [theme_positionr()]. The
#' map draws the probability with [scale_fill_probability()], the magma
#' palette of every positionR probability map, the receivers as open circles in the
#' viridis yellow `"#FDE725"` over a thin dark edge, the glyph of a detecting
#' receiver in [plot.loc_surface()], so they read on the dark and the pale
#' ends of magma (the mapped receiver larger), equal coordinates and a scale
#' bar at the bottom right. Land is left blank. The map spans the field's
#' water cells and its receivers, so a grid cut wider than the water draws
#' no empty margin. Axis titles name the grid's coordinate reference system,
#' as in [plot.loc_surface()].
#'
#' @param x A `de_model` or a `de_field`.
#' @param distance For `plot.de_model()`, the distances in metres to draw. By
#'   default 0 m to the distance where the last of the curves falls below
#'   0.01, at most 20 km, so every curve is drawn past its d05.
#' @param depth,tag_model For `plot.de_model()`, receiver depths and
#'   tag models; one line per combination. By default `tag_model`
#'   is every tag model the detection model knows.
#' @param receiver For `plot.de_field()`, the key of the receiver to map, or
#'   `NULL` to map the coverage of the whole array.
#' @param ... Not used. Change the figure by adding ggplot2 layers, scales or
#'   [ggplot2::labs()] to the plot it returns.
#'
#' @return A ggplot object. It prints when not assigned; assign it to add
#'   layers or to save it with [ggplot2::ggsave()].
#' @export
#' @examples
#' # The curve of the Stoney example at the shallowest, a middle and the
#' # deepest receiver.
#' m <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' plot(m, depth = c(1, 13, 25))
plot.de_model <- function(x, distance = NULL, depth = NA, tag_model = NULL, ...) {
  unused_dots(...)
  if (is.null(tag_model)) tag_model <- if (is.null(x$tag_models)) NA else x$tag_models
  if ("tag_model" %in% x$uses && anyNA(tag_model)) {
    stop("This detection model reads `tag_model`; name the tag models to draw.",
         call. = FALSE)
  }
  combos <- expand.grid(depth = depth, tag_model = tag_model, stringsAsFactors = FALSE)
  curve_at <- function(d, i) {
    stats::plogis(x$predict(d, depth_m = combos$depth[i], tag_model = combos$tag_model[i]))
  }
  if (is.null(distance)) {
    # The longest curve sets the range, so no curve is cut while still high.
    probe <- seq(0, 20000, by = 10)
    reach <- vapply(seq_len(nrow(combos)), function(i) {
      below <- which(curve_at(probe, i) < 0.01)
      if (length(below)) probe[below[1]] else 20000
    }, numeric(1))
    distance <- seq(0, max(reach), length.out = 201)
  }
  labels <- trimws(paste(ifelse(is.na(combos$tag_model), "", combos$tag_model),
                         ifelse(is.na(combos$depth), "", paste0(combos$depth, " m depth"))))
  curves <- data.frame(
    distance = rep(distance, nrow(combos)),
    p = unlist(lapply(seq_len(nrow(combos)), function(i) curve_at(distance, i))),
    curve = factor(rep(labels, each = length(distance)), levels = unique(labels))
  )
  keyed <- nrow(combos) > 1 || any(nzchar(labels))
  ggplot2::ggplot(curves, ggplot2::aes(distance, p, colour = curve)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_viridis_d(name = NULL, end = 0.85,
                                    guide = if (keyed) "legend" else "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Distance through water (m)",
                  y = "Detection probability per transmission") +
    theme_positionr()
}

#' @rdname plot.de_model
#' @export
#' @examples
#'
#' # The Stoney example: Stoney Lake (100 m cells, land NA) and 70 receivers.
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' # Beyond 2 km the Stoney curve is below 0.001 at every receiver; leaving
#' # those pairs out keeps the example quick.
#' field <- de_field(water_distances(grid, stoney_track$receivers,
#'                                   max_distance = 2000), curve)
#' plot(field, receiver = "S005")
#' plot(field)
plot.de_field <- function(x, receiver = NULL, ...) {
  unused_dots(...)
  if (is.null(receiver)) {
    values <- at_least_k(x$p, 1)
    title <- "Detection probability, any receiver"
    top <- 1
  } else {
    if (length(receiver) != 1 || !receiver %in% x$receivers) {
      stop(sprintf("`receiver` must be one of the field's receivers (%s).", x$key), call. = FALSE)
    }
    values <- x$p[, receiver]
    title <- sprintf("Detection probability, %s %s",
                     if (x$key == "station_id") "station" else "receiver", receiver)
    top <- if (max(values) > 0) max(values) else 1
  }
  grid <- as_water_grid(x$grid)
  water <- data.frame(terra::xyFromCell(grid, x$cells), p = values)
  rx <- x$receiver_xy
  hit <- rx[[x$key]] %in% receiver
  e <- field_extent(water, rx, terra::res(grid))
  titles <- axis_titles(grid)
  # Open yellow circles over a dark ring, as the detecting receivers of
  # plot.loc_surface(): the dark edge reads on pale magma, the yellow on
  # black, and the surface shows through. The mapped receiver is drawn larger.
  yellow <- positionr_ring_colours()[["heard"]]
  ggplot2::ggplot(water, ggplot2::aes(x, y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = p)) +
    ggplot2::geom_point(data = rx[!hit, ], shape = 21, colour = "grey10", fill = NA,
                        size = 2, stroke = 1.3) +
    ggplot2::geom_point(data = rx[!hit, ], shape = 21, colour = yellow, fill = NA,
                        size = 2, stroke = 0.6) +
    ggplot2::geom_point(data = rx[hit, ], shape = 21, colour = "grey10", fill = NA,
                        size = 3.5, stroke = 2.4) +
    ggplot2::geom_point(data = rx[hit, ], shape = 21, colour = yellow, fill = NA,
                        size = 3.5, stroke = 1.4) +
    scale_fill_probability(limits = c(0, top)) +
    ggplot2::coord_equal(xlim = e[c("xmin", "xmax")], ylim = e[c("ymin", "ymax")]) +
    ggplot2::labs(title = title, x = titles[["x"]], y = titles[["y"]]) +
    annotation_scale_bar(location = "br") +
    theme_positionr()
}

# Extent of a field map: the water cells to their edges, and each receiver
# with one cell around it. The whole grid's extent would draw kilometres of
# empty margin where a grid holds land beyond the outermost station.
field_extent <- function(water, rx, res) {
  xs <- c(water$x - res[1] / 2, water$x + res[1] / 2, rx$x - res[1], rx$x + res[1])
  ys <- c(water$y - res[2] / 2, water$y + res[2] / 2, rx$y - res[2], rx$y + res[2])
  c(xmin = min(xs, na.rm = TRUE), xmax = max(xs, na.rm = TRUE),
    ymin = min(ys, na.rm = TRUE), ymax = max(ys, na.rm = TRUE))
}

# The plot methods take `...` because plot() does. A ggplot object has no use
# for an argument given there, so say that it is ignored.
unused_dots <- function(...) {
  n <- ...length()
  if (n) {
    given <- names(list(...))
    given <- if (is.null(given)) rep("", n) else given
    warning(sprintf(paste0("Ignored %s; the plot is a ggplot object, so add layers, scales or ",
                           "labs() to it instead."),
                    paste(ifelse(nzchar(given), paste0("`", given, "`"), "an unnamed argument"),
                          collapse = ", ")),
            call. = FALSE)
  }
  invisible(NULL)
}
