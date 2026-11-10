#' Receivers as rings on a map
#'
#' `layer_receivers()` returns the ggplot2 layers that draw receivers in the
#' positionR style: a ring over a dark ring one step wider (`"#1A1A1A"`),
#' never filled, so the map shows through. A receiver that heard the animal
#' is a yellow ring (`"#FDE725"`); a listening receiver that heard none is a
#' white ring, which reads on the deep blue of a depth map and on the dark end
#' of magma.
#'
#' Without `count`, every receiver is a yellow ring of size `size`. With
#' `count`, the receivers with a count above 0 are yellow rings sized by the
#' count, area proportional to it from zero up to `max_size`, with a key
#' titled "Detections"; the receivers with a count of 0 are white rings of
#' size `size`, keyed "No detections".
#'
#' @param data A data frame with the receiver positions in map coordinates,
#'   `x` and `y`, such as `attr(distances, "receivers")` or
#'   `stoney_track$receivers`.
#' @param count `NULL`, or the name of a column of `data` that holds each
#'   receiver's detection count.
#' @param size Size of a ring of fixed size, as `size` in
#'   [ggplot2::geom_point()].
#' @param max_size Size of the ring of the largest count.
#' @param stroke Width of the coloured ring, in the units `ggplot2` uses for
#'   `stroke`; the dark ring under it is twice as wide. The default 1.1 suits
#'   one map; use about 0.5 on small multiples so the rings do not cover the
#'   cells beneath them.
#' @return A list of ggplot2 layers, and with `count` the size and shape
#'   scales of their keys, to add to a plot with `+`. The layers do not
#'   inherit the plot's aesthetics.
#' @seealso [theme_positionr()], [scale_fill_depth()] and
#'   [scale_fill_probability()], which draw the rest of a positionR map.
#' @export
#' @examples
#' # The 70 receivers of the Stoney example on the depth raster.
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' lake <- terra::as.data.frame(abs(terra::rast(tif)), xy = TRUE)
#' names(lake)[3] <- "depth"
#' receivers <- stoney_track$receivers
#' ggplot2::ggplot(lake, ggplot2::aes(x, y, fill = depth)) +
#'   ggplot2::geom_raster() +
#'   scale_fill_depth() +
#'   layer_receivers(receivers) +
#'   ggplot2::coord_equal() +
#'   theme_positionr()
#'
#' # Detections of one hour: two receivers heard the animal, the rest none.
#' receivers$detections <- 0
#' receivers$detections[receivers$station_id %in% c("S019", "S020")] <- c(5, 20)
#' ggplot2::ggplot(lake, ggplot2::aes(x, y, fill = depth)) +
#'   ggplot2::geom_raster() +
#'   scale_fill_depth() +
#'   layer_receivers(receivers, count = "detections") +
#'   ggplot2::coord_equal() +
#'   theme_positionr()
layer_receivers <- function(data, count = NULL, size = 2.6, max_size = 6,
                            stroke = 1.1) {
  if (!is.data.frame(data) || !all(c("x", "y") %in% names(data))) {
    stop("`data` must be a data frame with columns `x` and `y`.", call. = FALSE)
  }
  # An empty table would draw no receiver and say nothing, as a merge on the
  # wrong columns gives.
  if (!nrow(data)) stop("`data` has no rows, so there is no receiver to draw.", call. = FALSE)
  ring <- positionr_ring_colours()
  if (is.null(count)) return(ring_layers(data, ring[["heard"]], size, stroke = stroke))
  counts <- receiver_counts(data, count)
  heard <- data[counts > 0, , drop = FALSE]
  silent <- data[counts == 0, , drop = FALSE]
  heard$.count <- heard[[count]]
  silent$.key <- rep("No detections", nrow(silent))
  list(
    if (nrow(silent)) ring_layers(silent, ring[["silent"]], size, key = TRUE, stroke = stroke),
    if (nrow(heard)) ring_layers(heard, ring[["heard"]], NULL, sized = TRUE, stroke = stroke),
    if (nrow(heard)) {
      ggplot2::scale_size_area(name = "Detections", max_size = max_size,
                               limits = c(0, NA), breaks = count_breaks,
                               guide = ggplot2::guide_legend(order = 1))
    },
    if (nrow(silent)) {
      ggplot2::scale_shape_manual(name = NULL, values = c("No detections" = 21),
                                  guide = ggplot2::guide_legend(order = 2))
    }
  )
}

# The column `count` of `data`, checked to be counts of 0 or more.
receiver_counts <- function(data, count) {
  if (!is.character(count) || length(count) != 1 || !count %in% names(data)) {
    stop("`count` must be NULL or the name of a column of `data`.", call. = FALSE)
  }
  counts <- data[[count]]
  if (!is.numeric(counts) || anyNA(counts) || any(counts < 0)) {
    stop(sprintf("Column `%s` must hold counts of 0 or more.", count), call. = FALSE)
  }
  counts
}

# Two point layers: the dark ring, then the coloured ring over it, in the
# colours of positionr_ring_colours() (theme_positionr.R). A sized ring maps
# size to `.count`; a keyed ring maps shape to `.key`, so the shape scale
# keys it; any other ring has shape 21 and size `size`.
ring_layers <- function(rings, colour, size, sized = FALSE, key = FALSE, stroke = 1.1) {
  mapping <- if (sized) {
    ggplot2::aes(x, y, size = .count)
  } else if (key) {
    ggplot2::aes(x, y, shape = .key)
  } else {
    ggplot2::aes(x, y)
  }
  fixed <- list(fill = NA)
  if (!sized) fixed$size <- size
  if (!key) fixed$shape <- 21
  # The dark ring is twice the stroke of the coloured ring, so it shows one
  # step wider on each side.
  edge <- positionr_ring_colours()[["edge"]]
  lapply(list(c(edge, 2 * stroke), c(colour, stroke)), function(ring) {
    params <- c(fixed, list(colour = ring[1], stroke = as.numeric(ring[2])))
    ggplot2::layer(geom = "point", stat = "identity", position = "identity",
                   data = rings, mapping = mapping, inherit.aes = FALSE,
                   show.legend = if (sized || key) NA else FALSE,
                   params = params)
  })
}
