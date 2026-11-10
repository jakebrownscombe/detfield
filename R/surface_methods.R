# ---- Methods for location surfaces ----

#' @export
print.loc_surface <- function(x, ...) {
  b <- attr(x, "bins")
  od <- attr(x, "overdispersion")
  n_animals <- length(unique(b$animal_id))
  cat_wrapped(sprintf("Location surface: %s, %s, %s cells kept in all\n",
                      plural(n_animals, "animal"), plural(nrow(b), "bin"),
                      fmt_n(nrow(x))))
  if (nrow(b) && !is.null(b$n_cells)) {
    cat_wrapped(sprintf("Cells kept per bin: median %s of %s water cells (range %s to %s)\n",
                        fmt_n(stats::median(b$n_cells)),
                        fmt_n(attr(x, "n_water_cells")),
                        fmt_n(min(b$n_cells)),
                        fmt_n(max(b$n_cells))))
  }
  spread <- if (identical(attr(x, "likelihood"), "presence_only")) {
    "not used under the presence-only reading"
  } else if (od$form == "none") {
    "none (overdispersion = \"none\")"
  } else {
    sprintf("%s per receiver-hour (sigma$%s), integrated over %s",
            plural(signif(od$sigma, 4), "logit"), od$form, plural(od$nodes, "quadrature point"))
  }
  cat_wrapped(sprintf("Likelihood: %s\n", likelihood_label(attr(x, "likelihood"))))
  cat_wrapped(sprintf("Spread: %s\n", spread))
  cat_wrapped(sprintf("Each bin keeps the smallest set of cells holding %s of its probability\n",
                      format(attr(x, "kept")$mass)))
  if (!is.null(b$misfit_p)) {
    cat_wrapped(sprintf(
      paste0("Misfit: %s of %s bins have a detecting receiver whose count the ",
             "surface gives a chance below 0.01\n"),
      fmt_n(sum(b$misfit_p < 0.01, na.rm = TRUE)),
      fmt_n(nrow(b))
    ))
  }
  near_line(b, attr(x, "near_rule"))
  skipped <- attr(x, "skipped")
  if (!is.null(skipped) && nrow(skipped)) {
    cat_wrapped(sprintf("Left out: %s (attr(, \"skipped\"))\n",
                        plural(nrow(skipped), "degenerate bin")))
  }
  lost <- attr(x, "not_listening_bins")
  if (!is.null(lost) && nrow(lost)) {
    cat_wrapped(sprintf("Left out: %s heard only at receivers not listening (attr(, %s))\n",
                        plural(nrow(lost), "bin"), "\"not_listening_bins\""))
  }
  print(utils::head(data.frame(animal_id = x$animal_id, bin_start = x$bin_start,
                               cell = x$cell, prob = x$prob), n = 6))
  invisible(x)
}

likelihood_label <- function(likelihood) {
  if (identical(likelihood, "presence_only")) {
    return("presence-only (product of p over detections; misses not modelled)")
  }
  "binomial (detections and misses at every listening receiver)"
}

# Subsetting keeps the attributes, with the bin table cut to the bins left
# and each bin's stored mass and cell count recounted from the rows left.
#' @export
`[.loc_surface` <- function(x, i, j, ...) {
  keep <- attributes(x)[c("bins", "overdispersion", "likelihood", "kept", "grid", "n_water_cells",
                          "key",
                          "receivers", "receiver_xy", "near_rule", "skipped",
                          "not_listening_bins")]
  out <- NextMethod()
  if (!is.data.frame(out)) return(out)
  if (!all(c("animal_id", "bin_start", "cell", "prob") %in% names(out))) {
    class(out) <- "data.frame"
    return(out)
  }
  left <- paste(out$animal_id, as.numeric(out$bin_start))
  for (a in c("bins", "receivers")) {
    tab <- keep[[a]]
    if (!is.null(tab)) keep[[a]] <- tab[paste(tab$animal_id, as.numeric(tab$bin_start)) %in% left, ]
  }
  if (!is.null(keep$bins) && nrow(keep$bins)) {
    bkey <- paste(keep$bins$animal_id, as.numeric(keep$bins$bin_start))
    keep$bins$kept_mass <- as.numeric(tapply(out$prob, left, sum)[bkey])
    keep$bins$n_cells <- as.integer(table(left)[bkey])
  }
  attributes(out)[names(keep)] <- keep
  class(out) <- c("loc_surface", "data.frame")
  out
}

#' Map one location surface
#'
#' Draws the probability surface of one animal and bin on its grid, with its
#' credible regions outlined and the receivers listening in the bin marked.
#'
#' @details
#' The surface is drawn with [scale_fill_probability()], the magma scale of
#' every positionR probability map, over the stored cells. The scale runs from
#' 0 to the largest probability of a cell in the bin, which is the top label of
#' the key, because each cell holds a small share of the animal's position and
#' on 0 to 1 the whole map would be black. Other water cells are grey and land
#' is blank. The title names the animal and the start of the bin on one line
#' above the map and its keys. Each credible region is drawn as one outline of
#' its cells, in a colour outside the magma scale over a thin dark line so it
#' reads on the bright and the dark ends of the scale: with two levels, the
#' smaller (by default 50%) is sky blue and the larger (95%) bluish green; one
#' level is drawn bluish green. The key is titled "Credible region". When the
#' larger region spans under a tenth of the map's longer side, both rings are
#' drawn at half width, so the cells inside a small region still show.
#' The receivers listening in the bin are drawn over the regions by
#' [layer_receivers()], the receiver mark of every positionR map: a ring over
#' a dark ring one step wider, never filled, so the surface shows through. A
#' receiver that detected the animal is a yellow ring whose area is
#' proportional to its count, from zero, up to 6 mm across at the largest
#' count of the bin; the key "Detections" shows at most three counts. A
#' receiver that heard none is a white ring, keyed "No detections". Only the
#' receivers inside the map are drawn and keyed, so a key never names a mark
#' the map does not show.
#'
#' The map zooms to the stored cells and the receivers that detected the
#' animal, a quarter wider on each side; silent receivers near that box (within
#' its longer side) widen it, so the map shows where the animal was not heard
#' without zooming out to a whole array. On a grid in metres each side of the
#' map is at least 1.5 km and at least three quarters of the longer side. Give
#' `ext` to draw another extent. The figure is built with ggplot2 and
#' [theme_positionr()], with equal coordinates and a scale bar.
#'
#' @param x A location surface from [location_surface()].
#' @param animal_id The animal to draw; the first in the surface by default.
#' @param bin_start Start of the bin to draw (POSIXct, or text read as UTC); the
#'   animal's first bin by default.
#' @param level One or two levels of the credible regions outlined, each
#'   above 0 and at most 1, or `NULL` for none. The default draws the 50% and
#'   95% regions.
#' @param ... `ext`, the extent to draw, as a [terra::ext()] or as
#'   `c(xmin, xmax, ymin, ymax)`. Other arguments are ignored with a warning;
#'   change the figure by adding ggplot2 layers, scales or [ggplot2::labs()]
#'   to the plot it returns. On a grid in metres the axes are titled
#'   "Easting (m)" and "Northing (m)", with the grid's coordinate reference
#'   system named in the brackets when the grid carries one: "UTM 17N" for
#'   a UTM zone, otherwise its authority code, such as "EPSG:3161".
#' @return A ggplot object. It prints when not assigned; assign it to add
#'   layers, such as a track with [ggplot2::geom_path()], or to save it with
#'   [ggplot2::ggsave()].
#' @export
#' @examples
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Beyond 2 km the Stoney curve is below 0.001 at every receiver; leaving
#' # those pairs out keeps the example quick.
#' dist <- water_distances(grid, stoney_track$receivers,
#'                         max_distance = 2000)
#' model <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' field <- de_field(dist, de_model(model, sigma = list(station_hour = 0)))
#' # Position one animal in one hour, the hour the map draws
#' dets <- subset(stoney_track$detections, animal_id == "stoney-01" &
#'                  format(datetime, "%d %H", tz = "UTC") == "01 07")
#' counts <- bin_detections(dets, stoney_track$receivers, stoney_track$tags)
#' surface <- location_surface(counts, field)
#' plot(surface)
#' plot(surface, level = 0.95)
plot.loc_surface <- function(x, animal_id = NULL, bin_start = NULL, level = c(0.5, 0.95), ...) {
  grid <- sl_grid(x, strict = TRUE)
  level <- pl_levels(level)
  if (is.null(grid)) stop("The surface holds no grid to draw on.", call. = FALSE)
  args <- list(...)
  given <- if (is.null(names(args))) rep("", length(args)) else names(args)
  do.call(unused_dots, args[given != "ext"])
  pick <- pl_pick(x, animal_id, bin_start)
  one <- x[x$animal_id == pick$animal_id & x$bin_start == pick$bin_start, ]
  rx <- pl_receivers_of(x, pick)
  template <- terra::rast(grid, nlyrs = 1)
  e <- as.vector(if (is.null(args$ext)) pl_zoom(template, one$cell, rx) else terra::ext(args$ext))
  titles <- axis_titles(template)
  ggplot2::ggplot(pl_cells(template, grid, one, e), ggplot2::aes(x, y)) +
    ggplot2::geom_raster(fill = "grey85") +
    ggplot2::geom_raster(ggplot2::aes(fill = prob), na.rm = TRUE) +
    pl_region(one, template, level, e) +
    pl_points(rx, e) +
    pl_fill(max(one$prob)) +
    ggplot2::coord_equal(xlim = e[c("xmin", "xmax")], ylim = e[c("ymin", "ymax")]) +
    ggplot2::labs(title = sprintf("%s, %s UTC", pick$animal_id,
                                  format(pick$bin_start, "%Y-%m-%d %H:%M", tz = "UTC")),
                  x = titles[["x"]], y = titles[["y"]]) +
    annotation_scale_bar(width_hint = 0.25) +
    theme_positionr() +
    # The title sits on its own line above the panel and the keys, so the two
    # cannot overlap; keys on the right leave the height to a tall map.
    ggplot2::theme(legend.position = "right", legend.box.just = "left",
                   legend.spacing.y = grid::unit(2, "pt"),
                   plot.title.position = "plot")
}

# The fill of a surface map: magma from 0 to the bin's largest cell
# probability, with that maximum as the top label of the key, so every
# surface is read the same way and the key states the maximum.
pl_fill <- function(top) {
  if (!is.finite(top) || top <= 0) top <- 1
  scale_fill_probability(name = "Probability\nper cell", limits = c(0, top),
                         breaks = c(0, top / 2, top),
                         labels = format(signif(c(0, top / 2, top), 2), drop0trailing = TRUE),
                         guide = ggplot2::guide_colourbar(order = 4))
}

# Water cells in and near extent `e`, with the probability of the stored
# cells and NA elsewhere. Water beyond the map and a margin is not drawn, so a
# bin in a large lake does not carry every cell of the lake. Both raster
# layers share these cells, so a surface stored in one column of cells is
# still drawn at the grid's resolution.
pl_cells <- function(template, grid, one, e) {
  cells <- which(!is.na(terra::values(grid, mat = FALSE)))
  water <- data.frame(terra::xyFromCell(template, cells))
  water$prob <- one$prob[match(cells, one$cell)]
  pad <- 0.1 * c(e[["xmax"]] - e[["xmin"]], e[["ymax"]] - e[["ymin"]]) + terra::res(template)
  water[water$x >= e[["xmin"]] - pad[1] & water$x <= e[["xmax"]] + pad[1] &
          water$y >= e[["ymin"]] - pad[2] & water$y <= e[["ymax"]] + pad[2], ]
}

# Extent of the stored cells and the receivers that detected the animal, a
# quarter wider on each side, so one bin is not lost in a large lake. Silent
# receivers within the longer side of that box widen it, so the map shows
# where the animal was not heard; silent receivers further off are left out,
# because on a large array they would zoom the map out to the whole water body
# and draw the region a few pixels wide.
pl_zoom <- function(template, cells, rx) {
  xy <- terra::xyFromCell(template, cells)
  if (!is.null(rx)) {
    heard <- rx$y_count > 0
    xy <- rbind(xy, as.matrix(rx[heard, c("x", "y")]))
    reach <- max(diff(range(xy[, 1])), diff(range(xy[, 2])), 10 * terra::res(template))
    near <- !heard & rx$x >= min(xy[, 1]) - reach & rx$x <= max(xy[, 1]) + reach &
      rx$y >= min(xy[, 2]) - reach & rx$y <= max(xy[, 2]) + reach
    xy <- rbind(xy, as.matrix(rx[near, c("x", "y")]))
  }
  pad <- pmax(0.25 * c(diff(range(xy[, 1])), diff(range(xy[, 2]))), 5 * terra::res(template))
  lo <- c(min(xy[, 1]), min(xy[, 2])) - pad
  hi <- c(max(xy[, 1]), max(xy[, 2])) + pad
  # On a grid in metres each side is at least 1.5 km and at least three
  # quarters of the longer side, so a small region is not drawn in a strip.
  if (isTRUE(terra::linearUnits(template) == 1)) {
    side <- pmax(hi - lo, min_side_m, 0.75 * max(hi - lo))
    mid <- (lo + hi) / 2
    lo <- mid - side / 2
    hi <- mid + side / 2
  }
  terra::intersect(terra::ext(lo[1], hi[1], lo[2], hi[2]), terra::ext(template))
}

# Shortest side of a surface map on a grid in metres.
min_side_m <- 1500

pl_pick <- function(x, animal_id, bin_start) {
  b <- attr(x, "bins")
  if (is.null(animal_id)) animal_id <- b$animal_id[1]
  b <- b[b$animal_id == animal_id, ]
  if (!nrow(b)) stop(sprintf("Animal %s is not in the surface.", animal_id), call. = FALSE)
  if (is.null(bin_start)) return(list(animal_id = animal_id, bin_start = b$bin_start[1]))
  if (is.character(bin_start)) bin_start <- as.POSIXct(bin_start, tz = "UTC")
  hit <- which(as.numeric(b$bin_start) == as.numeric(bin_start))
  if (!length(hit)) {
    stop(sprintf("Animal %s has no bin starting %s UTC.", animal_id, ls_time(bin_start)),
         call. = FALSE)
  }
  list(animal_id = animal_id, bin_start = b$bin_start[hit[1]])
}

# The levels to outline: NULL, or one or two numbers in (0, 1], sorted.
pl_levels <- function(level) {
  if (is.null(level)) return(NULL)
  ok <- is.numeric(level) && length(level) %in% 1:2 && !anyNA(level) && all(level > 0) &&
    all(level <= 1)
  if (!ok) {
    stop("`level` must be one or two numbers above 0 and at most 1, or NULL.", call. = FALSE)
  }
  sort(unique(level))
}

# Region colours, from the Okabe-Ito set: neither hue is on the magma scale,
# and neither is the yellow of the detecting receivers. With two levels the
# smaller is sky blue and the larger bluish green; one level is bluish green.
region_colours <- function(n) if (n == 1) "#009E73" else c("#56B4E9", "#009E73")

# The outlines of the credible regions, one ring set per level, each drawn as
# a coloured line over a dark line twice as wide, so it reads on the
# near-black and the pale-yellow ends of magma. When the larger region spans
# under a tenth of the longer side of extent `e`, both lines are drawn at
# half width: on a dense array the default view holds every detecting
# receiver, and full-width lines would hide a region of about 30 cells under
# its own rings. NULL for none.
pl_region <- function(one, template, level, e = NULL) {
  if (is.null(level)) return(NULL)
  o <- order(one$prob, decreasing = TRUE)
  before <- cumsum(one$prob[o]) - one$prob[o]
  labels <- sprintf("%g%%", 100 * level)
  rings <- do.call(rbind, lapply(seq_along(level), function(i) {
    r <- cell_outline(template, one$cell[o][before < level[i] * (1 - 1e-12)])
    if (!nrow(r)) return(NULL)
    r$region <- labels[i]
    r$ring <- paste(labels[i], r$ring)
    r
  }))
  if (is.null(rings) || !nrow(rings)) return(NULL)
  rings$region <- factor(rings$region, levels = labels)
  width <- ring_width(rings, e)
  list(
    # A thin dark edge either side of the colour, so the outline reads on both
    # ends of magma without covering the cells.
    ggplot2::geom_path(data = rings, ggplot2::aes(x, y, group = ring), colour = "grey10",
                       linewidth = width + 0.5, inherit.aes = FALSE),
    ggplot2::geom_path(data = rings, ggplot2::aes(x, y, group = ring, colour = region),
                       linewidth = width, inherit.aes = FALSE),
    ggplot2::scale_colour_manual(name = "Credible region",
                                 values = stats::setNames(region_colours(length(level)), labels),
                                 breaks = labels, drop = FALSE,
                                 guide = ggplot2::guide_legend(order = 3))
  )
}

# Width of the coloured ring line: 0.6, or 0.3 when the rings span under a
# tenth of the longer side of the map extent `e`.
ring_width <- function(rings, e) {
  if (is.null(e)) return(0.6)
  span <- max(diff(range(rings$x)), diff(range(rings$y)))
  view <- max(e[["xmax"]] - e[["xmin"]], e[["ymax"]] - e[["ymin"]])
  if (span < 0.1 * view) 0.3 else 0.6
}

# Axis titles of a map: "Easting (m, UTM 17N)" and "Northing (m, UTM 17N)" on
# a grid in metres, with the coordinate reference system named as short as
# it reads unambiguously; "(m)" alone when the grid names none; "x"
# and "y" on a grid not in metres.
axis_titles <- function(template) {
  if (!isTRUE(terra::linearUnits(template) == 1)) return(c(x = "x", y = "y"))
  crs <- crs_label(template)
  unit <- if (is.null(crs)) "(m)" else sprintf("(m, %s)", crs)
  c(x = paste("Easting", unit), y = paste("Northing", unit))
}

# "UTM 17N" for a UTM zone, else "AUTHORITY:code", else the CRS name; NULL
# when the grid carries no CRS.
crs_label <- function(template) {
  d <- tryCatch(terra::crs(template, describe = TRUE), error = function(e) NULL)
  name <- if (is.null(d) || !nrow(d)) "" else paste0(d$name[1])
  if (name %in% c("", "NA", "unknown")) return(NULL)
  zone <- regmatches(name, regexec("UTM zone ([0-9]+[NS])", name))[[1]]
  if (length(zone) == 2) return(paste("UTM", zone[2]))
  code <- paste0(d$authority[1], ":", d$code[1])
  if (grepl("^[^:]+:[^:]+$", code) && !grepl("NA", code, fixed = TRUE)) code else name
}

# The boundary rings of a set of cells, one row per vertex, with `ring`
# naming each closed ring. The cells are dissolved into one shape. Cells that
# touch only at a corner would stay separate polygons, each with its own
# outline, so on a grid in metres the shape is widened by 1% of a cell, which
# joins them and moves the outline by less than a pixel.
cell_outline <- function(template, cells) {
  v <- rep(NA_real_, terra::ncell(template))
  v[cells] <- 1
  outline <- terra::as.polygons(terra::setValues(template, v), dissolve = TRUE)
  if (!nrow(outline)) return(data.frame(x = numeric(), y = numeric(), ring = character()))
  if (isTRUE(terra::linearUnits(template) == 1)) {
    shape <- sf::st_union(sf::st_as_sf(outline))
    shape <- sf::st_buffer(shape, 0.01 * min(terra::res(template)), joinStyle = "MITRE",
                           mitreLimit = 2)
    g <- as.data.frame(sf::st_coordinates(shape))
    parts <- intersect(c("L1", "L2", "L3"), names(g))
    return(data.frame(x = g$X, y = g$Y, ring = do.call(paste, g[parts])))
  }
  g <- as.data.frame(terra::geom(outline))
  data.frame(x = g$x, y = g$y, ring = paste(g$geom, g$part, g$hole))
}

# Listening receivers of one bin with their coordinates and counts, or NULL
# when the surface holds no coordinates.
pl_receivers_of <- function(x, pick) {
  xy <- attr(x, "receiver_xy")
  rows <- attr(x, "receivers")
  if (is.null(xy) || is.null(rows)) return(NULL)
  rows <- rows[rows$animal_id == pick$animal_id & rows$bin_start == pick$bin_start, ]
  at <- xy[match(rows$receiver, xy$receiver), ]
  data.frame(x = at$x, y = at$y, y_count = rows$y)
}

# The listening receivers of the bin inside extent `e`, drawn by
# layer_receivers(): yellow rings sized by count for those that heard the
# animal, white rings for those that heard none, each over a dark ring.
# Receivers outside the map are left out, so the keys name only marks the
# map shows. NULL when the surface holds no receiver coordinates.
pl_points <- function(rx, e) {
  if (is.null(rx)) return(NULL)
  rx <- rx[rx$x >= e[["xmin"]] & rx$x <= e[["xmax"]] &
             rx$y >= e[["ymin"]] & rx$y <= e[["ymax"]], , drop = FALSE]
  if (!nrow(rx)) return(NULL)
  layer_receivers(rx, count = "y_count")
}

# At most three whole counts above zero and within the limits, so each ring
# of the key differs from the next by a visible step.
count_breaks <- function(lim) {
  top <- max(lim, na.rm = TRUE)
  if (!is.finite(top) || top <= 0) return(numeric())
  b <- pretty(c(0, top), n = 3)
  b <- b[b > 0 & b <= top & b == round(b)]
  if (!length(b)) b <- ceiling(top)
  if (length(b) > 3) b <- b[seq(length(b) - 2L, length(b))]
  b
}

# How many bins put almost no probability near the receiver that heard the
# animal most, and how far the modes sit from the detecting receivers; the
# misfit count does not show a surface that sits in a ring.
near_line <- function(b, rule) {
  if (is.null(b$near_mass) || !nrow(b)) return(invisible(NULL))
  scored <- !is.na(b$near_mass)
  if (!any(scored)) {
    cat(paste0("Near: not computed; the field gives no d50 (no distance table, no receiver ",
               "positions, or no pair at p >= 0.5); set control$near_m\n"))
    return(invisible(NULL))
  }
  radius <- unique(b$near_m[scored])
  # Without control$near_m, the radius is the distance to the farthest water
  # cell centre at which the field gives that receiver p >= 0.5. It is read
  # off cell centres, so it falls short of the curve's d50 by up to a cell.
  label <- if (identical(rule, "control$near_m")) "control$near_m" else
    "per receiver, its farthest cell centre at p >= 0.5"
  within <- if (length(radius) == 1) sprintf("%s m", format(signif(radius, 4))) else
    sprintf("%s to %s m", format(signif(min(radius), 4)), format(signif(max(radius), 4)))
  cat_wrapped(sprintf(
    paste0("Near: %s of %s bins put under 0.001 of their probability within %s of ",
           "the receiver with the most detections (%s)\n"),
    fmt_n(sum(b$near_mass[scored] < 0.001)),
    fmt_n(sum(scored)), within, label
  ))
  d <- b$mode_dist_m[!is.na(b$mode_dist_m)]
  if (length(d)) {
    cat_wrapped(sprintf(
      paste0("Modes sit a median %s m (max %s m) from the nearest detecting ",
             "receiver, over %s\n"),
      fmt_n(round(stats::median(d))),
      fmt_n(round(max(d))), plural(length(d), "bin")
    ))
  }
  invisible(NULL)
}
