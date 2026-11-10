#' Distances through water from receivers to every water cell
#'
#' Measures, for every receiver and water cell, how far sound travels from the
#' receiver to the cell centre without crossing land. Where the straight line
#' stays in water, the distance is that line. Where land lies on it, the
#' distance is the shortest path through water: straight legs that bend only
#' at corners where water turns round land. The corners form a graph, with an
#' edge between any two that see each other across water, and on that graph
#' the algorithm of Dijkstra (1959) builds the minimal paths from the receiver
#' in order of increasing length. The path starts at the receiver itself, so
#' `distance_m` is never less than `straight_m`. In the Stoney example the
#' cell 640 m north-west of station S023 lies across a point of land; its
#' path is 1,955 m.
#'
#' The detection curve is read at this distance whether or not land lies on
#' the line; [de_model()] has a `blocked` argument for such pairs.
#'
#' **The land rule.** One rule tests every straight leg. A leg is blocked when
#' it comes more than a thousandth of a cell side (10 cm on 100 m cells) into
#' land. A point where two land cells meet diagonally counts as land, so a leg
#' through it is blocked. A leg that grazes land by less is clear. The margin
#' gives a receiver a few millimetres off its cell centre the same lines of
#' sight as one at the centre.
#'
#' **Sight one way.** A receiver within that margin of such a point can see
#' cells on both sides of it that cannot see it back. No distance changes.
#'
#' **The unreachable cell.** Two water cells that touch only where two land
#' cells meet are not joined. In the Stoney example one water cell joins the
#' lake only that way. It gets no row and is counted in the `unreachable`
#' attribute.
#'
#' **Memory.** Lines of sight from the corners are tested in fixed batches of
#' about 200,000 lines. The tests run in passes. With one core, each receiver
#' in turn makes one pass over the corners it reaches that no earlier receiver
#' reached. With several cores, each core makes one pass over its share of the
#' corners within `max_distance` of any receiver. On the Stoney example with
#' no `max_distance` and one core, where all lines are tested in one pass,
#' batches of 200,000 lines give a lower peak memory than batches ten times
#' larger, at about the same time. A pass of 200,000 lines or fewer makes the
#' same batches at either size. The script that measures both sizes on that
#' example, each run in its own R process, is
#' <https://github.com/jakebrownscombe/detfield/blob/main/validation/04_memory_time.R>.
#' Its log is
#' <https://github.com/jakebrownscombe/detfield/blob/main/validation/outputs/04_memory_time.log>.
#' Both are in the repository and not in the installed package.
#'
#' **Snapping.** Cells with a value are water and `NA` cells are land. Rows of
#' `receivers` that share a key are averaged to one position. A receiver on a
#' land cell moves to the centre of the nearest water cell when it lies within
#' one cell diagonal of the edge of some water cell; one farther from water is
#' an error. The move is at least half a cell side and at most one and a half
#' cell diagonals: 50 to 212 m on 100 m cells. The warning gives each move in
#' metres, and the `receivers` attribute keeps it as `snap_m`.
#'
#' @param grid A [terra::SpatRaster] in a projected CRS with metre units, such
#'   as the UTM zone of the lake, or a file path that [terra::rast()] reads.
#'   Only the first layer is used.
#' @param receivers A data frame with a `station_id` or `receiver_serial`
#'   column and either `x` and `y` in the CRS of `grid` or `deploy_lat` and
#'   `deploy_lon` in WGS 84; or an sf point object. When `x` and `y` carry an
#'   `epsg` attribute it must match the grid. Other columns pass through to
#'   the `receivers` attribute; with both `receiver_serial` and `station_id`
#'   the rows also form the `lookup` attribute.
#' @param max_distance The largest distance through water kept, in metres,
#'   the unit of the grid's projected CRS. A pair farther apart through water
#'   is left out. The default, `Inf`, keeps every pair connected by water.
#' @param control A list of rarely changed settings:
#'   * `cores`: number of processes for [parallel::mclapply()], one receiver
#'     each. The default is 1. Windows has no forked processes, so there 1
#'     is used, with a warning when more are asked for.
#'   * `path`: `"exact"`, the default, for the shortest path through water;
#'     `"grid"` to keep the straight line on clear pairs and take the
#'     8-neighbour grid path where land lies on the line (see NEWS).
#'   * `maxiter`: passed to [terra::gridDist()] when `path = "grid"`. The
#'     default is 1000.
#'
#' @return A data frame of class `water_distances`, one row per receiver and
#'   water cell within `max_distance`:
#'   * `cell`: the terra cell number of the water cell.
#'   * `station_id` or `receiver_serial`: the receiver key.
#'   * `distance_m`: the distance in metres used for detection; never less
#'     than `straight_m`. These are metres on the ground. A grid is refused
#'     when its CRS scale at the raster centre is more than 1 percent from 1
#'     on either axis or either diagonal, or cannot be measured, as in a local
#'     engineering CRS.
#'   * `straight_m`: the straight-line distance in metres.
#'   * `in_water`: `TRUE` where the straight line stays in water, so
#'     `distance_m` equals `straight_m`; `FALSE` where land lies on it.
#'
#'   Attributes: `grid` (the grid, packed by [terra::wrap()]), `cells` (all
#'   water cells), `key`, `receivers` (one row per key with the position used,
#'   `snap_m`, `n_rows` and `spread_m`), `lookup`, `max_distance` and
#'   `unreachable` (pairs that no water path connects, including water joined
#'   to the receiver's only through a point where two land cells meet).
#'
#' @references
#' Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
#' *Numerische Mathematik*, 1, 269-271. \doi{10.1007/BF01386390}
#' @seealso [de_field()], which turns the distances into detection
#'   probabilities.
#' @export
#' @examples
#' # The Stoney example: Stoney Lake (100 m cells, land NA) and 70 receivers.
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Pairs beyond 2 km are left out, which keeps the example fast.
#' d <- water_distances(grid, stoney_track$receivers, max_distance = 2000)
#' d
#' # The cell 640 m north-west of S023 lies across a point of land; the path
#' # goes round it, 1,955 m against 640 m in a straight line.
#' north_west <- terra::cellFromXY(grid, cbind(724054, 4936291))
#' d[d$station_id == "S023" & d$cell == north_west, ]
water_distances <- function(grid, receivers, max_distance = Inf, control = list()) {
  grid <- as_water_grid(grid, call = sys.call())
  if (is.data.frame(receivers) && !nrow(receivers)) {
    stop("`receivers` has no rows; give at least one receiver.", call. = FALSE)
  }
  ctrl <- distance_control(control)
  if (!is_one_number(max_distance) || max_distance <= 0) {
    stop("`max_distance` must be one positive number.", call. = FALSE)
  }
  check_land_coding(grid)
  info <- grid_info(grid)
  if (!any(info$water)) stop("`grid` has no water cells: every cell is NA.", call. = FALSE)
  rx_in <- receiver_xy(receivers, grid)
  key <- receiver_key(rx_in)
  rx <- snap_receivers(collapse_receivers(rx_in, key, info$res), info, key)
  cells <- which(info$water)
  packed <- terra::wrap(grid)
  wp <- NULL
  if (ctrl$path == "exact") {
    wp <- water_paths(info, max_distance)
    fill_sight(wp, info, rx, ctrl$cores)
  }
  one <- function(i) {
    one_receiver_distances(packed, info, cells, rx[i, ], ctrl, max_distance, wp)
  }
  parts <- parallel_apply(seq_len(nrow(rx)), one, ctrl$cores)
  n_unreachable <- sum(vapply(parts, attr, numeric(1), which = "unreachable"))
  out <- data.frame(
    cell = unlist(lapply(parts, `[[`, "cell"), use.names = FALSE),
    key = rep(rx[[key]], vapply(parts, function(p) length(p$cell), integer(1))),
    distance_m = unlist(lapply(parts, `[[`, "distance_m"), use.names = FALSE),
    straight_m = unlist(lapply(parts, `[[`, "straight_m"), use.names = FALSE),
    in_water = unlist(lapply(parts, `[[`, "in_water"), use.names = FALSE),
    stringsAsFactors = FALSE
  )
  names(out)[2] <- key
  structure(out, class = c("water_distances", "data.frame"), grid = packed, cells = cells,
            key = key, receivers = rx, lookup = receiver_lookup(rx_in),
            max_distance = max_distance, unreachable = n_unreachable)
}

distance_control <- function(control) {
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  unknown <- setdiff(names(control), c("cores", "maxiter", "path"))
  if (length(unknown)) {
    stop(sprintf("Unknown `control` %s: %s.", agree(length(unknown), "setting", "settings"),
                 paste(unknown, collapse = ", ")),
         call. = FALSE)
  }
  ctrl <- utils::modifyList(list(cores = 1L, maxiter = 1000, path = "exact"), control)
  ctrl$cores <- max(1L, as.integer(ctrl$cores))
  if (!is_one_string(ctrl$path) || !ctrl$path %in% c("exact", "grid")) {
    stop("`control$path` must be \"exact\" or \"grid\".", call. = FALSE)
  }
  ctrl
}

one_receiver_distances <- function(packed, info, cells, rx, ctrl, max_distance, wp = NULL) {
  xy <- xy_of_cell(info, cells)
  straight <- sqrt((xy[, "x"] - rx$x)^2 + (xy[, "y"] - rx$y)^2)
  near <- which(straight <= max_distance)
  clear <- rep(FALSE, length(cells))
  clear[near] <- line_clear(info, rx$x, rx$y, xy[near, "x"], xy[near, "y"])
  if (ctrl$path == "grid") {
    path <- grid_path(packed, info, cells, rx$cell, ctrl$maxiter)
    # The grid path starts at the centre of the receiver's cell and the
    # straight line at the receiver, so near the receiver it can come out shorter.
    distance <- ifelse(clear, straight, pmax(path, straight))
    unreachable <- sum(is.na(distance))
  } else {
    distance <- rep(NA_real_, length(cells))
    distance[clear] <- straight[clear]
    hidden <- near[!clear[near]]
    distance[hidden] <- water_path_from(wp, info, rx$x, rx$y)[hidden]
    distance[!is.finite(distance)] <- NA
    # Water no path reaches: another water body, joined at most through a
    # closed gap, that no line of sight reaches either.
    own <- wp$body[match(rx$cell, cells)]
    unreachable <- sum(wp$body != own & is.na(distance))
  }
  keep <- which(!is.na(distance) & distance <= max_distance)
  structure(
    list(cell = cells[keep], distance_m = distance[keep], straight_m = straight[keep],
         in_water = clear[keep]),
    unreachable = unreachable
  )
}

# Distance along the 8-neighbour grid from the centre of `from_cell` to every
# cell, by terra::gridDist(). NA where no grid path reaches.
grid_path <- function(packed, info, cells, from_cell, maxiter) {
  grid <- terra::unwrap(packed)
  v <- ifelse(info$water, 1, NA_real_)
  v[from_cell] <- 0
  target <- terra::setValues(terra::rast(grid), v)
  path <- terra::values(terra::gridDist(target, target = 0, maxiter = maxiter),
                        mat = FALSE)[cells]
  path[!is.finite(path)] <- NA
  path
}

# TRUE where the segment from (x0, y0) to (x, y) stays in water, by one
# geometric rule. Each land cell is shrunk by the
# tolerance on every side that faces water, and at every corner whose
# diagonal neighbour is water while the two cells beside it are land: what is
# left is the land a segment may not touch. Land cells that meet along an edge
# stay joined, so a segment along or just beside that edge is blocked. The
# point where two land cells meet diagonally (a closed gap) is a square of
# half-side the tolerance that a segment may not touch either, unless it
# starts inside that square: a receiver placed in the gap. So a segment that
# grazes land by less than the tolerance is clear, wherever and however it
# grazes.
#
# The segment is walked over every grid line it crosses; each land cell it
# enters, and the cell it starts in, is tested against its shrunk shape, and
# every closed-gap corner beside a crossing against its square. A segment
# leaves the walk at its first blocked cell.
line_clear <- function(info, x0, y0, x, y) {
  n <- length(x)
  x0 <- rep_len(x0, n)
  y0 <- rep_len(y0, n)
  # u runs east and v runs south, in cell widths from the grid's top-left
  # corner, so grid lines sit at whole numbers of u and v.
  seg <- list(u0 = (x0 - info$xmin) / info$res[1], v0 = (info$ymax - y0) / info$res[2])
  seg$du <- (x - x0) / info$res[1]
  seg$dv <- -(y - y0) / info$res[2]
  geo <- land_shape(info)
  clear <- !start_cell_blocked(geo, seg)
  clear <- crossing_walk(geo, clear, seg, "u")
  crossing_walk(geo, clear, seg, "v")
}

# A segment that passes within this distance of land, in metres, is clear: a
# thousandth of the cell side, 10 cm on 100 m cells. It is far above
# coordinate rounding (seven decimals of a degree is about 1 cm), so stations
# rounded off their cell centres keep their lines of sight, and below the
# 25 cm by which a line between two cell centres of a 100 m grid 200 cells
# across misses any corner it does not pass through.
touch_tolerance <- function(info) 1e-3 * min(info$res)

# What the land test needs, once per call: land as a vector over the grid
# with two rings of land round it, the closed-gap corners, which grid lines
# hold one, and the tolerance in cell widths along u and v.
land_shape <- function(info) {
  w <- info$ncol + 4
  pad <- matrix(TRUE, w, info$nrow + 4)
  pad[3:(info$ncol + 2), 3:(info$nrow + 2)] <- !info$water
  geo <- list(pad = as.vector(pad), w = w, ncol = info$ncol, nrow = info$nrow)
  # Corner (cu, cv), lines from 0, is element cv * (ncol + 1) + cu + 1.
  cu <- rep(0:info$ncol, times = info$nrow + 1)
  cv <- rep(0:info$nrow, each = info$ncol + 1)
  tl <- land_at(geo, cv, cu)
  tr <- land_at(geo, cv, cu + 1)
  bl <- land_at(geo, cv + 1, cu)
  br <- land_at(geo, cv + 1, cu + 1)
  gap <- (tl & br & !tr & !bl) | (tr & bl & !tl & !br)
  tol <- touch_tolerance(info)
  c(geo, list(gap = gap, any_gap = any(gap),
              gap_u = tabulate(cu[gap] + 1, info$ncol + 1) > 0,
              gap_v = tabulate(cv[gap] + 1, info$nrow + 1) > 0,
              tu = tol / info$res[1], tv = tol / info$res[2]))
}

# Land at row r and column c of the grid, from 1, for r and c from -1 to
# nrow + 2 and ncol + 2: the two rings round the grid are land.
land_at <- function(geo, r, c) geo$pad[(r + 1) * geo$w + c + 2]

# The part of each segment inside the box [ua, ub] x [va, vb], as segment
# parameters t0 to t1; empty where t0 > t1.
clip_box <- function(seg, i, ua, ub, va, vb) {
  t0 <- rep(0, length(i))
  t1 <- rep(1, length(i))
  for (ax in list(list(p = seg$u0[i], d = seg$du[i], lo = ua, hi = ub),
                  list(p = seg$v0[i], d = seg$dv[i], lo = va, hi = vb))) {
    flat <- ax$d == 0
    ta <- (ax$lo - ax$p) / ax$d
    tb <- (ax$hi - ax$p) / ax$d
    t0 <- ifelse(flat, t0, pmax(t0, pmin(ta, tb)))
    t1 <- ifelse(flat, t1, pmin(t1, pmax(ta, tb)))
    t0[flat & (ax$p < ax$lo | ax$p > ax$hi)] <- 2
  }
  list(t0 = t0, t1 = t1)
}

# TRUE where segment i touches the shrunk shape of land cell (r, c), rows and
# columns from 1. The shape is the cell less the tolerance on sides facing
# water, less a square of side the tolerance at each corner whose diagonal
# neighbour is water while both cells beside the corner are land.
cell_blocks <- function(geo, seg, i, r, c) {
  wet <- function(dr, dc) !land_at(geo, r + dr, c + dc)
  tu <- geo$tu
  tv <- geo$tv
  cut <- clip_box(seg, i, c - 1 + tu * wet(0, -1), c - tu * wet(0, 1),
                  r - 1 + tv * wet(-1, 0), r - tv * wet(1, 0))
  hit <- cut$t0 <= cut$t1
  if (!any(hit)) return(hit)
  pu <- function(t) seg$u0[i] + t * seg$du[i]
  pv <- function(t) seg$v0[i] + t * seg$dv[i]
  ends <- list(u = cbind(pu(cut$t0), pu(cut$t1)), v = cbind(pv(cut$t0), pv(cut$t1)))
  for (k in list(c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))) {
    notch <- wet(k[1], k[2]) & !wet(k[1], 0) & !wet(0, k[2])
    cu <- c - (k[2] < 0)
    cv <- r - (k[1] < 0)
    inside <- abs(ends$u - cu) <= tu & abs(ends$v - cv) <= tv
    hit <- hit & !(notch & inside[, 1] & inside[, 2])
  }
  hit
}

# TRUE where segment i touches the square of the closed-gap corner (cu, cv)
# and does not start inside it.
gap_blocks <- function(geo, seg, i, cu, cv) {
  cut <- clip_box(seg, i, cu - geo$tu, cu + geo$tu, cv - geo$tv, cv + geo$tv)
  starts_in <- abs(seg$u0[i] - cu) <= geo$tu & abs(seg$v0[i] - cv) <= geo$tv
  cut$t0 <= cut$t1 & !starts_in
}

# The cell each segment starts in, taken just past the start in its direction
# of travel, tested when it is land.
start_cell_blocked <- function(geo, seg) {
  eps <- 1e-9
  r <- floor(seg$v0 + eps * sign(seg$dv)) + 1
  c <- floor(seg$u0 + eps * sign(seg$du)) + 1
  out <- rep(FALSE, length(r))
  land <- which(land_at(geo, r, c))
  if (length(land)) out[land] <- cell_blocks(geo, seg, land, r[land], c[land])
  out
}

# One axis of line_clear(): every crossing of a grid line of `axis` strictly
# after the start and before the end. The active lines are compacted after
# each crossing.
crossing_walk <- function(geo, clear, seg, axis) {
  if (axis == "u") {
    a0 <- seg$u0
    da <- seg$du
    b0 <- seg$v0
    db <- seg$dv
    gap_line <- geo$gap_u
  } else {
    a0 <- seg$v0
    da <- seg$dv
    b0 <- seg$u0
    db <- seg$du
    gap_line <- geo$gap_v
  }
  s <- sign(da)
  first <- ifelse(s > 0, floor(a0) + 1, ceiling(a0) - 1)
  last <- ifelse(s > 0, ceiling(a0 + da) - 1, floor(a0 + da) + 1)
  count <- ifelse(s == 0, 0, pmax(0, s * (last - first) + 1))
  i <- which(clear & count >= 1)
  k <- first[i]
  si <- s[i]
  left <- count[i]
  a0i <- a0[i]
  b0i <- b0[i]
  slope <- (db / da)[i]
  # On the side of b the segment travels to, when it crosses exactly at a corner
  nudge <- 1e-9 * sign(db)[i]
  while (length(i)) {
    b <- b0i + (k - a0i) * slope
    lo <- floor(b + nudge)
    after <- k + (si > 0)
    land <- if (axis == "u") land_at(geo, lo + 1, after) else land_at(geo, after, lo + 1)
    blocked <- rep(FALSE, length(i))
    h <- which(land)
    if (length(h)) {
      blocked[h] <- if (axis == "u") cell_blocks(geo, seg, i[h], lo[h] + 1, after[h]) else
        cell_blocks(geo, seg, i[h], after[h], lo[h] + 1)
    }
    if (geo$any_gap) {
      h <- which(!blocked & gap_line[k + 1])
      if (length(h)) blocked[h] <- crossing_gaps(geo, seg, i[h], k[h], floor(b[h]), axis)
    }
    clear[i[blocked]] <- FALSE
    go <- !blocked & left > 1
    i <- i[go]
    si <- si[go]
    k <- k[go] + si
    left <- left[go] - 1
    a0i <- a0i[go]
    b0i <- b0i[go]
    slope <- slope[go]
    nudge <- nudge[go]
  }
  clear
}

# TRUE where a crossing of grid line k, at b between corners lo and lo + 1 of
# that line, touches the square of a closed-gap corner at either.
crossing_gaps <- function(geo, seg, i, k, lo, axis) {
  out <- rep(FALSE, length(i))
  for (at in list(lo, lo + 1)) {
    cu <- if (axis == "u") k else at
    cv <- if (axis == "u") at else k
    ok <- cu >= 0 & cu <= geo$ncol & cv >= 0 & cv <= geo$nrow
    h <- which(ok & !out)
    h <- h[geo$gap[cv[h] * (geo$ncol + 1) + cu[h] + 1]]
    if (length(h)) out[h] <- gap_blocks(geo, seg, i[h], cu[h], cv[h])
  }
  out
}

#' @export
print.water_distances <- function(x, ...) {
  rx <- attr(x, "receivers")
  key <- attr(x, "key")
  if (is.null(rx) || is.null(key)) return(print(as.data.frame(unclass(x)), ...))
  cat_wrapped(sprintf("<water_distances> %s: %s (by %s) x %s\n", plural(nrow(x), "pair"),
                      plural(nrow(rx), "receiver"), key,
                      plural(length(attr(x, "cells")), "water cell")))
  cat_wrapped(sprintf("Straight line used for %s of %s pairs; longest distance %s m\n",
                      fmt_n(sum(x$in_water)), fmt_n(nrow(x)),
                      fmt_n(round(if (nrow(x)) max(x$distance_m) else 0))))
  print(utils::head(as.data.frame(unclass(x), stringsAsFactors = FALSE)), ...)
  invisible(x)
}
