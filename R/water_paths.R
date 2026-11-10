# ---- Shortest paths through water ----
# The shortest path between two points of water that cannot see each other
# bends only at corners where the water turns round land: grid corners with
# exactly one land cell of the four round them (a cell off the grid counts
# as land). water_paths() finds those corners and which of them see each
# other, once per grid. For one receiver, water_path_from() runs Dijkstra
# from the receiver over the corners, then gives each cell the shortest
# distance to a corner it sees plus the straight leg from that corner. Line
# of sight is line_clear() throughout, so the paths obey the same land rule
# as the in_water flag, closed gaps included.

# The corners, the legs between corners that see each other, the water cells,
# their 4-connected water bodies, and a store for which cells each corner
# sees, filled in when a corner is first needed.
water_paths <- function(info, max_distance = Inf) {
  cn <- shore_corners(info)
  n <- nrow(cn)
  legs <- list(to = vector("list", n), len = vector("list", n))
  if (n > 1) {
    pairs <- corner_pairs(n)
    long <- (cn[pairs[, 1], 1] - cn[pairs[, 2], 1])^2 + (cn[pairs[, 1], 2] - cn[pairs[, 2], 2])^2
    pairs <- pairs[long <= max_distance^2, , drop = FALSE]
    seen <- line_clear(info, cn[pairs[, 1], 1], cn[pairs[, 1], 2], cn[pairs[, 2], 1],
                       cn[pairs[, 2], 2])
    pairs <- pairs[seen, , drop = FALSE]
    len <- sqrt((cn[pairs[, 1], 1] - cn[pairs[, 2], 1])^2 +
                  (cn[pairs[, 1], 2] - cn[pairs[, 2], 2])^2)
    from <- factor(c(pairs[, 1], pairs[, 2]), levels = seq_len(n))
    legs$to <- split(c(pairs[, 2], pairs[, 1]), from)
    legs$len <- split(c(len, len), from)
  }
  cells <- which(info$water)
  sight <- new.env(parent = emptyenv())
  sight$cells <- vector("list", n)
  list(corners = cn, legs = legs, cells = cells, xy = xy_of_cell(info, cells),
       body = water_bodies(info)[cells], sight = sight, max_distance = max_distance)
}

# Grid corners with exactly one land cell among the four round them: x, y and
# the quadrant (qx, qy) the land cell lies in, east and north positive.
shore_corners <- function(info) {
  land <- matrix(!info$water, info$nrow, info$ncol, byrow = TRUE)
  pad <- matrix(TRUE, info$nrow + 2, info$ncol + 2)
  pad[2:(info$nrow + 1), 2:(info$ncol + 1)] <- land
  i <- rep(0:info$nrow, times = info$ncol + 1)
  j <- rep(0:info$ncol, each = info$nrow + 1)
  tl <- pad[cbind(i + 1, j + 1)]
  tr <- pad[cbind(i + 1, j + 2)]
  bl <- pad[cbind(i + 2, j + 1)]
  br <- pad[cbind(i + 2, j + 2)]
  k <- tl + tr + bl + br == 1
  cbind(x = info$xmin + j[k] * info$res[1], y = info$ymax - i[k] * info$res[2],
        qx = ifelse(tr[k] | br[k], 1, -1), qy = ifelse(tl[k] | tr[k], 1, -1))
}

# Every unordered pair of 1..n, as a two-column matrix.
corner_pairs <- function(n) {
  a <- rep(seq_len(n - 1), times = (n - 1):1)
  b <- unlist(lapply(2:n, seq.int, to = n), use.names = FALSE)
  cbind(a, b)
}

# Water bodies joined through cell edges only: two water cells that meet at a
# corner and nowhere else are a closed gap and do not join.
water_bodies <- function(info) {
  r <- terra::rast(nrows = info$nrow, ncols = info$ncol, xmin = 0, xmax = info$ncol,
                   ymin = 0, ymax = info$nrow)
  r <- terra::setValues(r, ifelse(info$water, 1, NA_real_))
  terra::values(terra::patches(r, directions = 4), mat = FALSE)
}

# The cells (rows of wp$cells) that each corner of `vs` sees, with the leg to
# each, made on first use and kept in wp$sight. Only cells a shortest path can
# reach through a corner are tested: a path bends round the corner's land
# cell, so it leaves the corner into one of the two quadrants beside the
# land, never into the land's quadrant or the one opposite it. Each call of
# this function is one pass. A pass tests the lines of all its corners
# together, vectorised, in batches of about `chunk` lines. On the Stoney
# example with no max_distance and one core, where all lines are tested in
# one pass, batches of 200,000 lines give a lower peak memory than batches
# ten times larger, at about the same time;
# <https://github.com/jakebrownscombe/detfield/blob/main/validation/04_memory_time.R>
# measures both. A pass of 200,000 lines or fewer makes the same batches at
# either size.
corner_sight <- function(wp, info, vs, chunk = 2e5) {
  todo <- vs[vapply(wp$sight$cells[vs], is.null, logical(1))]
  cn <- wp$corners
  while (length(todo)) {
    batch <- list()
    n_lines <- 0
    while (length(todo) && (n_lines == 0 || n_lines < chunk)) {
      v <- todo[1]
      todo <- todo[-1]
      dx <- wp$xy[, 1] - cn[v, "x"]
      dy <- wp$xy[, 2] - cn[v, "y"]
      ex <- dx * cn[v, "qx"]
      ey <- dy * cn[v, "qy"]
      leg <- sqrt(dx^2 + dy^2)
      k <- which(!(ex > 0 & ey > 0) & !(ex < 0 & ey < 0) & leg <= wp$max_distance)
      batch[[length(batch) + 1]] <- list(v = v, k = k, leg = leg[k])
      n_lines <- n_lines + length(k)
    }
    v_of <- rep(vapply(batch, `[[`, numeric(1), "v"), vapply(batch, function(b) length(b$k), 1L))
    k_all <- unlist(lapply(batch, `[[`, "k"), use.names = FALSE)
    ok <- line_clear(info, cn[v_of, "x"], cn[v_of, "y"], wp$xy[k_all, 1], wp$xy[k_all, 2])
    ok <- split(ok, factor(v_of, levels = vapply(batch, `[[`, numeric(1), "v")))
    for (b in seq_along(batch)) {
      seen <- ok[[b]]
      wp$sight$cells[[batch[[b]]$v]] <- list(k = batch[[b]]$k[seen], leg = batch[[b]]$leg[seen])
    }
  }
  invisible(wp$sight$cells[vs])
}

# Shortest distance through water from (x0, y0) to every corner: Dijkstra
# over the legs, starting from the corners the point sees.
corner_distances <- function(wp, info, x0, y0) {
  cn <- wp$corners
  n <- nrow(cn)
  d <- rep(Inf, n)
  if (!n) return(d)
  sees <- line_clear(info, x0, y0, cn[, "x"], cn[, "y"])
  d[sees] <- sqrt((cn[sees, "x"] - x0)^2 + (cn[sees, "y"] - y0)^2)
  d[d > wp$max_distance] <- Inf
  done <- rep(FALSE, n)
  repeat {
    open <- which(!done & is.finite(d))
    if (!length(open)) break
    u <- open[which.min(d[open])]
    done[u] <- TRUE
    to <- wp$legs$to[[u]]
    if (length(to)) d[to] <- pmin(d[to], d[u] + wp$legs$len[[u]])
  }
  d
}

# Shortest distance through water from (x0, y0) to every water cell by way of
# at least one corner: the cell's best corner plus the straight leg from it.
# Inf where no corner path within wp$max_distance reaches the cell.
water_path_from <- function(wp, info, x0, y0) {
  dv <- corner_distances(wp, info, x0, y0)
  best <- rep(Inf, length(wp$cells))
  vs <- which(dv < wp$max_distance)
  sights <- corner_sight(wp, info, vs)
  for (i in seq_along(vs)) {
    s <- sights[[i]]
    if (length(s$k)) best[s$k] <- pmin(best[s$k], dv[vs[i]] + s$leg)
  }
  best
}

# With several cores, the sight of every corner within reach of a receiver is
# made up front, split across the cores, so the workers that then take one
# receiver each find it made.
fill_sight <- function(wp, info, rx, cores) {
  cn <- wp$corners
  if (cores <= 1 || !nrow(cn)) return(invisible(wp))
  reach <- vapply(seq_len(nrow(cn)), function(v) {
    any((rx$x - cn[v, "x"])^2 + (rx$y - cn[v, "y"])^2 <= wp$max_distance^2)
  }, logical(1))
  vs <- which(reach)
  if (!length(vs)) return(invisible(wp))
  parts <- split(vs, cut(seq_along(vs), min(cores, length(vs)), labels = FALSE))
  made <- parallel_apply(parts, function(p) corner_sight(wp, info, p), cores)
  for (i in seq_along(parts)) wp$sight$cells[parts[[i]]] <- made[[i]]
  invisible(wp)
}
