# ---- Independent reference for distances through water ----
# Adapted from an independent reference implementation, with the package's
# own closed-gap square. It therefore tests the code against the rule as
# detfield states it, and is not independent evidence that the rule is
# right. Land is everything outside
# the water cells, the grid border included. A segment is blocked when it
# meets the land eroded by the tolerance (GEOS, mitre joins, so square corners
# stay square), or meets the square of half-side the tolerance round a point
# where two land cells meet diagonally, unless it starts inside that square.
# Bends: every grid corner touching both land and water, except closed-gap
# points and corners along a straight shore. Shortest paths: Floyd-Warshall.
ref_build <- function(r, tol = 1e-3 * min(terra::res(r))) {
  m <- matrix(!is.na(terra::values(r, mat = FALSE)), terra::nrow(r), terra::ncol(r), byrow = TRUE)
  e <- unname(as.vector(terra::ext(r)))
  res <- terra::res(r)
  nr <- nrow(m)
  nc <- ncol(m)
  land_pad <- matrix(TRUE, nr + 2, nc + 2)
  land_pad[2:(nr + 1), 2:(nc + 1)] <- !m
  ci <- expand.grid(i = 0:nr, j = 0:nc)
  q <- cbind(land_pad[cbind(ci$i + 1, ci$j + 1)], land_pad[cbind(ci$i + 1, ci$j + 2)],
             land_pad[cbind(ci$i + 2, ci$j + 1)], land_pad[cbind(ci$i + 2, ci$j + 2)])
  nl <- rowSums(q)
  pinch <- nl == 2 & ((q[, 1] & q[, 4]) | (q[, 2] & q[, 3]))
  straight <- nl == 2 & !pinch
  xy <- function(k) cbind(e[1] + ci$j[k] * res[1], e[4] - ci$i[k] * res[2])
  polys <- terra::as.polygons(terra::ifel(is.na(r), NA, 1), dissolve = TRUE)
  water <- sf::st_union(sf::st_make_valid(sf::st_as_sf(polys)))
  water <- sf::st_set_crs(water, NA)
  box <- sf::st_as_sfc(sf::st_bbox(c(xmin = e[1] - 2 * res[1], ymin = e[3] - 2 * res[2],
                                     xmax = e[2] + 2 * res[1], ymax = e[4] + 2 * res[2])))
  land <- sf::st_buffer(sf::st_difference(box, water), -tol, joinStyle = "MITRE", mitreLimit = 10)
  gap <- xy(which(pinch))
  squares <- lapply(seq_len(nrow(gap)), function(k) {
    sf::st_polygon(list(rbind(gap[k, ] + c(-tol, -tol), gap[k, ] + c(tol, -tol),
                              gap[k, ] + c(tol, tol), gap[k, ] + c(-tol, tol),
                              gap[k, ] + c(-tol, -tol))))
  })
  list(land = land, gap = gap, squares = sf::st_sfc(squares),
       verts = xy(which(nl >= 1 & nl <= 3 & !pinch & !straight)), tol = tol)
}

# TRUE where the segment from row i of a to row i of b is clear.
ref_visible <- function(ref, a, b) {
  if (!nrow(a)) return(logical())
  seg <- sf::st_sfc(lapply(seq_len(nrow(a)), function(i) sf::st_linestring(rbind(a[i, ], b[i, ]))))
  hit <- lengths(sf::st_intersects(seg, ref$land)) > 0
  if (nrow(ref$gap)) {
    touch <- sf::st_intersects(seg, ref$squares)
    for (i in which(lengths(touch) > 0)) {
      g <- ref$gap[touch[[i]], , drop = FALSE]
      starts_in <- abs(g[, 1] - a[i, 1]) <= ref$tol & abs(g[, 2] - a[i, 2]) <= ref$tol
      if (any(!starts_in)) hit[i] <- TRUE
    }
  }
  !hit
}

# Distance through water from each row of src to each row of tgt; Inf where
# no path.
ref_distances <- function(ref, src, tgt) {
  v <- ref$verts
  nv <- nrow(v)
  dmat <- matrix(Inf, nv, nv)
  diag(dmat) <- 0
  pr <- which(upper.tri(dmat), arr.ind = TRUE)
  ok <- ref_visible(ref, v[pr[, 1], , drop = FALSE], v[pr[, 2], , drop = FALSE])
  len <- sqrt(rowSums((v[pr[, 1], , drop = FALSE] - v[pr[, 2], , drop = FALSE])^2))
  dmat[pr[ok, , drop = FALSE]] <- len[ok]
  dmat[pr[ok, 2:1, drop = FALSE]] <- len[ok]
  for (k in seq_len(nv)) dmat <- pmin(dmat, outer(dmat[, k], dmat[k, ], "+"))
  tv <- expand.grid(t = seq_len(nrow(tgt)), v = seq_len(nv))
  seen <- ref_visible(ref, tgt[tv$t, , drop = FALSE], v[tv$v, , drop = FALSE])
  leg <- matrix(Inf, nrow(tgt), nv)
  leg[cbind(tv$t, tv$v)[seen, , drop = FALSE]] <-
    sqrt(rowSums((tgt[tv$t, , drop = FALSE] - v[tv$v, , drop = FALSE])^2))[seen]
  t(vapply(seq_len(nrow(src)), function(s) {
    from <- matrix(src[s, ], nv, 2, byrow = TRUE)
    first <- rep(Inf, nv)
    sv <- ref_visible(ref, from, v)
    first[sv] <- sqrt(rowSums((v[sv, , drop = FALSE] - from[sv, , drop = FALSE])^2))
    dv <- apply(first + dmat, 2, min)
    via <- apply(sweep(leg, 2, dv, "+"), 1, min)
    direct <- ref_visible(ref, matrix(src[s, ], nrow(tgt), 2, byrow = TRUE), tgt)
    ifelse(direct, sqrt((tgt[, 1] - src[s, 1])^2 + (tgt[, 2] - src[s, 2])^2), via)
  }, numeric(nrow(tgt))))
}

# A 100 m lake of 14 x 14 cells with each kind of place a receiver can sit
# beside: a one-cell island, two land cells meeting along an edge (north
# and south of each other, and east and west), a closed gap, and a shore.
feature_lake <- function() {
  g <- terra::rast(nrows = 14, ncols = 14, xmin = 0, xmax = 1400, ymin = 0, ymax = 1400,
                   crs = "EPSG:32617", vals = 1)
  land <- rbind(c(4, 4), c(4, 8), c(4, 9), c(8, 8), c(9, 8), c(10, 3), c(11, 4),
                cbind(14, 1:14), cbind(13, 10:14), c(12, 12))
  g[terra::cellFromRowCol(g, land[, 1], land[, 2])] <- NA
  g
}

# Points on and beside the features of feature_lake(): the target point itself
# and points 1 mm, 1 cm, 11 cm and 0.5 m from it in eight directions.
placed_points <- function() {
  targets <- rbind(
    island_corner = c(300, 1100), island_edge = c(350, 1100),
    shared_edge = c(800, 1050), shared_edge_end = c(800, 1100),
    shared_edge_h = c(750, 600), closed_gap = c(300, 400),
    shore_corner = c(1100, 300), shore_edge = c(600, 100),
    open_corner = c(1100, 900), open_line = c(1150, 700),
    island_corner_se = c(400, 1000), shared_edge_h_end = c(700, 600),
    open_corner_2 = c(500, 500), open_line_v = c(1000, 1250)
  )
  dirs <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1), c(-1, 1), c(1, -1), c(-1, -1))
  dirs <- dirs / sqrt(rowSums(dirs^2))
  out <- do.call(rbind, lapply(rownames(targets), function(nm) {
    p <- targets[nm, ]
    r <- c(0, rep(c(1e-3, 1e-2, 0.11, 0.5), each = 8))
    d <- rbind(c(0, 0), dirs[rep(1:8, 4), ])
    data.frame(place = nm, x = p[1] + r * d[, 1], y = p[2] + r * d[, 2], offset_m = r)
  }))
  out$station_id <- sprintf("P%03d", seq_len(nrow(out)))
  out
}
