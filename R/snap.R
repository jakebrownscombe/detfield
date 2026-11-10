# ---- Place receivers on water cells ----
# A receiver close to shore can fall in a land cell of a coarse grid. Within
# one cell diagonal of water that is a rasterisation effect, so it moves to
# the nearest water cell; farther away the coordinates or the CRS are more
# likely wrong. The limit is measured to the edge of the nearest water cell,
# so every point of a land cell that borders water qualifies. The move goes
# to that cell's centre, so it can reach the limit plus half a diagonal.
snap_receivers <- function(rx, info, key) {
  cell <- cell_of_xy(info, rx$x, rx$y)
  outside <- is.na(cell)
  if (any(outside)) {
    stop(sprintf("%s %s outside `grid`: %s.", plural(sum(outside), "receiver"),
                 agree(sum(outside), "lies", "lie"),
                 paste(utils::head(rx[[key]][outside], 10), collapse = ", ")), call. = FALSE)
  }
  rx$snap_m <- 0
  on_land <- which(!info$water[cell])
  if (length(on_land)) {
    water <- which(info$water)
    wxy <- xy_of_cell(info, water)
    limit <- sqrt(sum(info$res^2))
    for (i in on_land) {
      to_edge <- sqrt(pmax(abs(wxy[, "x"] - rx$x[i]) - info$res[1] / 2, 0)^2 +
                        pmax(abs(wxy[, "y"] - rx$y[i]) - info$res[2] / 2, 0)^2)
      if (min(to_edge) > limit) {
        stop(sprintf("Receiver %s is %s m from the nearest water cell; %s",
                     rx[[key]][i], fmt_n(round(min(to_edge))), "check its coordinates and CRS."),
             call. = FALSE)
      }
      d <- sqrt((wxy[, "x"] - rx$x[i])^2 + (wxy[, "y"] - rx$y[i])^2)
      j <- which.min(d)
      cell[i] <- water[j]
      rx$snap_m[i] <- d[j]
      rx$x[i] <- wxy[j, "x"]
      rx$y[i] <- wxy[j, "y"]
    }
    moves <- sprintf("%s (%s m)", rx[[key]][on_land], fmt_n(round(rx$snap_m[on_land])))
    warning(sprintf(paste("%s sat on %s and moved to the centre of the nearest water cell,",
                          "by the distance in brackets (snap_m): %s."),
                    plural(length(on_land), "receiver"),
                    agree(length(on_land), "a land cell", "land cells"),
                    paste(moves, collapse = ", ")), call. = FALSE)
  }
  rx$cell <- cell
  rx
}
