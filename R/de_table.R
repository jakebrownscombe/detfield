# ---- Curves from a table of values ----

table_curve <- function(tab) {
  has_tag <- "tag_model" %in% names(tab)
  has_depth <- "depth_m" %in% names(tab)
  tab <- check_curve_table(tab)
  predict <- function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
    if (!is.numeric(distance_m)) stop("`distance_m` must be numeric.", call. = FALSE)
    n <- length(distance_m)
    tm <- if (has_tag) rep_len(as.character(tag_model), n) else rep("", n)
    dp <- if (has_depth) rep_len(as.numeric(depth_m), n) else rep(0, n)
    known <- unique(tab$tag_model)
    if (has_tag && any(!tm %in% known)) {
      stop(sprintf("Unknown tag_model %s. The table holds: %s.",
                   paste(unique(tm[!tm %in% known]), collapse = ", "),
                   paste(known, collapse = ", ")), call. = FALSE)
    }
    if (has_depth && anyNA(dp)) stop("This curve needs `depth_m`.", call. = FALSE)
    eta <- numeric(n)
    combos <- unique(data.frame(tm = tm, dp = dp, stringsAsFactors = FALSE))
    for (i in seq_len(nrow(combos))) {
      at <- which(tm == combos$tm[i] & dp == combos$dp[i])
      eta[at] <- table_eta(tab[tab$tag_model == combos$tm[i], ], distance_m[at], combos$dp[i])
    }
    check_eta(eta, n)
  }
  uses <- c("distance_m", if (has_depth) "depth_m", if (has_tag) "tag_model")
  list(predict = predict, uses = uses, tag_models = if (has_tag) unique(tab$tag_model),
       source = sprintf("table of %s", plural(nrow(tab), "curve value")))
}

check_curve_table <- function(tab) {
  tab <- as.data.frame(tab)
  if (!"distance_m" %in% names(tab)) stop("The curve table needs `distance_m`.", call. = FALSE)
  if (!"eta" %in% names(tab)) {
    if (!"p" %in% names(tab)) stop("The curve table needs `eta` or `p`.", call. = FALSE)
    if (!is.numeric(tab$p) || any(!(tab$p > 0 & tab$p < 1))) {
      stop("`p` must lie strictly between 0 and 1; give `eta` for the end values.", call. = FALSE)
    }
    tab$eta <- stats::qlogis(tab$p)
  }
  check_numeric(tab, c("distance_m", "eta"), "the curve table")
  if ("depth_m" %in% names(tab)) check_numeric(tab, "depth_m", "the curve table")
  tab$tag_model <- if ("tag_model" %in% names(tab)) as.character(tab$tag_model) else ""
  key <- paste(tab$tag_model, if ("depth_m" %in% names(tab)) tab$depth_m else "")
  sizes <- table(key)
  if (any(sizes < 2) || anyDuplicated(paste(key, tab$distance_m))) {
    stop("Each curve in the table needs at least two rows with distinct distances.", call. = FALSE)
  }
  if (!"depth_m" %in% names(tab)) tab$depth_m <- 0
  tab[order(tab$tag_model, tab$depth_m, tab$distance_m), ]
}

# Linear in depth between the two nearest tabulated depths, held outside them.
table_eta <- function(tab, distance, depth) {
  depths <- sort(unique(tab$depth_m))
  if (length(depths) == 1) return(one_curve_eta(tab, distance))
  depth <- min(max(depth, depths[1]), depths[length(depths)])
  j <- findInterval(depth, depths, rightmost.closed = TRUE)
  w <- (depth - depths[j]) / (depths[j + 1] - depths[j])
  lower <- if (w < 1) one_curve_eta(tab[tab$depth_m == depths[j], ], distance) else 0
  upper <- if (w > 0) one_curve_eta(tab[tab$depth_m == depths[j + 1], ], distance) else 0
  (1 - w) * lower + w * upper
}

one_curve_eta <- function(curve, distance) {
  d <- curve$distance_m
  e <- curve$eta
  out <- stats::approx(d, e, xout = distance, rule = 2, ties = "ordered")$y
  m <- length(d)
  slope <- (e[m] - e[m - 1]) / (d[m] - d[m - 1])
  far <- distance > d[m]
  if (slope < 0 && any(far)) out[far] <- e[m] + slope * (distance[far] - d[m])
  out
}
