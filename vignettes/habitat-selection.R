# Code of the vignette "Habitat selection: depth use from hourly positions".
# Extracted from vignettes/habitat-selection.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Depth shift ----
# Attach detfield, survival for the conditional logit, and ggplot2
library(detfield)
library(survival)
library(ggplot2)
# Hour by hour on 70 receivers: the mean bottom depth of the 20 position
# draws minus that of the 20 true points (every third minute)
trial_hours <- stoney_selection$hours
trial_hours$shift_m <- trial_hours$drawn_m - trial_hours$true_m
trial_hours$band <- cut(trial_hours$true_m, c(-Inf, 1, 2.5, 5, 10, Inf),
                        labels = c("1 m or less", "1 to 2.5 m", "2.5 to 5 m",
                                   "5 to 10 m", "over 10 m"))
shift_by_band <- aggregate(shift_m ~ block + band, data = trial_hours,
                           FUN = function(m) c(mean = mean(m), hours = length(m)))
shift_by_band <- do.call(data.frame, shift_by_band)
names(shift_by_band) <- c("block", "band", "mean_m", "hours")
shift_by_band$mean_m <- round(shift_by_band$mean_m, 2)
shift_by_band
# Hours beyond 6 m either way, which the figure's axis cuts off
sum(abs(trial_hours$shift_m) > 6)

# ---- Plot depth shift ----
# Draw the hourly differences by true depth band, one panel per block
ggplot(trial_hours, aes(band, shift_m)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, fill = "grey92",
               colour = "grey35") +
  stat_summary(aes(colour = "Mean"), fun = mean, geom = "point", size = 2.2) +
  scale_colour_manual(name = NULL, values = c(Mean = "#2C5AA0")) +
  facet_wrap(~ block, labeller = as_labeller(c(`1` = "Block 1",
                                               `2` = "Block 2"))) +
  coord_cartesian(ylim = c(-6, 6)) +
  labs(x = "Bottom depth at the true position",
       y = "Draws minus true position (m)") +
  theme_positionr() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- Load ----
# The parts of the included study, and the animals by preference
study <- stoney_selection
names(study)
table(beta = study$tags$beta)

# ---- Presence ----
# One row per animal, detected hour and draw
presence <- study$presence
head(presence, 3)
writeLines(strwrap(attr(presence, "estimand"), 72))

# ---- Depth used ----
# Bottom depth in metres from the lake raster
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "detfield"))
depth <- abs(lake)
names(depth) <- "depth_m"
# The true cells at minutes 0, 3, ..., 57 of an animal's detected hours
true_points <- function(id) {
  minutes <- study$track[study$track$animal_id == id, ]
  minutes$bin_start <- as.POSIXct(trunc(minutes$datetime, "hours"))
  minute <- as.numeric(format(minutes$datetime, "%M"))
  keep <- minutes$bin_start %in% presence$bin_start[presence$animal_id == id] &
    minute %% 3 == 0
  xy <- terra::xyFromCell(lake, minutes$cell[keep])
  data.frame(animal_id = id, bin_start = minutes$bin_start[keep],
             draw = minute[keep] / 3 + 1, cell = minutes$cell[keep],
             x = xy[, 1], y = xy[, 2])
}
truth <- do.call(rbind, lapply(study$tags$animal_id, true_points))
# Depth of every true point and every draw, with the animal's preference
depth_of_cell <- terra::values(depth, mat = FALSE)
depth_used <- rbind(
  data.frame(animal_id = truth$animal_id, set = "True positions",
             depth_m = depth_of_cell[truth$cell]),
  data.frame(animal_id = presence$animal_id, set = "Position draws",
             depth_m = depth_of_cell[presence$cell])
)
depth_used$beta <- study$tags$beta[match(depth_used$animal_id,
                                         study$tags$animal_id)]
table(depth_used$set)

# ---- Plot depth used ----
# Order the animals by their mean true depth, one panel per preference
true_mean <- tapply(depth_used$depth_m[depth_used$set == "True positions"],
                    depth_used$animal_id[depth_used$set == "True positions"],
                    mean)
depth_used$animal_id <- factor(depth_used$animal_id,
                               names(sort(true_mean)))
depth_used$set <- factor(depth_used$set, c("True positions",
                                           "Position draws"))
depth_used$group <- ifelse(depth_used$beta == 0, "No preference (beta 0)",
                           "Shallow preference (beta -0.3)")
ggplot(depth_used, aes(depth_m, animal_id, fill = set)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3,
               position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ group, scales = "free_y") +
  scale_fill_manual(values = c("#2ECC40", "#F28E2B"), name = NULL) +
  scale_x_continuous(trans = "sqrt", breaks = c(0, 1, 2, 5, 10, 20, 30)) +
  labs(x = "Bottom depth (m, square-root scale)",
       y = "Animal, ordered by the mean depth of its true positions") +
  theme_positionr() +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom")

# ---- Deeper count ----
# Mean depth per animal: true positions against position draws
by_animal <- aggregate(depth_m ~ animal_id + beta + set, data = depth_used,
                       FUN = mean)
by_animal <- reshape(by_animal, idvar = c("animal_id", "beta"),
                     timevar = "set", direction = "wide")
names(by_animal)[names(by_animal) == "depth_m.True positions"] <- "true_m"
names(by_animal)[names(by_animal) == "depth_m.Position draws"] <- "drawn_m"
c(animals = nrow(by_animal),
  draws_deeper = sum(by_animal$drawn_m > by_animal$true_m))
aggregate(cbind(true_m, drawn_m) ~ beta, data = by_animal,
          FUN = function(m) round(mean(m), 2))
# Share of points in water 1 m deep or less, by preference and set
shallow_share <- aggregate(depth_m ~ beta + set, data = depth_used,
                           FUN = function(m) round(mean(m <= 1), 3))
shallow_share

# ---- Coverage ----
# Coverage per cell and day: per hour of 24 transmissions, and per
# transmission
coverage <- study$coverage
head(coverage$hour, 3)
max(abs(coverage$hour$p - (1 - (1 - coverage$transmission$p)^24)))
# Share of cell-days from which an hour is detected with p >= 0.95
mean(coverage$hour$p >= 0.95)

# ---- Available ----
# Draw 20 available cells per presence draw, with the animal's second seed
available <- lapply(split(study$tags, study$tags$animal_id), function(tag) {
  sample_available(presence[presence$animal_id == tag$animal_id, ],
                   coverage$hour, lake, n = 20, seed = 100 + tag$replicate)
})
head(available[["selection-01"]], 3)
writeLines(strwrap(attr(available[["selection-01"]], "estimand"), 72))

# ---- Available uniform ----
# Uniform over cells heard with p >= 0.1 per transmission, one animal
uniform <- sample_available(presence[presence$animal_id == "selection-01", ],
                            coverage$transmission, lake, n = 20,
                            uniform_above = 0.1, seed = 101)
writeLines(strwrap(attr(uniform, "estimand"), 72))
# On 1 June: the water cells, those this mode draws from, and those 1.x drew from
june1 <- coverage$transmission[coverage$transmission$date ==
                                 as.Date("2024-06-01"), ]
june1_cells <- c(water = nrow(june1), at_or_above = sum(june1$p >= 0.1),
                 strictly_above = sum(june1$p > 0.1))
june1_cells

# ---- Map data ----
# One animal with a preference for shallow water: its true track, its
# presence draws and the available cells of its first draw
animal <- "selection-05"
animal_track <- study$track[study$track$animal_id == animal, ]
animal_presence <- presence[presence$animal_id == animal, ]
first_draw <- available[[animal]][available[[animal]]$draw == 1, ]
# Water as a table of cells, depth in metres
depth_cells <- as.data.frame(depth, xy = TRUE, na.rm = TRUE)
c(hours = length(unique(animal_presence$bin_start)),
  presence = nrow(animal_presence), available = nrow(first_draw))

# ---- Plot map lake ----
# The fill scale holds depth, so the key takes the available cells' fill
# from override.aes
point_key <- list(shape = c(21, 19), fill = c("#E7298A", NA), size = c(2, 1.6))
# Map depth, the available cells of draw 1, then the presence draws
ggplot(depth_cells, aes(x, y)) +
  geom_raster(aes(fill = depth_m)) +
  scale_fill_depth() +
  geom_point(data = first_draw, aes(colour = "Available cells, draw 1"),
             shape = 21, fill = "#E7298A", stroke = 0.2, size = 1.2) +
  geom_point(data = animal_presence, aes(colour = "Presence draws"), size = 0.9) +
  scale_colour_manual(name = NULL,
                      values = c("Available cells, draw 1" = "black",
                                 "Presence draws" = "#F28E2B"),
                      guide = guide_legend(order = 2, override.aes = point_key)) +
  coord_equal(expand = FALSE) +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr() +
  theme(legend.position = "bottom")

# ---- Plot map animal ----
# Count the animal's detections per receiver
heard <- table(study$detections$receiver_serial[study$detections$animal_id ==
                                                  animal])
animal_rx <- study$receivers
animal_rx$detections <- as.integer(heard[animal_rx$receiver_serial])
animal_rx$detections[is.na(animal_rx$detections)] <- 0L
# Zoom to the animal's track, with 600 m round it
view <- list(x = range(animal_track$x) + c(-600, 600),
             y = range(animal_track$y) + c(-600, 600))
in_view <- animal_rx$x >= view$x[1] & animal_rx$x <= view$x[2] &
  animal_rx$y >= view$y[1] & animal_rx$y <= view$y[2]
# Map depth, the presence draws, the true track over them, the receivers,
# then the track's start
ggplot(depth_cells, aes(x, y)) +
  geom_raster(aes(fill = depth_m)) +
  scale_fill_depth() +
  geom_point(data = animal_presence, aes(colour = "Presence draws"), size = 1.3,
             alpha = 0.6, stroke = 0) +
  geom_path(data = animal_track, aes(colour = "True track"), linewidth = 0.3,
            alpha = 0.85) +
  layer_receivers(animal_rx[in_view, ], count = "detections") +
  geom_point(data = animal_track[1, ], aes(colour = "Start"), shape = 21,
             fill = "white", size = 2.6, stroke = 0.8) +
  # A constant mapped to colour gives each mark an entry in the key
  scale_colour_manual(name = NULL,
                      values = c("True track" = "grey10", "Start" = "black",
                                 "Presence draws" = "#F28E2B"),
                      breaks = c("True track", "Start", "Presence draws"),
                      guide = guide_legend(order = 3)) +
  coord_equal(xlim = view$x, ylim = view$y, expand = FALSE) +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr()

# ---- Stack ----
# Stack used and available points with bottom depth
stacked <- lapply(names(available), function(id) {
  used_available(presence[presence$animal_id == id, ], available[[id]],
                 depth)
})
names(stacked) <- names(available)
head(stacked[["selection-01"]], 3)
table(case = stacked[["selection-01"]]$case)[c("1", "0")]

# ---- Fit ----
# Fit each draw, then pool: mean estimate; within plus between variance
pool <- function(points) {
  fits <- sapply(1:20, function(m) {
    fit <- clogit(case ~ depth_m + strata(stratum),
                  data = points[points$draw == m, ])
    c(coef(fit), vcov(fit))
  })
  c(estimate = mean(fits[1, ]),
    se = sqrt(mean(fits[2, ]) + (1 + 1 / 20) * var(fits[1, ])))
}
from_positions <- t(sapply(stacked, pool))
from_truth <- t(sapply(names(available), function(id) {
  tag <- study$tags[study$tags$animal_id == id, ]
  used <- truth[truth$animal_id == id, ]
  pool(used_available(used, sample_available(used, coverage$hour, lake,
                                             n = 20,
                                             seed = 100 + tag$replicate),
                      depth))
}))
head(round(cbind(truth = from_truth[, "estimate"],
                 positions = from_positions[, "estimate"]), 3), 3)

# ---- Plot animals ----
# Coefficient per animal: true track against position draws
tags <- study$tags[match(rownames(from_truth), study$tags$animal_id), ]
animals <- data.frame(tags[c("animal_id", "beta")],
                      truth = from_truth[, "estimate"],
                      truth_se = from_truth[, "se"],
                      positions = from_positions[, "estimate"],
                      positions_se = from_positions[, "se"])
animals$group <- ifelse(animals$beta == 0, "beta 0", "beta -0.3")
ggplot(animals, aes(truth, positions, shape = group)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40") +
  geom_errorbar(aes(ymin = positions - 1.96 * positions_se,
                    ymax = positions + 1.96 * positions_se),
                colour = "#2C5AA0", alpha = 0.35, width = 0) +
  geom_errorbar(aes(xmin = truth - 1.96 * truth_se,
                    xmax = truth + 1.96 * truth_se),
                colour = "grey45", alpha = 0.35, width = 0,
                orientation = "y") +
  geom_point(colour = "#2C5AA0", fill = "#2C5AA0", size = 2) +
  scale_shape_manual(values = c(`beta -0.3` = 21, `beta 0` = 1),
                     name = NULL) +
  labs(x = "Coefficient from the true track (per m)",
       y = "Coefficient from position draws (per m)") +
  theme_positionr() +
  theme(legend.position = "bottom")

# ---- Trial table ----
# Population mean of per-animal coefficients with a t interval over animals
population <- function(rows) {
  half <- qt(0.975, nrow(rows) - 1) * sd(rows$estimate) / sqrt(nrow(rows))
  c(mean = mean(rows$estimate), lower = mean(rows$estimate) - half,
    upper = mean(rows$estimate) + half, animals = nrow(rows))
}
# Ratio of mean positions to mean truth, with a paired bootstrap over animals
ratio_ci <- function(rows) {
  both <- merge(rows[rows$source == "truth", c("animal", "estimate")],
                rows[rows$source == "positions", c("animal", "estimate")],
                by = "animal", suffixes = c("_truth", "_positions"))
  boot <- withr::with_seed(4, replicate(2000, {
    k <- sample.int(nrow(both), replace = TRUE)
    mean(both$estimate_positions[k]) / mean(both$estimate_truth[k])
  }))
  c(ratio = mean(both$estimate_positions) / mean(both$estimate_truth),
    quantile(boot, c(0.025, 0.975), names = FALSE))
}
trial <- study$trial
scenarios <- unique(trial[c("block", "receivers", "availability", "beta")])
results <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(k) {
  s <- scenarios[k, ]
  rows <- merge(trial, s)
  data.frame(s,
             t(setNames(population(rows[rows$source == "truth", ]),
                        c("truth", "truth_lower", "truth_upper", "animals"))),
             t(setNames(population(rows[rows$source == "positions", ])[1:3],
                        c("positions", "positions_lower", "positions_upper"))),
             t(setNames(ratio_ci(rows), c("ratio", "ratio_lower",
                                          "ratio_upper"))))
}))
# A ratio is shown where the true track's interval excludes zero
detected <- results$truth_lower > 0 | results$truth_upper < 0
results[!detected | results$beta == 0,
        c("ratio", "ratio_lower", "ratio_upper")] <- NA
nrow(results)
