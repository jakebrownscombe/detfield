# ---- The shipped draws and coverage, and the article's analysis ----
# The 40 animals of stoney_selection are the trial's at beta -0.3 and 0 on 70
# receivers (block 1). Their shipped position draws and coverage must be what
# the functions give from the study, and the article's analysis of them must
# give their rows of stoney_selection$trial.
test_that("the shipped draws and coverage are what the functions give", {
  skip_on_cran()
  lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
  study <- stoney_selection
  curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
  field <- suppressMessages(de_field(water_distances(lake, study$receivers),
                                     de_model(curve, sigma = list(station_hour = 0))))
  counts <- bin_detections(study$detections, study$receivers, study$tags)
  surface <- location_surface(counts, field)
  presence <- do.call(rbind, lapply(seq_len(nrow(study$tags)), function(k) {
    hours <- surface[surface$animal_id == study$tags$animal_id[k], ]
    sample_locations(hours, n = 20, seed = study$tags$replicate[k])
  }))
  rownames(presence) <- NULL
  expect_equal(as.data.frame(presence), as.data.frame(study$presence), ignore_attr = TRUE)
  listening <- expand.grid(station_id = study$receivers$station_id,
                           date = as.Date("2024-06-01") + 0:2, stringsAsFactors = FALSE)
  expect_equal(system_de(field, listening, transmissions = 24)[c("date", "cell", "p")],
               study$coverage$hour, tolerance = 1e-12)
  expect_equal(system_de(field, listening)[c("date", "cell", "p")], study$coverage$transmission,
               tolerance = 1e-12)
})

test_that("the article's default analysis matches the trial rows", {
  skip_on_cran()
  skip_if_not_installed("survival")
  # clogit() reads a bare strata() term only with survival attached
  withr::local_package("survival")
  lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
  depth <- abs(lake)
  names(depth) <- "depth_m"
  study <- stoney_selection
  hourly <- study$coverage$hour
  pool <- function(used, seed) {
    available <- sample_available(used, hourly, lake, n = 20, seed = seed)
    points <- used_available(used, available, depth)
    fits <- sapply(1:20, function(m) {
      fit <- clogit(case ~ depth_m + strata(stratum), data = points[points$draw == m, ])
      c(stats::coef(fit), stats::vcov(fit))
    })
    c(mean(fits[1, ]), sqrt(mean(fits[2, ]) + (1 + 1 / 20) * stats::var(fits[1, ])))
  }
  true_points <- function(id, hours) {
    m <- study$track[study$track$animal_id == id, ]
    m$bin_start <- as.POSIXct(trunc(m$datetime, "hours"))
    minute <- as.numeric(format(m$datetime, "%M"))
    keep <- m$bin_start %in% hours & minute %% 3 == 0
    xy <- terra::xyFromCell(lake, m$cell[keep])
    data.frame(animal_id = id, bin_start = m$bin_start[keep], draw = minute[keep] / 3 + 1,
               cell = m$cell[keep], x = xy[, 1], y = xy[, 2])
  }
  live <- study$trial[study$trial$block == 1 & study$trial$receivers == 70 &
                        study$trial$availability == "weighted" &
                        study$trial$beta %in% c(-0.3, 0), ]
  for (k in seq_len(nrow(study$tags))) {
    tag <- study$tags[k, ]
    positions <- study$presence[study$presence$animal_id == tag$animal_id, ]
    truth <- true_points(tag$animal_id, unique(positions$bin_start))
    for (src in c("positions", "truth")) {
      got <- pool(if (src == "positions") positions else truth, 100 + tag$replicate)
      want <- live[live$beta == tag$beta & live$animal == tag$replicate & live$source == src, ]
      expect_equal(got, c(want$estimate, want$se), tolerance = 1e-8,
                   info = paste(tag$animal_id, src))
    }
  }
})

test_that("the shipped hourly depths are those of the shipped draws and track", {
  lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
  depth_of_cell <- abs(terra::values(lake, mat = FALSE))
  study <- stoney_selection
  ex <- study$hours[study$hours$block == 1 & study$hours$beta %in% c(-0.3, 0), ]
  ids <- study$tags$animal_id[match(paste(ex$beta, ex$animal),
                                    paste(study$tags$beta, study$tags$replicate))]
  drawn <- tapply(depth_of_cell[study$presence$cell],
                  paste(study$presence$animal_id, as.numeric(study$presence$bin_start)), mean)
  expect_equal(as.vector(drawn[paste(ids, as.numeric(ex$bin_start))]), ex$drawn_m,
               tolerance = 1e-12)
  minutes <- study$track
  minutes$bin_start <- as.numeric(trunc(minutes$datetime, "hours"))
  minutes <- minutes[as.numeric(format(minutes$datetime, "%M")) %% 3 == 0, ]
  true_m <- tapply(depth_of_cell[minutes$cell], paste(minutes$animal_id, minutes$bin_start), mean)
  expect_equal(as.vector(true_m[paste(ids, as.numeric(ex$bin_start))]), ex$true_m,
               tolerance = 1e-12)
})
