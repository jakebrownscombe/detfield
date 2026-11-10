# ---- Fixtures ----
utc <- function(x) as.POSIXct(x, tz = "UTC")

rt_receivers <- function() {
  data.frame(receiver_serial = c("480101", "480102"), station_id = c("RT1", "RT2"),
             deploy_datetime = utc("2024-05-20 09:30:00"),
             recover_datetime = utc("2024-05-20 13:10:00"),
             deploy_lat = c(44.5500, 44.5500), deploy_lon = c(-78.2500, -78.2400),
             deploy_z = c(3, 12))
}

rt_tags <- function() {
  data.frame(transmitter = c("A69-1602-9001", "A69-1602-9002"), tag_model = "V13-1x",
             deploy_datetime = utc("2024-05-20 09:45:00"),
             recover_datetime = utc("2024-05-20 13:05:00"),
             deploy_lat = c(44.5509, 44.5599), deploy_lon = -78.2500,
             min_delay = 90, max_delay = 150)
}

rt_detections <- function() {
  data.frame(datetime = utc(c("2024-05-20 10:02:00", "2024-05-20 10:31:00",
                              "2024-05-20 11:15:00")),
             receiver_serial = "480101", transmitter = "A69-1602-9001")
}

# ---- Every pair and every whole bin ----
test_that("every receiver, tag and whole hour is a row, with zero for silence", {
  out <- range_test_counts(rt_detections(), rt_receivers(), rt_tags())
  # Whole hours 10:00, 11:00 and 12:00; 2 receivers x 2 tags x 3 hours.
  expect_identical(nrow(out), 12L)
  expect_identical(sum(out$heard), 3L)
  expect_identical(sum(out$heard == 0), 10L)
  expect_true(all(out$sent == 30L))
  expect_identical(names(out), c("datetime", "receiver_serial", "station_id", "depth_m",
                                 "transmitter", "tag_model", "distance_m", "sent", "heard"))
  expect_identical(attr(out$datetime, "tzone"), "UTC")
  one <- out[out$receiver_serial == "480101" & out$transmitter == "A69-1602-9001", ]
  expect_identical(one$heard, c(2L, 1L, 0L))
  expect_identical(unique(out$depth_m[out$receiver_serial == "480102"]), 12)
})

test_that("distances match the deployment positions", {
  out <- range_test_counts(rt_detections(), rt_receivers(), rt_tags())
  d <- unique(out[, c("receiver_serial", "transmitter", "distance_m")])
  near <- d$distance_m[d$receiver_serial == "480101" & d$transmitter == "A69-1602-9001"]
  # 0.0009 degrees of latitude is about 100 m at 44.55 N.
  expect_equal(near, 100, tolerance = 0.01)
  pts <- sf::st_as_sf(data.frame(lon = c(-78.24, -78.25), lat = c(44.55, 44.5599)),
                      coords = c("lon", "lat"), crs = 4326)
  far <- d$distance_m[d$receiver_serial == "480102" & d$transmitter == "A69-1602-9002"]
  expect_equal(far, as.numeric(sf::st_distance(pts[1, ], pts[2, ])), tolerance = 0.005)
})

test_that("a single receiver keeps the tag-hours it did not hear", {
  # With one receiver, counting only the hours with a detection would drop
  # every silent hour and bias the curve outward.
  out <- range_test_counts(rt_detections(), rt_receivers()[1, ], rt_tags())
  expect_identical(nrow(out), 6L)
  expect_identical(sum(out$heard > 0), 2L)
})

test_that("a count above the transmissions raises sent, with a message", {
  det <- rt_detections()
  det <- det[rep(1, 31), ]
  det$datetime <- utc("2024-05-20 10:00:00") + seq(10, 3500, length.out = 31)
  expect_message(out <- range_test_counts(det, rt_receivers(), rt_tags()), "raised")
  row <- out[out$heard == 31L, ]
  expect_identical(row$sent, 31L)
})

test_that("detections that cannot be counted are reported", {
  det <- rbind(rt_detections(),
               data.frame(datetime = utc("2024-05-20 09:50:00"), receiver_serial = "480101",
                          transmitter = "A69-1602-9001"),
               data.frame(datetime = utc("2024-05-20 10:50:00"), receiver_serial = "480101",
                          transmitter = "A69-9001-1"))
  expect_message(range_test_counts(det, rt_receivers(), rt_tags()), "2 of 5 detections")
})

test_that("a deployment with no recovery ends at the last detection", {
  rx <- rt_receivers()
  rx$recover_datetime <- utc(NA)
  det <- rbind(rt_detections(), data.frame(datetime = utc("2024-05-20 12:59:00"),
                                           receiver_serial = "480101",
                                           transmitter = "A69-1602-9001"))
  out <- range_test_counts(det, rx, rt_tags())
  expect_identical(sort(unique(format(out$datetime, "%H"))), c("10", "11"))
})

test_that("inputs that cannot be counted stop with what to fix", {
  local <- rt_detections()
  local$datetime <- as.POSIXct(format(local$datetime), tz = "America/Toronto")
  expect_error(range_test_counts(local, rt_receivers(), rt_tags()), "must be in UTC")
  no_delay <- rt_tags()
  no_delay$min_delay <- NULL
  no_delay$max_delay <- NULL
  expect_error(range_test_counts(rt_detections(), rt_receivers(), no_delay), "`delay` column")
  part <- rt_tags()
  part$max_delay[2] <- NA
  expect_error(range_test_counts(rt_detections(), rt_receivers(), part), "no delay for")
  no_lat <- rt_receivers()
  no_lat$deploy_lat <- NULL
  expect_error(range_test_counts(rt_detections(), no_lat, rt_tags()), "no column `deploy_lat`")
  late <- rt_tags()
  late$deploy_datetime <- utc("2024-06-01 00:00:00")
  late$recover_datetime <- utc("2024-06-02 00:00:00")
  expect_error(range_test_counts(rt_detections(), rt_receivers(), late), "for a whole bin")
  expect_error(range_test_counts(rt_detections(), rt_receivers(), rbind(rt_tags(), rt_tags())),
               "one row per transmitter")
})
