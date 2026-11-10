# nolint start: object_usage_linter.
# The fixtures come from helper-fixtures.R, which the linter cannot see.

t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")

toy_detections <- function() {
  data.frame(datetime = t0 + c(100, 400, 700, 3700), receiver_serial = c("R1", "R1", "R2", "R1"),
             transmitter = "A69-9001-1", animal_id = "a1", stringsAsFactors = FALSE)
}
toy_days <- function() {
  data.frame(receiver_serial = c("R1", "R2", "R3"), date = as.Date("2024-06-01"))
}

test_that("counts carry a zero row for every silent listening receiver", {
  out <- bin_detections(toy_detections(), toy_days(), tag_delay = 180)
  expect_s3_class(out, "data.frame")
  expect_false(inherits(out, "data.table"))
  expect_named(out, c("animal_id", "bin_start", "bin_end", "receiver_serial", "y", "n",
                      "listening"))
  expect_equal(nrow(out), 6)
  expect_equal(out$y, c(2L, 1L, 0L, 1L, 0L, 0L))
  expect_equal(out$n, rep(20L, 6))
  expect_true(all(out$listening))
  expect_equal(attr(out, "tzone"), NULL)
  expect_equal(attr(out$bin_start, "tzone"), "UTC")
  expect_equal(unique(as.numeric(out$bin_end - out$bin_start, units = "secs")), 3600)
  expect_equal(attr(out, "bin"), 3600)
})

test_that("codes on one schedule share n; independent codes add theirs", {
  det <- toy_detections()
  det$transmitter[4] <- "A69-9001-2"
  delays <- data.frame(transmitter = c("A69-9001-1", "A69-9001-2"), delay = c(90, 120))
  out <- bin_detections(det, toy_days(), delays, codes = "independent")
  expect_equal(unique(out$n), 3600 / 90 + 3600 / 120)
  out <- bin_detections(det, toy_days(), delays, codes = "shared")
  expect_equal(unique(out$n), round(3600 / 105))
  expect_equal(attr(out, "codes")$rule, "shared")
  mm <- data.frame(transmitter = c("A69-9001-1", "A69-9001-2"), min_delay = c("60", "60"),
                   max_delay = c("120", "180"), stringsAsFactors = FALSE)
  out <- bin_detections(det, toy_days(), mm, codes = "independent")
  expect_equal(unique(out$n), as.integer(round(3600 / 90 + 3600 / 120)))
  # A `codes` column sets the rule per transmitter and beats the argument.
  delays$codes <- "shared"
  out <- bin_detections(det, toy_days(), delays, codes = "independent")
  expect_equal(unique(out$n), round(3600 / 105))
})

test_that("an animal with two codes and no rule stops instead of doubling n", {
  # Two-code tags send their codes in turn (most consecutive detections
  # switch code, none within 60 s),
  # so summing the codes' rates set n = 60 where the tag sent 30.
  det <- toy_detections()
  det$transmitter[4] <- "A69-9001-2"
  expect_error(bin_detections(det, toy_days(), 120), "1 animal carries more than one transmitter")
  expect_error(bin_detections(det, toy_days(), 120, codes = "turns"), "`codes` must be NULL")
  both <- data.frame(transmitter = c("A69-9001-1", "A69-9001-2"), delay = 120,
                     codes = c("shared", "independent"))
  expect_error(bin_detections(det, toy_days(), both), "both rules")
  odd <- data.frame(transmitter = "A69-9001-1", delay = 120, codes = "turns")
  expect_error(bin_detections(det, toy_days(), odd), "holds \"turns\"")
  twice <- data.frame(transmitter = rep("A69-9001-1", 2), delay = 120,
                      codes = c("shared", "independent"))
  expect_error(bin_detections(det, toy_days(), twice), "two delays or two rules")
  # One code needs no rule.
  expect_equal(unique(bin_detections(toy_detections(), toy_days(), 120)$n), 30L)
})

test_that("the rate check flags n doubled for a tag whose codes take turns", {
  # A tag beside R1 sends one transmission every 120 s, codes alternating,
  # for 150 hours; R1 logs 95% of them.
  set.seed(4)
  ts <- seq(0, 150 * 3600 - 1, by = 120)
  heard <- stats::runif(length(ts)) < 0.95
  det <- data.frame(datetime = as.POSIXct("2024-06-01", tz = "UTC") + ts[heard],
                    receiver_serial = "R1",
                    transmitter = c("A69-1601-1", "A69-1601-2")[seq_along(ts) %% 2 + 1][heard],
                    animal_id = "bass-3")
  days <- data.frame(receiver_serial = "R1", date = as.Date("2024-06-01") + 0:7)
  delays <- data.frame(transmitter = c("A69-1601-1", "A69-1601-2"), min_delay = 80,
                       max_delay = 160)
  expect_warning(wrong <- bin_detections(det, days, delays, codes = "independent"),
                 "For 1 animal .* below 0.6 or above 1.5 times n.*bass-3")
  expect_lt(attr(wrong, "rate_check")$ratio, 0.6)
  right <- expect_silent(bin_detections(det, days, delays, codes = "shared"))
  chk <- attr(right, "rate_check")
  expect_equal(chk$bins, 150L)
  expect_equal(chk$n, 30)
  expect_gt(chk$ratio, 0.9)
  # Under 100 detected bins there is no ratio and no warning.
  short <- bin_detections(det[det$datetime < min(det$datetime) + 50 * 3600, ], days, delays,
                          codes = "independent")
  expect_true(is.na(attr(short, "rate_check")$ratio))
  # Animals with odd and even numbers of bins in one call: the median of an
  # odd count of integers is an integer, and data.table needs one type.
  two <- toy_detections()
  two$animal_id[4] <- "a2"
  chk <- attr(bin_detections(rbind(toy_detections(), two), toy_days(), 180), "rate_check")
  expect_equal(nrow(chk), 2)
})

test_that("a code listed for an animal counts even when never detected", {
  delays <- data.frame(animal_id = c("a1", "a1", "other"),
                       transmitter = c("A69-9001-1", "A69-9001-9", "A69-9001-5"),
                       delay = c(180, 180, 60))
  out <- bin_detections(toy_detections(), toy_days(), delays, codes = "independent")
  expect_equal(unique(out$n), 40L)
  out <- bin_detections(toy_detections(), toy_days(), delays, codes = "shared")
  expect_equal(unique(out$n), 20L)
  expect_equal(attr(out, "codes")$transmitters, 2L)
})

test_that("the floor raises n to the largest count and is reported", {
  det <- data.frame(datetime = t0 + seq(10, 3500, length.out = 30), receiver_serial = "R1",
                    transmitter = "A69-9001-1", animal_id = "a1")
  expect_message(out <- bin_detections(det, toy_days(), tag_delay = 600),
                 "raised to the largest receiver count in 1 of 1 bins")
  expect_equal(unique(out$n), 30L)
  expect_equal(attr(out, "n_raised"), 1L)
})

test_that("detections at receivers outside the listening table are kept and flagged", {
  det <- toy_detections()
  det$receiver_serial[3] <- "R9"
  expect_warning(out <- bin_detections(det, toy_days(), 180), "1 detection in 1 bin was logged")
  expect_equal(out$listening[out$receiver_serial == "R9"], FALSE)
  expect_equal(attr(out, "not_listening"), 1L)
})

test_that("the listening column of a receiver-day table removes days", {
  # Receivers failing the reference-tag check are marked
  # not listening in the listening table, and no row is built for them.
  days <- toy_days()
  days$listening <- c(TRUE, FALSE, TRUE)
  out <- suppressWarnings(bin_detections(toy_detections(), days, 180))
  expect_false(any(out$receiver_serial == "R2" & out$listening))
})

test_that("a deployment listens in a bin it covers for at least min_cover", {
  # A receiver deployed at 10:59 covers 1 minute of the 10:00 bin. The format
  # is given because as.POSIXct() takes one format for a whole vector and
  # would read 10:59 as midnight beside a date with no time.
  at <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M", tz = "UTC")
  dep <- data.frame(receiver_serial = c("R1", "R2", "R3"), station_id = c("S1", "S2", "S3"),
                    deploy_datetime = at(c("2024-01-01 00:00", "2024-01-01 00:00",
                                           "2024-06-01 10:59")),
                    recover_datetime = at(c(NA, "2024-06-01 10:30", NA)))
  out <- bin_detections(toy_detections(), dep, 180)
  first <- out[out$bin_start == t0, ]
  second <- out[out$bin_start == t0 + 3600, ]
  expect_setequal(first$receiver_serial, c("R1", "R2"))
  expect_setequal(second$receiver_serial, c("R1", "R3"))
  expect_equal(attr(out, "cover_dropped"), 1L)
  expect_equal(out$station_id[out$receiver_serial == "R1"][1], "S1")
  strict <- bin_detections(toy_detections()[-3, ], dep, 180, min_cover = 0.6)
  expect_setequal(strict$receiver_serial[strict$bin_start == t0], "R1")
  loose <- bin_detections(toy_detections(), dep, 180, min_cover = 0)
  expect_setequal(loose$receiver_serial[loose$bin_start == t0], c("R1", "R2", "R3"))
  expect_error(bin_detections(toy_detections(), dep, 180, min_cover = 2), "`min_cover` must")
  # A detection at a receiver that covers too little of its bin is kept as
  # not listening, with a warning.
  late <- toy_detections()
  late$receiver_serial[1] <- "R3"
  late$datetime[1] <- t0 + 3590
  expect_warning(out <- bin_detections(late, dep, 180), "in the water for less than")
  expect_false(out$listening[out$receiver_serial == "R3" & out$bin_start == t0])
})

test_that("deployments and receiver-days together drop failed days and dry hours", {
  # Receiver-days alone list every hour of a deployment day, so a
  # receiver recovered at 10:05 would count as silent until midnight.
  at <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M", tz = "UTC")
  dep <- data.frame(receiver_serial = c("R1", "R2", "R3"), station_id = c("S1", "S2", "S3"),
                    deploy_datetime = at(rep("2024-01-01 00:00", 3)),
                    recover_datetime = at(c(NA, NA, "2024-06-01 10:05")))
  days <- data.frame(receiver_serial = c("R1", "R2", "R3"), date = as.Date("2024-06-01"),
                     listening = c(TRUE, FALSE, TRUE))
  only_days <- suppressWarnings(bin_detections(toy_detections(), days, 180))
  expect_true("R3" %in% only_days$receiver_serial[only_days$bin_start == t0])
  both <- suppressWarnings(bin_detections(toy_detections(), list(deployments = dep, days = days),
                                          180))
  expect_setequal(both$receiver_serial[both$listening], "R1")
  expect_false(any(both$receiver_serial == "R3"))
  expect_equal(both$station_id[both$receiver_serial == "R1"][1], "S1")
  expect_error(bin_detections(toy_detections(), list(dep, days), 180),
               "must hold `deployments` and `days`")
  expect_error(bin_detections(toy_detections(), list(deployments = dep, days = days), 180,
                              bin = 7 * 3600), "must divide a day")
})

test_that("a receiver swap at one station keeps one row with both receivers' detections", {
  # A swap day gives duplicate rows per station unless they are merged.
  dep <- data.frame(receiver_serial = c("A", "B", "C"), station_id = c("S1", "S1", "S2"),
                    deploy_datetime = as.POSIXct(c("2024-01-01 00:00", "2024-06-01 10:30",
                                                   "2024-01-01 00:00"), tz = "UTC"),
                    recover_datetime = as.POSIXct(c("2024-06-01 10:20", NA, NA), tz = "UTC"))
  # A covers a third of the bin and B half: the station covers 0.83 and
  # listens, although neither receiver alone covers more than half.
  det <- data.frame(datetime = t0 + c(100, 2000, 2100), receiver_serial = c("A", "B", "B"),
                    transmitter = "A69-1", animal_id = "a1")
  out <- bin_detections(det, dep, 60)
  expect_equal(out$receiver_serial, c("B", "C"))
  expect_equal(out$y, c(3L, 0L))
  expect_equal(attr(out, "swaps_merged"), 1L)
  days <- data.frame(receiver_serial = c("A", "A", "B", "C"),
                     station_id = c("S1", "S1", "S1", "S2"),
                     date = as.Date(c("2024-05-31", "2024-06-01", "2024-06-01", "2024-06-01")))
  out2 <- bin_detections(det, days, 60)
  expect_equal(out2$receiver_serial, c("B", "C"))
  expect_equal(out2$y, c(3L, 0L))
})

test_that("bins are aligned to UTC and take any length", {
  out <- bin_detections(toy_detections(), toy_days(), 180, bin = "day")
  expect_equal(unique(format(out$bin_start, tz = "UTC")), "2024-06-01")
  expect_equal(unique(out$n), 480L)
  out <- bin_detections(toy_detections(), toy_days(), 180, bin = 1800)
  expect_equal(as.numeric(unique(out$bin_start) - t0, units = "secs"), c(0, 3600))
  expect_error(bin_detections(toy_detections(), toy_days(), 180, bin = 7 * 3600),
               "must divide a day")
  expect_error(bin_detections(toy_detections(), toy_days(), 180, bin = "week"), "`bin` must be")
})

test_that("each bin's counts do not depend on the rest of the record", {
  # No 10-day windows are needed, because bins are
  # counted one at a time.
  full <- bin_detections(toy_detections(), toy_days(), 180)
  part <- bin_detections(toy_detections()[4, ], toy_days(), 180)
  expect_equal(part$y, full$y[full$bin_start == t0 + 3600])
})

test_that("rows with no animal_id are left out with a message", {
  det <- toy_detections()
  det$animal_id[1] <- NA
  expect_message(out <- bin_detections(det, toy_days(), 180), "1 of 4 detections")
  expect_equal(sum(out$y), 3)
  det$animal_id <- NA_character_
  expect_error(suppressMessages(bin_detections(det, toy_days(), 180)), "no detection")
})

test_that("bad inputs stop with the column or value named", {
  det <- toy_detections()
  expect_error(bin_detections(det[, -2], toy_days(), 180), "no column `receiver_serial`")
  bad <- det
  bad$receiver_serial <- factor(bad$receiver_serial)
  expect_error(bin_detections(bad, toy_days(), 180), "must be character, not factor")
  bad <- det
  attr(bad$datetime, "tzone") <- "America/Toronto"
  expect_error(bin_detections(bad, toy_days(), 180), "must be in UTC")
  expect_error(bin_detections(det, toy_days(), -1), "positive number")
  expect_error(bin_detections(det, toy_days(), data.frame(transmitter = "x", delay = 1)),
               "no delay for 1 transmitter: A69-9001-1")
  expect_error(bin_detections(det, toy_days(), data.frame(transmitter = rep("A69-9001-1", 2),
                                                          delay = c(1, 2))), "two delays")
  expect_error(bin_detections(det, toy_days(), data.frame(transmitter = "A69-9001-1")),
               "needs a `delay` column")
  expect_error(bin_detections(det, toy_days(), data.frame(transmitter = "A69-9001-1",
                                                          delay = -5)), "must be positive")
  expect_error(bin_detections(det, toy_days(), list(1)), "data frame with a `transmitter`")
  bad <- det
  bad$datetime[1] <- NA
  expect_error(bin_detections(bad, toy_days(), 180), "must have no NA")
  expect_error(bin_detections(det, "x", 180), "must be a data frame")
})

test_that("a keyed data.table gives the same counts as a data frame", {
  # A stale data.table key can make a join return partial matches.
  det <- data.table::as.data.table(toy_detections())
  data.table::setkey(det, receiver_serial)
  days <- data.table::as.data.table(toy_days())
  data.table::setkey(days, receiver_serial, date)
  expect_equal(bin_detections(det, days, 180), bin_detections(toy_detections(), toy_days(), 180))
})

# ---- Deployments without serials ----
otn_like <- function() {
  # A deployment sheet with stations but no serials, as glatos's OTN
  # example; keyed by serial alone, none of its detections would count.
  deps <- data.frame(receiver_serial = NA_character_, station_id = c("HFX1", "HFX2", "HFX3"),
                     deploy_datetime = t0 - 86400, recover_datetime = t0 + 86400,
                     stringsAsFactors = FALSE)
  det <- toy_detections()
  det$receiver_serial <- c("VR2W-1", "VR2W-1", "VR2W-2", "VR2W-1")
  det$station_id <- c("HFX1", "HFX1", "HFX2", "HFX1")
  list(deps = deps, det = det)
}

test_that("deployments without serials are matched to detections by station", {
  x <- otn_like()
  expect_message(out <- bin_detections(x$det, x$deps, tag_delay = 180),
                 "3 stations have listening rows with no receiver_serial; 4 detections")
  expect_true(all(out$listening))
  expect_equal(nrow(out), 6L)
  expect_equal(out$receiver_serial[out$bin_start == t0], paste0("station:HFX", 1:3))
  expect_equal(out$y, c(2L, 1L, 0L, 1L, 0L, 0L))
  expect_equal(attr(out, "not_listening"), 0L)
  both <- list(deployments = x$deps,
               days = data.frame(receiver_serial = NA_character_, station_id = "HFX1",
                                 date = as.Date("2024-06-01")))
  expect_warning(out <- suppressMessages(bin_detections(x$det, both, tag_delay = 180)),
                 "1 detection in 1 bin was logged by receivers not listening")
  expect_equal(unique(out$receiver_serial[out$listening]), "station:HFX1")
  # A days table already keyed by station needs no change.
  keyed <- list(deployments = x$deps,
                days = data.frame(receiver_serial = paste0("station:HFX", 1:3),
                                  date = as.Date("2024-06-01")))
  out <- suppressMessages(bin_detections(x$det, keyed, tag_delay = 180))
  expect_true(all(out$listening))
  x$det$station_id <- NULL
  w <- capture_warnings(bin_detections(x$det, x$deps, tag_delay = 180))
  expect_match(w, "`detections` has no station_id column", all = FALSE)
})

test_that("tag_delay takes the shortest and longest delay as two numbers", {
  expect_equal(bin_detections(toy_detections(), toy_days(), tag_delay = c(60, 180))$n,
               bin_detections(toy_detections(), toy_days(), tag_delay = 120)$n)
  expect_error(bin_detections(toy_detections(), toy_days(), tag_delay = c(180, 60)),
               "shortest and longest")
  expect_error(bin_detections(toy_detections(), toy_days(), tag_delay = 1:3),
               "two\\s+numbers\\s+\\(shortest\\s+and\\s+longest\\s+delay\\)")
})

test_that("empty listening tables stop with the cause named, and no base warning leaks", {
  deps <- data.frame(receiver_serial = character(), station_id = character(),
                     deploy_datetime = as.POSIXct(character(), tz = "UTC"),
                     recover_datetime = as.POSIXct(character(), tz = "UTC"))
  expect_error(bin_detections(toy_detections(), deps, tag_delay = 180),
               "`listening` has no rows")
  expect_error(bin_detections(toy_detections(), toy_days()[0, ], tag_delay = 180),
               "`listening` has no rows")
  early <- data.frame(receiver_serial = "R1", station_id = "S1", deploy_datetime = t0 - 86400,
                      recover_datetime = t0 - 3600)
  # Without the check, base R warns "no non-missing arguments to max".
  w <- character()
  withCallingHandlers(bin_detections(toy_detections(), early, tag_delay = 180),
                      warning = function(cnd) {
                        w <<- c(w, conditionMessage(cnd))
                        invokeRestart("muffleWarning")
                      })
  expect_length(grep("not listening", w), 1)
  expect_length(grep("non-missing", w), 0)
})

test_that("the help defines the delay and leaves the V16AT rule open", {
  src <- system.file("man", "bin_detections.Rd", package = "detfield")
  rd <- if (nzchar(src)) tools::parse_Rd(src) else tools::Rd_db("detfield")[["bin_detections.Rd"]]
  txt <- gsub("\\s+", " ", paste(utils::capture.output(tools::Rd2txt(rd)), collapse = " "))
  expect_match(txt, "V16AT tags the evidence is mixed and the rule is unresolved", fixed = TRUE)
  expect_match(txt, "start of one transmission to the start of the next", fixed = TRUE)
})

# nolint end
