test_that("counts print with the thousands separator and agree with their noun", {
  expect_identical(detfield:::fmt_n(100000), "100,000")
  expect_identical(detfield:::fmt_n(2078L), "2,078")
  expect_identical(detfield:::plural(1, "row"), "1 row")
  expect_identical(detfield:::plural(4240, "row"), "4,240 rows")
  expect_identical(detfield:::agree(1, "has", "have"), "has")
  expect_identical(detfield:::agree(2, "has", "have"), "have")
})

test_that("a listening = FALSE message prints its counts with separators", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  counts <- counts[rep(1, 1200), ]
  counts$bin_start <- counts$bin_start + 3600 * seq_len(1200)
  counts$listening <- FALSE
  counts <- rbind(counts, counts_for(field, c(3, 1, 0, 0), 12))
  expect_message(suppressWarnings(location_surface(counts, field)),
                 "1,200 rows \\(3,600 detections\\) with listening = FALSE were left out")
})
