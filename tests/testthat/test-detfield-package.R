test_that("the contract columns, data.table symbols and plot columns are declared as globals", {
  contract <- c(
    ".I", ".N", ":=", "datetime", "receiver_serial", "transmitter",
    "animal_id", "station_id", "sensor_raw", "sensor_value", "sensor_unit",
    "deploy_datetime", "recover_datetime", "deploy_lat", "deploy_lon",
    "deploy_z", "x", "y"
  )
  plot_columns <- c("p", "curve", "prob", "region", "ring", "y_count", ".count", ".key")
  expect_setequal(utils::globalVariables(package = "detfield"), c(contract, plot_columns))
})

# ---- Help pages and vignettes ----
rd_all <- function() {
  man <- system.file("man", package = "detfield")
  if (nzchar(man) && length(list.files(man, "\\.Rd$"))) {
    files <- list.files(man, "\\.Rd$", full.names = TRUE)
    return(stats::setNames(lapply(files, tools::parse_Rd), basename(files)))
  }
  tools::Rd_db("detfield")
}

test_that("no example builds a de_field by hand", {
  db <- rd_all()
  ex <- vapply(db, function(rd) {
    tmp <- tempfile()
    tools::Rd2ex(rd, tmp)
    if (file.exists(tmp)) paste(readLines(tmp), collapse = "\n") else ""
  }, character(1))
  expect_gt(sum(nzchar(ex)), 10)
  expect_false(any(grepl("class = \"de_field\"", ex, fixed = TRUE)))
})

test_that("two vignettes run the chain", {
  # The Get started page and the positioning article each run the whole
  # chain, from distances to draws.
  for (page in c("detfield.Rmd", "positioning.Rmd")) {
    vig <- c(system.file("doc", page, package = "detfield"),
             test_path("..", "..", "vignettes", page))
    vig <- vig[nzchar(vig) & file.exists(vig)]
    if (!length(vig)) next
    src <- paste(readLines(vig[1]), collapse = "\n")
    for (fn in c("range_test_counts(", "water_distances(", "de_field(", "bin_detections(",
                 "location_surface(", "credible_region(", "sample_locations(")) {
      if (fn == "range_test_counts(" && page != "detfield.Rmd") next
      expect_match(src, fn, fixed = TRUE, info = page)
    }
  }
})

# ---- DESCRIPTION ----
test_that("a suite package in Suggests is named in Remotes", {
  # detqc is not on CRAN, so the Remotes field tells remotes::install_github()
  # where to find it. No r-universe serves the suite, so the DESCRIPTION
  # names no Additional_repositories; add one with this test when it does.
  desc <- utils::packageDescription("detfield")
  suggests <- trimws(sub("\\(.*$", "", strsplit(desc$Suggests, ",")[[1]]))
  siblings <- intersect(suggests, c("detqc", "arraydesign", "arraydesignR"))
  expect_true("detqc" %in% siblings)
  expect_null(desc$Additional_repositories)
  remotes <- trimws(strsplit(if (is.null(desc$Remotes)) "" else desc$Remotes, ",")[[1]])
  expect_setequal(remotes[nzchar(remotes)], paste0("jakebrownscombe/", siblings))
})

test_that("the version floors cover linewidth and local_mocked_bindings", {
  desc <- utils::packageDescription("detfield")
  deps <- gsub("\\s+", " ", paste(desc$Imports, desc$Suggests, sep = ","))
  expect_match(deps, "ggplot2 (>= 3.5.0)", fixed = TRUE)
  expect_match(deps, "testthat (>= 3.1.7)", fixed = TRUE)
})
