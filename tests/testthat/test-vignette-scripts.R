# ---- Vignette scripts ----
# Each vignette's code is also a plain R script, vignettes/<name>.R, written
# by tools/purl_vignettes.R and committed, with a copy in inst/scripts/ that
# installs. The pkgdown-only articles in vignettes/articles/ get theirs from
# tools/purl_articles.R in the same way. A script must match its vignette
# and run on the files the package installs. In a source checkout the
# scripts are read from vignettes/ and vignettes/articles/; in R CMD check,
# from the scripts folder of the installed package.

# Scripts that took more than 3 s on the build machine run only off CRAN.
slow_scripts <- c("from-detqc.R", "habitat-selection.R", "positioning.R", "range-model.R")

source_vignettes <- function() {
  testthat::test_path("..", "..", "vignettes")
}

vignette_scripts <- function() {
  src <- source_vignettes()
  dirs <- if (dir.exists(src)) {
    c(src, file.path(src, "articles"))
  } else {
    system.file("scripts", package = "detfield")
  }
  normalizePath(list.files(dirs, pattern = "\\.R$", full.names = TRUE))
}

test_that("every vignette has a script", {
  scripts <- vignette_scripts()
  expect_gt(length(scripts), 0)
  src <- source_vignettes()
  skip_if_not(dir.exists(src), "vignette sources are not beside the tests")
  rmds <- list.files(c(src, file.path(src, "articles")), pattern = "\\.Rmd$")
  expect_setequal(sub("\\.Rmd$", ".R", rmds), basename(scripts))
})

test_that("each committed script matches its vignette", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")
  src <- source_vignettes()
  tools <- test_path("..", "..", "tools", c("purl_vignettes.R", "purl_articles.R"))
  skip_if_not(dir.exists(src) && all(file.exists(tools)),
              "vignette sources are not beside the tests")
  env <- new.env()
  for (tool in tools) sys.source(tool, envir = env)
  out_dir <- withr::local_tempdir()
  article_dir <- withr::local_tempdir()
  fresh <- c(env$purl_vignettes(dirname(src), out_dir = out_dir, install_dir = NULL),
             env$purl_articles(dirname(src), out_dir = article_dir, install_dir = NULL))
  installed <- file.path(dirname(src), "inst", "scripts")
  for (f in fresh) {
    home <- if (dirname(f) == article_dir) file.path(src, "articles") else src
    for (committed in file.path(c(home, installed), basename(f))) {
      expect_true(file.exists(committed), info = committed)
      expect_identical(readLines(committed), readLines(f), info = committed)
    }
  }
})

for (script in vignette_scripts()) {
  test_that(paste("vignette script", basename(script), "runs without error or warning"), {
    if (basename(script) %in% slow_scripts) skip_on_cran()
    code <- readLines(script)
    if (any(grepl("glmmTMB::", code, fixed = TRUE))) {
      # glmmTMB warns on load when its TMB build differs from the installed
      # TMB; the warning concerns the installation, not the script.
      suppressWarnings(skip_if_not_installed("glmmTMB"))
    }
    if (any(grepl("detqc::", code, fixed = TRUE))) skip_if_not_installed("detqc")
    withr::local_dir(withr::local_tempdir())
    withr::local_pdf(NULL)
    # The scripts show the errors a reader will meet inside try(); their text
    # goes to a null connection, not the test log. Top-level values print, as
    # they do for a reader, so each figure is drawn and its warnings count.
    sink_to <- file(nullfile(), open = "w")
    withr::defer(close(sink_to))
    withr::local_options(try.outFile = sink_to)
    warned <- character()
    expect_error(
      withCallingHandlers(
        utils::capture.output(source(script, local = new.env(parent = globalenv()),
                                     print.eval = TRUE)),
        warning = function(w) {
          warned <<- c(warned, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage")
      ),
      NA
    )
    expect_identical(warned, character())
  })
}
