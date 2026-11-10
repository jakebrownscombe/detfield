#!/usr/bin/env Rscript
# ---- Article scripts ----
# Usage: Rscript tools/purl_articles.R [package-dir]
# Writes vignettes/articles/<name>.R from each vignettes/articles/<name>.Rmd,
# with the purl_one() of tools/purl_vignettes.R, and copies it to
# inst/scripts/. pkgdown builds these pages; R CMD build leaves them out.
# purl_vignettes.R is shared byte for byte across the suite and reads only
# vignettes/, so this detfield file covers the articles folder. Run it after
# editing an article and commit both files;
# tests/testthat/test-vignette-scripts.R fails when they differ.

purl_articles <- function(pkg_dir = ".",
                          out_dir = file.path(pkg_dir, "vignettes", "articles"),
                          install_dir = file.path(pkg_dir, "inst", "scripts")) {
  shared <- new.env()
  sys.source(file.path(pkg_dir, "tools", "purl_vignettes.R"), envir = shared)
  rmds <- list.files(file.path(pkg_dir, "vignettes", "articles"),
                     pattern = "\\.Rmd$", full.names = TRUE)
  outs <- file.path(out_dir, sub("\\.Rmd$", ".R", basename(rmds)))
  for (i in seq_along(rmds)) {
    shared$purl_one(rmds[i], outs[i])
    # The shared header names its own tool and vignettes/; name this one.
    lines <- readLines(outs[i], warn = FALSE, encoding = "UTF-8")
    lines[1] <- sub("^# Code of the vignette ", "# Code of the article ", lines[1])
    lines[2] <- sprintf("# Extracted from vignettes/articles/%s by tools/purl_articles.R, a repository file",
                        basename(rmds[i]))
    lines[3] <- "# the installed package does not carry. The article is the source; this file follows it."
    writeLines(lines, outs[i], useBytes = TRUE)
  }
  if (!is.null(install_dir)) {
    dir.create(install_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(outs, install_dir, overwrite = TRUE)
    outs <- c(outs, file.path(install_dir, basename(outs)))
  }
  invisible(outs)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  outs <- purl_articles(if (length(args)) args[1] else ".")
  cat(sprintf("wrote %s\n", outs), sep = "")
}
