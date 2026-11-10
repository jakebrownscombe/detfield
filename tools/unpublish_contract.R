#!/usr/bin/env Rscript
# ---- Remove the maintainers' pages from the built site ----
# Usage: Rscript tools/unpublish_contract.R [package-dir]
# CONTRACT.md is the suite's shared table contract, byte-identical across the
# packages, and .github/CONTRIBUTING.md lists the development checks; both are
# written for the maintainers. pkgdown 2.0.9 renders every top-level and
# .github .md file as a page (pkgdown:::package_mds() has no exclusion and
# does not read .Rbuildignore), so run this after pkgdown::build_site(). It
# deletes docs/CONTRACT.html and docs/CONTRIBUTING.html, then rebuilds the search index and sitemap from
# the pages left.
# Last, it rewrites rdrr.io help links to the suite packages' own reference
# pages, a step that runs only in a checkout of the whole suite.

args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- normalizePath(if (length(args)) args[1] else ".")
docs <- file.path(pkg_dir, "docs")

pages <- file.path(docs, c("CONTRACT.html", "CONTRIBUTING.html"))
cat("present before:", paste(basename(pages), file.exists(pages), collapse = ", "), "\n")
unlink(pages)
# build_search() rebuilds search.json from the pages left; sitemap.xml keeps
# its entries, so the lines naming the two pages are dropped here.
pkgdown::build_search(pkg_dir)
sitemap <- file.path(docs, "sitemap.xml")
if (file.exists(sitemap)) {
  xml <- paste(readLines(sitemap, warn = FALSE), collapse = "\n")
  xml <- gsub("\\s*<url>\\s*<loc>[^<]*(CONTRACT|CONTRIBUTING)\\.html</loc>\\s*</url>", "", xml)
  xml <- gsub("\\s*<url>\\s*</url>", "", xml)
  writeLines(xml, sitemap)
}

mentions <- function(file) {
  sum(grepl("CONTRACT.html|CONTRIBUTING.html", readLines(file, warn = FALSE)))
}
cat("sitemap.xml lines naming either page:", mentions(file.path(docs, "sitemap.xml")), "\n")
cat("search.json lines naming either page:", mentions(file.path(docs, "search.json")), "\n")
html <- list.files(docs, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
linking <- html[vapply(html, mentions, numeric(1)) > 0]
cat("pages linking either page:", length(linking), "of", length(html), "\n")

# ---- Point suite help links at the suite's own reference pages ----
# downlit falls back to rdrr.io when it cannot resolve the suite's sites; a
# script of the suite checkout rewrites those links.
relink <- file.path(dirname(pkg_dir), "tools", "relink_reference.py")
status <- system2("python3", c("-I", shQuote(relink), shQuote(docs)))
if (!identical(status, 0L)) stop("relink_reference.py failed with status ", status, call. = FALSE)
