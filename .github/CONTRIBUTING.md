# Development

`tools/check_all.sh` runs roxygen, `R CMD build`, `R CMD check --as-cran`, lintr and the code
gates of `tools/check_gates.R`: at most 80 lines per function, at most 6 formals per exported
function, a tarball of at most 2 MB, examples under 5 s, no exported name that masks base R, and
UK spelling. To lint against an installed copy of the package, run it as
`R_LIBS=/path/to/library tools/check_all.sh`.
