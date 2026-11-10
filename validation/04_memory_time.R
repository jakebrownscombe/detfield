# Title: Peak memory and time of water_distances()
# Description: Measures the peak memory and time that ?water_distances
#   compares under Memory, as does the comment above corner_sight() in
#   R/water_paths.R, and the times of the exact and the grid path that
#   NEWS.md gives. Every run is its own R process on the shipped example
#   (stoney_depth.tif and the 70 receivers of stoney_track, no max_distance,
#   one core). Its peak resident memory is the
#   peak of the whole process: R start-up and loading detfield and terra as
#   well as the call, so a ratio of two peaks is closer to 1 than the ratio
#   of the calls alone. Its time is the call alone, without R start-up.
#   Part 1: the lines of sight from the corners are tested in batches of
#   about 200,000 lines (the default) or 2,000,000, in ten pairs of runs;
#   each pair gives the ratio of the peaks and of the times. The order
#   within a pair alternates, so neither size always runs first.
#   Part 2: the exact path (the default) and the grid path
#   (control = list(path = "grid")) on the 100 m raster, five runs each, and
#   on the same raster split into 25 m cells by terra::disagg(), three runs
#   each. Medians are printed.
# Needs /usr/bin/time, which reports the peak with -l on macOS and -v on
#   Linux. Memory and time depend on the machine and its load; the load is
#   printed at the start and the end.
# Run from the package root:
#   Rscript -e 'source("validation/04_memory_time.R", echo = TRUE)' \
#     > validation/outputs/04_memory_time.log 2>&1
# Output: validation/outputs/04_memory_time.log

# ---- Libraries ----
library(detfield)

# ---- One run in its own process ----
# The child takes the batch size, the path and the disaggregation factor of
# the raster, and prints rows, a checksum of the distances and the elapsed
# time of the water_distances() call.
child <- tempfile(fileext = ".R")
writeLines(c(
  "a <- commandArgs(TRUE)",
  "chunk <- as.numeric(a[1])",
  "suppressPackageStartupMessages(library(detfield))",
  "f <- get('corner_sight', asNamespace('detfield'))",
  "formals(f)$chunk <- chunk",
  "utils::assignInNamespace('corner_sight', f, 'detfield')",
  "g <- terra::rast(system.file('extdata', 'stoney_depth.tif', package = 'detfield'))",
  "if (as.numeric(a[3]) > 1) g <- terra::disagg(g, as.numeric(a[3]))",
  "ctrl <- list(path = a[2])",
  "t <- system.time(d <- suppressWarnings(water_distances(g, stoney_track$receivers,",
  "                                                       control = ctrl)))",
  "cat(nrow(d), sum(d$distance_m), t[['elapsed']], '\\n')"
), child)
time_flag <- if (Sys.info()[["sysname"]] == "Darwin") "-l" else "-v"
rscript <- file.path(R.home("bin"), "Rscript")
if (nzchar(Sys.which("uptime"))) system("uptime", intern = TRUE)  # load at the start
parallel::detectCores()  # cores of this machine; each run uses one

run_once <- function(chunk = 2e5, path = "exact", fact = 1) {
  out <- system2("/usr/bin/time",
                 c(time_flag, shQuote(rscript), shQuote(child), chunk, path, fact),
                 stdout = TRUE, stderr = TRUE)
  vals <- scan(text = out[grepl("^[0-9]+ [0-9.e+]+ [0-9.]+ $", out)], quiet = TRUE)
  rss_line <- out[grepl("maximum resident set size|Maximum resident set size", out)]
  rss <- as.numeric(regmatches(rss_line, regexpr("[0-9]+", rss_line)))
  # macOS reports bytes and Linux kilobytes.
  if (time_flag == "-v") rss <- rss * 1024
  data.frame(chunk = chunk, path = path, cell_m = 100 / fact, rows = vals[1],
             checksum = vals[2], elapsed_s = vals[3], peak_mb = rss / 2^20)
}

# ---- Part 1: batch size, ten pairs ----
n_pairs <- 10
batches <- do.call(rbind, lapply(seq_len(n_pairs), function(pair) {
  sizes <- if (pair %% 2) c(2e6, 2e5) else c(2e5, 2e6)
  cbind(pair = pair, do.call(rbind, lapply(sizes, function(s) run_once(chunk = s))))
}))
print(batches, digits = 6)
# The batch size changes memory and time only, so every run gives one table.
length(unique(batches$rows)) == 1 && length(unique(batches$checksum)) == 1

big <- batches[batches$chunk == 2e6, ]
small <- batches[batches$chunk == 2e5, ]
small <- small[match(big$pair, small$pair), ]
ratio <- data.frame(pair = big$pair, memory_ratio = small$peak_mb / big$peak_mb,
                    time_ratio = small$elapsed_s / big$elapsed_s)
print(ratio, digits = 3)
nrow(ratio)  # pairs
range(ratio$memory_ratio)  # peak of 200,000 over the peak of 2,000,000
range(ratio$time_ratio)  # time of 200,000 over the time of 2,000,000
cat(sprintf("Memory ratio %.2f to %.2f, time ratio %.3f to %.3f, over %d pairs\n",
            min(ratio$memory_ratio), max(ratio$memory_ratio),
            min(ratio$time_ratio), max(ratio$time_ratio), nrow(ratio)))

# ---- Part 2: exact and grid path ----
# Runs alternate between the two paths, so load drifts over both alike.
path_runs <- function(fact, runs) {
  do.call(rbind, lapply(seq_len(runs), function(r) {
    cbind(run = r, rbind(run_once(path = "exact", fact = fact),
                         run_once(path = "grid", fact = fact)))
  }))
}
paths <- rbind(path_runs(1, 5), path_runs(4, 3))
print(paths, digits = 6)
path_med <- aggregate(cbind(elapsed_s, peak_mb, rows) ~ path + cell_m, data = paths,
                      FUN = stats::median)
path_med$runs <- aggregate(elapsed_s ~ path + cell_m, data = paths, FUN = length)$elapsed_s
print(path_med, digits = 4)  # medians per path and cell size
if (nzchar(Sys.which("uptime"))) system("uptime", intern = TRUE)  # load at the end

# ---- Session ----
sessionInfo()
