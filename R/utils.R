# ---- Classed errors ----
# Raises an error of classes detfield_error_<kind> and detfield_error, so a
# caller can catch one kind of failure by its class.
# Kinds in use: "degenerate" (every bin of a location_surface() run is
# degenerate), "degenerate_bin" (one bin is, and the run stops), "worker"
# (a parallel worker returned nothing), "grid" (a grid that terra cannot
# read, or that has no CRS, a CRS in longitude and latitude, a CRS in
# units other than metres, or a CRS whose scale at the raster centre is more
# than 1 percent from 1 or cannot be measured). `call` is the user's call,
# where the caller passes it.
detfield_abort <- function(message, kind, call = NULL) {
  stop(errorCondition(message, class = c(paste0("detfield_error_", kind), "detfield_error"),
                      call = call))
}

# ---- Counts in messages ----
# A count with the thousands separator and never in scientific notation:
# "2,078", "100,000".
fmt_n <- function(k) format(k, big.mark = ",", scientific = FALSE, trim = TRUE)

# "1 animal", "2 animals", "1,690 cells"; a value of exactly 1 takes the
# singular.
plural <- function(k, word) {
  sprintf("%s %s%s", fmt_n(k), word, if (isTRUE(k == 1)) "" else "s")
}

# A message picks its verb by the count, so that "1 row was" and "2 rows
# were" both read right.
agree <- function(k, one, many) if (isTRUE(k == 1)) one else many

# ---- Parallel apply ----
# lapply() over forked processes when `cores` is above 1. Windows has no
# forks, so there it runs on one process with a warning. A worker's error is
# raised again with its own class; a worker that returns NULL, as one killed
# for memory does, stops the run.
parallel_apply <- function(x, fun, cores) {
  if (cores > 1 && .Platform$OS.type == "windows") {
    warning("control$cores > 1 needs forked processes, which Windows lacks; using 1.",
            call. = FALSE)
    cores <- 1L
  }
  if (cores <= 1 || length(x) <= 1) return(lapply(x, fun))
  # mclapply() warns when a worker fails; the failure itself is raised below.
  res <- withCallingHandlers(
    parallel::mclapply(x, fun, mc.cores = min(cores, length(x)), mc.preschedule = FALSE),
    warning = function(w) {
      if (grepl("errors in user code|resulted in an error", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  failed <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(failed)) {
    first <- res[[which(failed)[1]]]
    if (inherits(first, "try-error")) stop(attr(first, "condition"))
    detfield_abort("A worker returned nothing, possibly because it ran out of memory.", "worker")
  }
  res
}
