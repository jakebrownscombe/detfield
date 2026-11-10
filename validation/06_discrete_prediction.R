# Title: Discrete and per-row prediction of a bam() detection model
# Description: The number that NEWS.md states about de_model() and a
#   bam(discrete = TRUE) fit. The model of the range-model article is fitted
#   to stoney_range_covariates and given to de_model() twice: as de_model()
#   gives it, which predicts each pair on its own (discrete = FALSE), and
#   with discrete = TRUE, which discretises all the pairs of the field
#   together. Both build the field of the article (the 70 receivers of
#   stoney_track on stoney_depth.tif, a V13-1x tag) on the article's two
#   days, and the logits of each receiver-cell pair are compared.
# Run from the package root:
#   Rscript -e 'source("validation/06_discrete_prediction.R", echo = TRUE)' \
#     > validation/outputs/06_discrete_prediction.log 2>&1
# Output: validation/outputs/06_discrete_prediction.log

# ---- Libraries ----
library(detfield)

# ---- The article's model ----
range_test <- stoney_range_covariates
first_day <- as.Date("2024-05-22")
range_test$day <- as.numeric(as.Date(range_test$datetime) - first_day)
range_test$receiver_hour <- factor(paste(range_test$receiver_serial,
                                         format(range_test$datetime)))
range_test$receiver_serial <- factor(range_test$receiver_serial)
range_test$tag_model <- factor(range_test$tag_model)
nrow(range_test)  # receiver-tag-hour rows
time_fit <- mgcv::bam(
  cbind(heard, sent - heard) ~ distance_m * depth_m + tag_model +
    s(day, k = 5) + s(receiver_serial, bs = "re") +
    s(receiver_hour, bs = "re"),
  family = binomial, data = range_test, discrete = TRUE
)
sds <- 1 / sqrt(time_fit$sp[c("s(receiver_serial)", "s(receiver_hour)")])
spread <- sqrt(sum(sds^2))
spread  # logit spread of one receiver-hour

day_of <- function(date) as.numeric(as.Date(date) - first_day)
add_levels <- function(pairs) {
  pairs$receiver_serial <- levels(range_test$receiver_serial)[1]
  pairs$receiver_hour <- levels(range_test$receiver_hour)[1]
  pairs
}
model_of <- function(...) {
  de_model(time_fit, sigma = list(station_hour = spread),
           modulators = list(day = day_of), covariates = add_levels,
           exclude = c("s(receiver_serial)", "s(receiver_hour)"), ...)
}
per_row <- model_of()
together <- model_of(discrete = TRUE)

# ---- Fields on the article's two days ----
lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
distances <- water_distances(lake, stoney_track$receivers)
days <- as.Date(c("2024-05-22", "2024-05-24"))
diffs <- lapply(days, function(date) {
  a <- suppressMessages(de_field(distances, per_row, tag_model = "V13-1x", date = date))
  b <- suppressMessages(de_field(distances, together, tag_model = "V13-1x", date = date))
  la <- stats::qlogis(a$p)
  lb <- stats::qlogis(b$p)
  ok <- is.finite(la) & is.finite(lb)
  data.frame(date = date, pairs = sum(ok), max_abs_logit = max(abs(la[ok] - lb[ok])),
             median_abs_logit = stats::median(abs(la[ok] - lb[ok])),
             pairs_p_01 = sum(ok & a$p >= 0.01),
             max_abs_logit_p_01 = max(abs(la - lb)[ok & a$p >= 0.01]))
})
res <- do.call(rbind, diffs)
print(res, digits = 3)
cat(sprintf("Largest logit difference, discrete against per-row: %.3f over %s pairs a day\n",
            max(res$max_abs_logit), format(res$pairs[1], big.mark = ",")))

# ---- Session ----
sessionInfo()
