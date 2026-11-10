#' Five simulated animals in Stoney Lake with known tracks
#'
#' The example study of the positionR packages: Stoney Lake, Ontario, on a
#' 100 m depth raster, 70 receivers on a grid of 700 m steps with some of
#' 600 m (nearest neighbour 600 to 922 m, median 700 m), and
#' five animals walked for 48 hours from 00:00 UTC on 1 June 2024 with the
#' correlated random walk of the arraydesign simulator (30 m a minute on
#' average, step SD 15 m, turns of SD 40 degrees), each from a random water
#' cell. Each tag transmits with a delay drawn uniformly from 120 to 180 s.
#' Every transmission was detected or missed by each receiver with the
#' probability of the Stoney example's range model: the logistic curve of
#' [de_logistic()] with `d50` 250 m and `d05` 500 m at a receiver depth of
#' 1 m and 750 m and 1,500 m at 35 m, read at the distance through water,
#' with no variation between receivers or hours. A detection field built from
#' that model with a spread of 0 is therefore the model the detections came
#' from. The array is dense enough that 17% of the 212 detected animal-hours
#' are heard by a single receiver.
#'
#' @format A list of six elements:
#' \describe{
#'   \item{`receivers`}{The receiver deployments, one row per receiver: `receiver_serial`,
#'     `station_id`, `deploy_datetime`, `recover_datetime`, `deploy_lat`, `deploy_lon`,
#'     `deploy_z` (receiver depth, 1.1 to 16.1 m, the depth of its cell), and the
#'     projected `x` and `y` in UTM zone 17N (EPSG:32617), with `attr(, "epsg")`.}
#'   \item{`detections`}{3,641 detections: `datetime`, `receiver_serial`,
#'     `transmitter` and `animal_id`; 803 of `stoney-01`, 1,063 of `stoney-02`,
#'     487 of `stoney-03`, 273 of `stoney-04` and 1,015 of `stoney-05`.}
#'   \item{`tags`}{One row per transmitter: `transmitter`, `animal_id`,
#'     `min_delay` (120 s) and `max_delay` (180 s).}
#'   \item{`track`}{The true position each minute: `animal_id`, `datetime`,
#'     `x` and `y`, 2,881 rows per animal from 00:00 UTC on 1 June 2024 to
#'     00:00 UTC on 3 June.}
#'   \item{`transmissions`}{Every transmission, heard or not: `animal_id`,
#'     `transmitter`, `datetime`, `x`, `y`, `cell` (the raster cell) and
#'     `n_detections` (receivers that heard it); 1,148 to 1,153 per animal.}
#'   \item{`truth`}{The values the study was simulated with: the range model
#'     (`d50`, `d05`, `depth`, `station_hour`), the walk (`step_mean`,
#'     `step_sd`, `turn_sd`, `time_step`), `tag_delay` and the seeds.}
#' }
#' The depth raster is included as `stoney_depth.tif`, and the study as the files
#' `stoney_deployments.csv`, `stoney_detections.csv` and `stoney_tags.csv`, in
#' `system.file("extdata", package = "detfield")`. The raster stores depth as
#' metres of zero or below, with land `NA`; take `abs()` before depth
#' arithmetic or a depth map. 56 of its 4,818 water cells hold 0.
#' `vignette("detfield")` reads the files, and `vignette("positioning")`
#' scores the positions against the tracks.
#'
#' @source Simulated with the arraydesign simulator by
#'   <https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_track.R>,
#'   which is in the repository and not in the installed package. The raster
#'   has the same grid and depths as `depth_raster` of positionR 1.7.0, a
#'   100 m depth raster of Stoney Lake, derived from the bathymetry
#'   made before 2024, with the CRS set to EPSG:32617. The arraydesignR app carries the same file.
#' @seealso [stoney_range_test] for the range test of the same range model.
#' @examples
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Beyond 2 km the Stoney curve is below 0.001 at every receiver; leaving
#' # those pairs out keeps the example quick.
#' dist <- water_distances(grid, stoney_track$receivers, max_distance = 2000)
#' truth <- stoney_track$truth
#' model <- de_model(de_logistic(truth$d50, truth$d05, depth = truth$depth),
#'                   sigma = list(station_hour = 0))
#' field <- de_field(dist, model)
#' # Position one animal, to keep the example short
#' one <- stoney_track$detections[stoney_track$detections$animal_id == "stoney-01", ]
#' counts <- bin_detections(one, stoney_track$receivers, stoney_track$tags)
#' surface <- location_surface(counts, field, overdispersion = "none")
#' head(attr(credible_region(surface), "summary"))
"stoney_track"

#' A simulated range test of the Stoney Lake range model
#'
#' Hourly counts from a range test on two lines, one in 1 m of water and one
#' in 25 m, the depths that bracket the receivers of [stoney_track]. Each line
#' held three receivers and five range-test tags, so its 15 pairs span 100 to
#' 1,600 m, for 48 hours: the shallow line from 20 May 2024 and the deep line
#' from 23 May. The detections were simulated with the Stoney example's range model
#' of [stoney_track], with no variation between receivers or hours, so a
#' logistic curve fitted at each depth reproduces the example's curve at that
#' depth. The detections, the receivers and the tags of the test are included as files
#' (`stoney_range_test_detections.csv`, `stoney_range_test_receivers.csv` and
#' `stoney_range_test_tags.csv` in `system.file("extdata", package =
#' "detfield")`), and `stoney_range_test` is what [range_test_counts()] makes
#' of them.
#'
#' @format A data frame with one row per receiver, tag and hour, 1,440 rows:
#' \describe{
#'   \item{`datetime`}{Start of the hour, POSIXct in UTC.}
#'   \item{`receiver_serial`, `station_id`}{The receiver and its station.}
#'   \item{`depth_m`}{Depth of the receiver, 1 or 25 m.}
#'   \item{`transmitter`, `tag_model`}{The range-test tag and its model.}
#'   \item{`distance_m`}{Distance from the receiver to the tag, in metres.}
#'   \item{`sent`}{Transmissions the tag sent in the hour.}
#'   \item{`heard`}{Detections of the tag the receiver logged in the hour.}
#' }
#' @source Simulated by
#'   <https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_range_test.R>,
#'   which is in the repository and not in the installed package.
#' @seealso [stoney_range_covariates] for a range test with tag models, wind
#'   and variation between receivers and hours.
#' @examples
#' aggregate(cbind(heard, sent) ~ depth_m + distance_m, data = stoney_range_test,
#'           FUN = sum)
"stoney_range_test"

#' A simulated range test with two tag models, two depths and daily wind
#'
#' Hourly counts from a range test with six receivers, three at 1 m and three
#' at 25 m, and eight range-test tags of two models, over six days. For a
#' V13-1x tag at a wind of 4 m/s the curve is the Stoney example's range model of
#' [stoney_track] at each depth. A V9-2x tag is 0.9 logits lower, each m/s of
#' the day's mean wind above 4 lowers the logit by 0.25, and the counts carry
#' an effect for each receiver (standard deviation 0.5) and one for each
#' receiver and hour (standard deviation 0.9). These effects, absent from
#' [stoney_track], give a range model fitted to these counts a spread to
#' estimate. Use it to fit a range model with covariates; the counts do not
#' match the detections of [stoney_track].
#'
#' @format A data frame with one row per receiver, tag and hour, 6,912 rows,
#'   from 00:00 UTC on 22 May 2024:
#' \describe{
#'   \item{`datetime`}{Start of the hour, POSIXct in UTC.}
#'   \item{`receiver_serial`, `station_id`}{The receiver and its station.}
#'   \item{`depth_m`}{Depth of the receiver, 1 or 25 m.}
#'   \item{`transmitter`, `tag_model`}{The range-test tag and its model,
#'     `"V9-2x"` or `"V13-1x"`.}
#'   \item{`distance_m`}{Distance from the receiver to the tag, in metres.}
#'   \item{`wind_ms`}{Mean wind speed of the UTC day, in metres per second.}
#'   \item{`sent`}{Transmissions the tag sent in the hour.}
#'   \item{`heard`}{Detections of the tag the receiver logged in the hour.}
#' }
#' @source Simulated by
#'   <https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_range_test.R>,
#'   which is in the repository and not in the installed package.
#' @seealso [stoney_range_test].
#' @examples
#' aggregate(cbind(heard, sent) ~ tag_model + depth_m,
#'           data = stoney_range_covariates, FUN = sum)
"stoney_range_covariates"

#' A simulated habitat-selection study on Stoney Lake
#'
#' Forty animals tracked for 48 hours from 00:00 UTC on 1 June 2024 on the 70
#' receivers of [stoney_track], with a known preference for bottom depth. The
#' animals walk the correlated random walk of [stoney_track] (steps of mean
#' 30 m and SD 15 m a minute, turns of SD 40 degrees), and each step is
#' accepted in proportion to exp(beta x bottom depth in m), so beta is a
#' selection coefficient per step. `selection-01` to `selection-20` have beta
#' = -0.3 per m, a preference for shallow water; `selection-21` to
#' `selection-40` have beta = 0, no preference. Both groups start from the
#' same 20 water cells. Tags send every 120 to 180 s, and detections are drawn
#' from the range model of [stoney_track] with no variation between receivers
#' or hours.
#'
#' The study is part of a trial of 20 animals at each of five betas (-0.3,
#' -0.1, 0, 0.1 and 0.3) on two arrays, the 70 receivers of [stoney_track]
#' and 20 receivers on a regular grid (the array of arraydesign's examples),
#' in which the depth coefficient from position draws was set
#' against the same coefficient from the true track. The trial was run on two
#' disjoint blocks of animals; block 2 has its own starts, tracks and
#' detections. Within a block the
#' five betas share start cells. The forty animals here are
#' block 1's at beta -0.3 and 0 on 70 receivers, so the analysis of the
#' article `vignette("habitat-selection")` reproduces their rows of `trial`.
#'
#' @format A list:
#' \describe{
#'   \item{`receivers`}{The 70 receivers of [stoney_track].}
#'   \item{`detections`}{27,442 detections in the 48 hours: `datetime`,
#'     `receiver_serial`, `transmitter`, `animal_id`.}
#'   \item{`tags`}{One row per animal, 40 rows: `transmitter`, `animal_id`,
#'     `beta`, `replicate` (the animal's number, 1 to 20 within each beta; the
#'     same number marks the same start cell), `min_delay`, `max_delay`.}
#'   \item{`track`}{The true position each minute, 115,240 rows (2,881 per
#'     animal): `animal_id`, `datetime`, `cell`, the cell of the lake raster
#'     the animal was in, and `x` and `y`, the simulated position in metres
#'     (UTM 17N), kept to within 10 m of it, toward the centre of its cell.}
#'   \item{`presence`}{The position draws of the trial, 35,700 rows: 20 cells
#'     per detected hour (1,785 hours) drawn with [sample_locations()] from the
#'     hourly surfaces of the field the detections were drawn from: the draws
#'     the trial analysed for these animals. Columns and attributes as
#'     [sample_locations()] returns them.}
#'   \item{`coverage`}{A list of two [system_de()] tables, `date`, `cell` and
#'     `p`, for the three UTC days of the study and the 4,818 water cells:
#'     `transmission`, the probability that one transmission is detected, and
#'     `hour`, the probability that an hour of 24 transmissions is.}
#'   \item{`truth`}{The settings: layer, betas, hours, movement, tag delay,
#'     transmissions per hour (24), position draws and available points per
#'     stratum (20 each), and the seeds of both blocks, for reproducing them.}
#'   \item{`trial`}{3,144 rows, one per block, animal, beta, array,
#'     availability and source, for detected hours: `block` (1 or 2); `beta`;
#'     `receivers` (70 or 20); `availability`, drawn from the whole detection
#'     surface: "weighted" (in proportion to the hourly detection probability),
#'     "per_transmission" (in proportion to the per-transmission probability),
#'     "uniform_0.1" or "uniform_0.2" (uniform over cells heard with that
#'     probability or more per transmission); `source` ("truth", every third
#'     minute of the true track, or "positions", the position draws); `animal`
#'     (the replicate); `estimate` and `se`, the conditional logit coefficient
#'     for bottom depth per m pooled over 20 imputations with Rubin's rules;
#'     `n_hours`, the animal's detected hours; `single_share`, the share of
#'     them heard by one receiver; and `warned`, `TRUE` when an imputation fit
#'     warned (separation). Animals never heard have no rows: 4 of 100 on 20
#'     receivers in block 1 and 3 in block 2.}
#'   \item{`hours`}{8,872 rows, one per detected hour of the trial's animals
#'     on 70 receivers, both blocks: `block`, `beta`, `animal`, `bin_start`,
#'     `receivers` (receivers that heard the animal in the hour), `true_m`
#'     and `drawn_m`, the mean bottom depth in m of the hour's 20 true points
#'     and of its 20 position draws.}
#' }
#' @source Simulated with the walk, tag and detection simulators of the
#'   **arraydesign** R package on the Stoney raster and the range model of
#'   [stoney_track] (the raster's source is under [stoney_track]), then
#'   positioned and analysed with this package's own functions. The script
#'   that builds it is
#'   <https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_selection.R>, and
#'   <https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_selection_track_xy.R>
#'   adds the positions of `track`. Both are in the repository and not in the
#'   installed package. The scripts rebuild the study and the draws; they copy
#'   `trial` and `hours` from the outputs of the trial, which was run outside
#'   the package, and those outputs are not in the repository.
#' @seealso [sample_available()], [used_available()], [stoney_track]
#' @examples
#' table(stoney_selection$tags$beta)
#' # Population mean of the per-animal coefficients in block 1 of the trial
#' block1 <- stoney_selection$trial[stoney_selection$trial$block == 1 &
#'                                    stoney_selection$trial$availability ==
#'                                      "weighted", ]
#' aggregate(estimate ~ receivers + source + beta, data = block1, FUN = mean)[1:10, ]
#' # Position draws minus true position, bottom depth (m), per detected hour
#' hours <- stoney_selection$hours
#' tapply(hours$drawn_m - hours$true_m, hours$block, mean)
"stoney_selection"
