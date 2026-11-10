# detfield 0.1.1

First release.

## Range tests and detection models

- `range_test_counts()` counts a range test from its files (detections, receivers, and
  range-test tags with positions and delays) into hourly `sent` and `heard` per receiver and
  tag. Hours in which a receiver heard nothing from a tag are kept with `heard = 0`, because a
  curve fitted to the heard hours alone reaches too far.
- `de_logistic()` builds a detection curve from d50 and d05, the distances at which half and
  5% of transmissions are heard, at one depth or interpolated between depths. It refuses
  `d05 <= d50`. Its print gives the distances to four significant digits.
- `de_model()` builds a detection model from a fitted `glm` or `glmmTMB` model, a table of
  curve values, a function or a `de_logistic()` curve. The spread of one receiver's hourly
  detection probability about the curve, `sigma$station_hour`, is carried as given and never
  fitted. A non-`NA` `sigma$station_day` stops. Covariates that change by day, such as wind,
  enter as modulators.
- `de_model()` warns when it is given a glmmTMB fit whose formula holds a smooth (`s()`,
  `te()`, `t2()` or `ti()`) and no `re.form`: glmmTMB fits the wiggly part of a smooth as a
  random effect, so with the default `re.form = NA` the smooth reaches the field as a straight
  line. The warning names the term and the two ways out: fit the smooth in mgcv, or pass
  `re.form = NULL`.
- `de_model()` predicts a `bam(discrete = TRUE)` fit with `discrete = FALSE`, so each pair is
  predicted on its own. A discrete prediction discretises all the rows it is given together,
  which moves a pair's detection probability with the other pairs of the field: by up to 0.12
  logit on the range-model article's model
  (<https://github.com/jakebrownscombe/detfield/blob/main/validation/06_discrete_prediction.R>).
- `de_model(blocked = )` sets a logit offset for receiver-cell pairs whose line of sight
  crosses land. The default, 0, reads the curve at the distance through water.

## Detection fields

- `water_distances()` measures the distance through water from every receiver to every water
  cell, going round land, and records whether the straight line stays in water. A receiver
  just on land snaps to the centre of the nearest water cell, by at most one and a half cell
  diagonals, and the warning gives each move in metres. A grid with no land cells, no CRS or
  no water, and a serial found at two positions, stop or warn with a plain message.
- `water_distances()` reads the grid's units as metres on the ground. It refuses a grid in
  longitude and latitude, in a unit other than the metre, or in a CRS whose scale at the
  raster centre differs from 1 by more than 1 percent on its two axes or its diagonals. The
  scale is measured on four 1 km segments, and the refusal names it and the UTM zone of the
  centre. The Stoney example grid in Web Mercator has a scale of about 1.40 there and is
  refused. A CRS whose scale cannot be measured, such as a local engineering CRS, is refused
  too.
- The land test of `water_distances()` decides a line that passes within a thousandth of a
  cell side of a cell corner at the corner, whichever side it passes: grazing one land cell's
  corner leaves the line clear, and passing where two land cells meet diagonally blocks it.
  A station a few millimetres off its cell centre, as a station file with latitude and
  longitude to seven decimals gives, has the lines of sight of one at the centre.
- Where land lies on the straight line, `water_distances()` gives the exact shortest path
  through water: straight legs that bend only at the corners where water turns round land,
  found by Dijkstra's algorithm over those corners with the land test as line of sight. The
  exact path starts at the receiver itself, so a receiver moved by millimetres moves
  `distance_m` by millimetres.
- `control = list(path = "grid")` keeps the straight line where it stays in water and, where
  land lies on it, takes the path between cell centres in steps to the eight neighbouring
  cells, by `terra::gridDist()`. On the 70 receivers of `stoney_track` it gives the straight
  line on all 92,798 clear pairs. On the 244,392 blocked pairs it gives 0.664 to 1.164 times
  the exact path, 1.059 on average, and it reaches the 70 pairs of the one cell the exact path
  counts as unreachable
  (<https://github.com/jakebrownscombe/detfield/blob/main/validation/05_grid_path.R>).
  positionR 1.0.0 and 1.1.0 took the 8-neighbour grid path for every pair; from 1.2.0 the
  straight line is used where no land crosses it. On open water the 8-neighbour path is up to
  8.2% longer than the straight line, 5.5% on average over the 92,728 clear pairs off the
  receiver's own cell.
- On the 70 receivers of `stoney_track`, `water_distances()` takes 3.9 s on one core and the
  grid option 1.3 s. Over the same lake split into 25 m cells it takes 88 s and 28 s. These
  are medians of 5 and 3 runs on one machine
  (<https://github.com/jakebrownscombe/detfield/blob/main/validation/04_memory_time.R>).
- Two water cells that touch only where two land cells meet diagonally are not joined, since the
  land test closes that point; water reached only through such a point is counted in the
  `unreachable` attribute. One Stoney cell is cut off this way. A receiver within the corner
  tolerance of such a point is tested on the cells ahead of it.
- The land test is one geometric rule for every leg: a segment is blocked when it comes deeper
  than a thousandth of a cell side into land, or that close to a closed gap it does not start
  in. A receiver within 10 cm of a grid line is blocked by land lying on both sides of that
  line, and a receiver just outside the tolerance of a closed gap sees its own cell.
  positionR from 1.2.0 tested 50 points along each line whatever its length, which misses land
  narrower than the gap between points.
- `de_field()` gives every water cell a detection probability for every receiver, read at the
  distance through water and the receiver's depth, and counts in a message the pairs with land
  on the line of sight.
- `system_de()` gives, per day, the probability that at least `k` listening receivers detect a
  transmission from each cell. It accepts the listening table that
  `detqc::as_detfield_input()` returns. A day on which every receiver is `"off"` or `"deaf"`
  gives 0 in every cell. A day with no listening receiver for any other reason, such as
  `"not_downloaded"`, is left out, because its receivers may have listened with their data
  missing. On a day with some listening receivers, the coverage comes from those receivers.
  A receiver the field does not hold is an error on any row of `listening`, listening or not.
- `system_de(transmissions = )` gives the same probability for a bin in which the tag sent that
  many transmissions: the probability that the bin is detected from each cell. The default, 1,
  gives the per-transmission result unchanged.

## Location surfaces

- `bin_detections()` counts each animal's detections per hour (or day) at every listening
  receiver, silent receivers included, and the transmissions the tag sent from its delay. It
  reads listening as receiver-days, deployments or both; a receiver listens in a bin when its
  deployment covers at least `min_cover` of it. Animals with two or more codes need
  `codes = "shared"` or `"independent"`. The attribute `rate_check` flags animals whose
  transmission count looks doubled or halved.
- `location_surface()` turns each animal's counts in a bin into the probability that it was in
  each water cell. `likelihood = "binomial"`, the default, uses detections and misses at every
  listening receiver; `"presence_only"` uses the detections alone, as a check when the curve's
  level is not trusted. With `overdispersion = "station_hour"` the spread is integrated by
  Gauss-Hermite quadrature. The bin table gives a misfit test per bin (`misfit_p`), the
  distance from the most probable cell to the detecting receivers (`mode_dist_m`) and the
  probability near the receiver with the most detections (`near_mass`).
- A run of `location_surface()` in which every bin is degenerate stops with an error of class
  `detfield_error_degenerate`; a stop on one bin has class `detfield_error_degenerate_bin`.
  Both also have class `detfield_error`, so a caller can catch them by class.
- `credible_region()` gives the smallest set of cells holding a level of each bin's
  probability, with its area, as a table or a raster. `sample_locations()` draws cells in
  proportion to their probability; its attribute `estimand` says that the draws describe
  detected bins only.
- The print of a surface counts bins with misfit, bins with almost no probability near the
  receiver that heard the animal most, and the distance of the modes from the detecting
  receivers. The radius of the near count is, per receiver, the distance to the farthest cell
  centre at which the field gives a detection probability of 0.5 or more.

## Used and available points

- `sample_available()` draws availability points for each animal, detected bin and presence
  draw from the whole detection surface of the bin's date: in proportion to the probability that
  the bin is detected from each cell (`system_de(transmissions = )`), or uniformly over the cells
  at or above a probability (`uniform_above`). positionR 1.x drew availability from the same
  surface but weighted by the per-transmission probability, kept cells strictly above a minimum
  of 0.001 and below an optional maximum, and drew one set per fish and time period. Here the
  weight is the coverage table given, with no minimum or maximum, the
  uniform mode keeps cells at or above the threshold, and each presence draw gets its own set,
  for multiple imputation. Pass the hourly table: an hour is a presence when the animal was heard
  at least once in it. With the per-transmission weight, simulated animals with no preference
  showed one (-0.091 per m of bottom depth on their true tracks, 95% CI -0.161 to -0.020).
  Strata with no cell to draw are dropped, counted in a warning and listed in the attribute
  `dropped`.
- `used_available()` stacks the used and available points into one table for a conditional
  logistic regression: one stratum per animal, bin and draw with its used point (`case = 1`)
  and its available points (`case = 0`), and one column per covariate layer read at each
  point's cell. Rows with a missing covariate are dropped and counted by reason in the
  attribute `dropped`. It fits no model.

## Figures

- `plot()` of a detection model draws the curves from 0 m to the distance where the last curve
  falls below 0.01.
- `plot()` of a detection field maps one receiver, or the probability that any receiver
  detects a transmission, on the magma scale, with receivers as yellow rings and the scale bar
  at the bottom right.
- `plot()` of a location surface maps one animal and bin with the outlines of its 50% and 95%
  credible regions and the receivers listening in the bin, drawn by `layer_receivers()`, on a
  map at least 1.5 km across. Only receivers inside the map are drawn and keyed. The colour
  runs from 0 to the bin's largest cell probability, which is the top label of the key. The
  title gives the animal and bin on one line.
- `layer_receivers()` draws receivers on a map as rings over a dark ring: yellow for a receiver
  that heard the animal, sized by its count when a count column is given, and white for a
  listening receiver that heard none.
- `theme_positionr()`, `scale_fill_depth()`, `scale_fill_probability()` and
  `annotation_scale_bar()` draw every figure of the positionR packages in one style.
  `scale_fill_depth()` reads depth stored as positive or negative metres.

## Example data

- `stoney_track`: five simulated animals in Stoney Lake, Ontario, for 48 hours from 1 June
  2024, with 70 receivers on a grid, every transmission heard or missed, and the true tracks.
  The array is dense enough that most detected hours are heard by two or more receivers.
- `stoney_selection`: forty simulated animals on the same 70 receivers with a known preference
  for bottom depth (twenty for shallow water, twenty with none), their true tracks, the
  position draws of every detected hour, the coverage of the array per day, and the per-animal
  results and hourly depths of a trial over five preferences, two arrays and two blocks of
  animals.
- `stoney_range_test`: a range test at receiver depths of 1 m and 25 m, counted by
  `range_test_counts()` from the files included in `inst/extdata`.
- `stoney_range_covariates`: a range test with two tag models, two depths and daily wind.
- The Stoney Lake depth raster and the study's deployments, detections and tags are included as
  files, so every example reads the same inputs a user's own study gives.

## Articles

- "Habitat selection: depth use from hourly positions" sets the depth of position draws against
  the true depth of 40 simulated animals, then draws availability points and reports a trial of
  the depth coefficient from positions against the true track. "Hourly positions against the
  true track" links to it for depth use.
- The articles are "Range models: bring your own, or start from two distances", "Hourly
  positions against the true track", "Habitat selection: depth use from hourly positions" and
  "Position fish from a cleaned detection set".
- The walk through the package's steps is the Get started vignette, `vignette("detfield")`;
  the home page holds one example. "Position fish from a cleaned detection set" is on the
  website only, as it needs detqc installed.
