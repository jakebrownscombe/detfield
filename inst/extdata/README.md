# Example files

`stoney_depth.tif` has the same grid and depths as `depth_raster` of
positionR 1.7.0, a 100 m depth raster of Stoney Lake, Ontario,
derived from the bathymetry made before 2024, with the CRS set to EPSG:32617
(UTM zone 17N on WGS 84). Depth is zero or negative metres and land is `NA`.
4,818 cells are water, and 56 of them hold 0.

Every other file is simulated. `stoney_deployments.csv`,
`stoney_detections.csv` and `stoney_tags.csv` are the study of
`?stoney_track`, built by
<https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_track.R>.
The three `stoney_range_test_*.csv` files are the range test of
`?stoney_range_test`, built by
<https://github.com/jakebrownscombe/detfield/blob/main/data-raw/stoney_range_test.R>.
Both scripts are in the repository and not in the installed package. The
receivers, stations, tags and animals in them are invented.
