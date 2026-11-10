# Code of the article "Position fish from a cleaned detection set".
# Extracted from vignettes/articles/from-detqc.Rmd by tools/purl_articles.R, a repository file
# the installed package does not carry. The article is the source; this file follows it.

# ---- Column maps ----
# Attach detfield, and ggplot2 for the figures
library(detfield)
library(ggplot2)
# A column map that reads each named column as itself, with its type
same_names <- function(types) {
  data.frame(format = "glatos", target = names(types), source = names(types),
             sep = "", type = unname(types),
             tz = ifelse(types == "datetime", "UTC", ""), note = "")
}
same_names(c(datetime = "datetime", receiver_serial = "character"))

# ---- Read files ----
# Find a file included with detfield
in_file <- function(name) system.file("extdata", name, package = "detfield")
# Read the detections, the receiver deployments and the tags with detqc
detections <- detqc::read_detections(
  in_file("stoney_detections.csv"), format = "glatos",
  col_map = same_names(c(datetime = "datetime", receiver_serial = "character",
                         transmitter = "character", animal_id = "character"))
)
deployments <- detqc::read_deployments(
  in_file("stoney_deployments.csv"), format = "glatos",
  col_map = same_names(c(receiver_serial = "character",
                         station_id = "character",
                         deploy_datetime = "datetime",
                         recover_datetime = "datetime",
                         deploy_lat = "numeric", deploy_lon = "numeric",
                         deploy_z = "numeric"))
)
animals <- detqc::read_animals(
  in_file("stoney_tags.csv"), format = "glatos",
  col_map = same_names(c(transmitter = "character", animal_id = "character",
                         min_delay = "numeric", max_delay = "numeric"))
)
c(detections = nrow(detections), deployments = nrow(deployments),
  animals = nrow(animals))

# ---- Detset ----
# Build the detection set
dets <- detqc::detset(detections, deployments, animals)
dets

# ---- Hand off ----
# Export the detections, the listening table and the tag delays for detfield
detfield_input <- detqc::as_detfield_input(dets)
str(detfield_input, max.level = 1)

# ---- Bin ----
# Count detections per animal, hour and listening receiver
hourly <- do.call(bin_detections, detfield_input)
table(unique(hourly[c("animal_id", "bin_start")])$animal_id)

# ---- Field ----
# Read the lake and measure distances through water from the deployments
lake <- terra::rast(in_file("stoney_depth.tif"))
distances <- water_distances(lake, dets$deployments)
# Build the field from the stand-in curve, with no spread between hours
curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
field <- de_field(distances, de_model(curve, sigma = list(station_hour = 0)))

# ---- Surfaces ----
# Compute one surface per animal and detected hour
surface <- location_surface(hourly, field)
surface

# ---- Heard by ----
# Count the receivers that heard each animal in each hour, three or more
# together
heard <- hourly[hourly$y > 0, ]
heard_by <- aggregate(list(receivers = heard$station_id,
                           detections = heard$y),
                      by = heard[c("animal_id", "bin_start")],
                      FUN = function(v) if (is.numeric(v)) sum(v) else length(v))
heard_by$heard_by <- factor(pmin(heard_by$receivers, 3), 1:3,
                            c("1", "2", "3 or more"))
# Tabulate the hours per animal
hours_table <- table(heard_by$animal_id, heard_by$heard_by)
hours_table
round(mean(heard_by$receivers == 1), 2)

# ---- Heard by plot ----
# Plot the detected hours of each animal by the receivers that heard it
ggplot(heard_by, aes(animal_id, fill = heard_by)) +
  geom_bar() +
  scale_fill_viridis_d(name = "Receivers\nthat heard\nthe animal",
                       option = "magma", begin = 0.3, end = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "Detected hours") +
  theme_positionr()

# ---- Pick hours ----
# Pick stoney-01's hour with the most detections for each count of receivers
st01 <- heard_by[heard_by$animal_id == "stoney-01", ]
picked <- do.call(rbind, lapply(split(st01, st01$heard_by), function(kind) {
  kind[which.max(kind$detections), ]
}))
picked[c("bin_start", "receivers", "detections")]
# The true path of each picked hour, one position a minute
true_path <- function(hour_start) {
  track <- detfield::stoney_track$track
  track[track$animal_id == "stoney-01" & track$datetime >= hour_start &
          track$datetime <= hour_start + 3600, ]
}
# The map window: 1.7 km square on the cell holding the hour's mean true
# position, with its edges on cell edges
around <- function(path, half = 850) {
  centre <- terra::xyFromCell(lake, terra::cellFromXY(lake,
                                                      cbind(mean(path$x),
                                                            mean(path$y))))
  c(centre[1] + c(-half, half), centre[2] + c(-half, half))
}

# ---- Hour many ----
# Map the hour heard by three or more receivers, with its true path
many_hour <- picked$bin_start[picked$heard_by == "3 or more"]
many_path <- true_path(many_hour)
plot(surface, animal_id = "stoney-01", bin_start = many_hour,
     ext = around(many_path)) +
  geom_path(data = many_path, aes(x, y), colour = "white", linewidth = 1.1,
            inherit.aes = FALSE) +
  geom_path(data = many_path, aes(x, y), colour = "grey45", linewidth = 0.5,
            inherit.aes = FALSE)

# ---- Hour two ----
# Map the hour heard by two receivers, with its true path
two_hour <- picked$bin_start[picked$heard_by == "2"]
two_path <- true_path(two_hour)
plot(surface, animal_id = "stoney-01", bin_start = two_hour,
     ext = around(two_path)) +
  geom_path(data = two_path, aes(x, y), colour = "white", linewidth = 1.1,
            inherit.aes = FALSE) +
  geom_path(data = two_path, aes(x, y), colour = "grey45", linewidth = 0.5,
            inherit.aes = FALSE)

# ---- Hour one ----
# Map the single-receiver hour, with its true path
one_hour <- picked$bin_start[picked$heard_by == "1"]
one_path <- true_path(one_hour)
plot(surface, animal_id = "stoney-01", bin_start = one_hour,
     ext = around(one_path)) +
  geom_path(data = one_path, aes(x, y), colour = "white", linewidth = 1.1,
            inherit.aes = FALSE) +
  geom_path(data = one_path, aes(x, y), colour = "grey45", linewidth = 0.5,
            inherit.aes = FALSE)
