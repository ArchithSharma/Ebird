sampling <- read_sampling("/Users/archith/Downloads/amekes_FLORIDA_SAMPLING.txt")
kestrels <- read_ebd("/Users/archith/Downloads/amekes_FLORIDA.txt")

totaldata = auk_zerofill(x = kestrels, sampling_events = sampling, collapse = TRUE)

sampling <- sampling |> 
  filter(all_species_reported,
         between(year(observation_date), 1976, 2022),
         month(observation_date) == 6)

# filter the observation data
kestrels <- kestrels |> 
  filter(all_species_reported,
         between(year(observation_date), 1976, 2022),
         month(observation_date) == 6)
sampling_sf <- sampling |> 
  dplyr::select(checklist_id, latitude, longitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# study_region_buffered <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", layer = "ne_states") |>
#   filter(state_code %in% c("US-DC", "US-NY", "US-PA", "US-CT", "US-DE", "US-FL", "US-GA", "US-ME", "US-MD", "US-MA", "US-NH", "US-NJ", "US-NC", "US-RI", "US-SC", "US-VT", "US-VA", "US-WV")) |>
#   st_transform(crs = st_crs(sampling_sf)) |>
#   st_buffer(dist = 1000)

study_region_buffered <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", layer = "ne_states") |>
 filter(state_code %in% c("US-DC", "US-NY", "US-PA", "US-CT", "US-DE", "US-FL", "US-GA", "US-ME", "US-MD", "US-MA", "US-NH", "US-NJ", "US-NC", "US-RI", "US-SC", "US-VT", "US-VA", "US-WV")) |>
st_transform(crs = st_crs(sampling_sf)) |>
 st_buffer(dist = 1000)
in_region <- sampling_sf[study_region_buffered, ]
kestrels <- semi_join(kestrels, in_region, by = "checklist_id")
sampling <- semi_join(sampling, in_region, by = "checklist_id")


time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
totaldata <- totaldata |> 
  dplyr::mutate(
    # convert count to integer and X to NA
    # ignore the warning "NAs introduced by coercion"
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for stationary counts
    effort_distance_km = if_else(protocol_type == "Stationary", 
                                 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

totaldata <- totaldata |> 
  filter(protocol_type %in% c("Stationary", "Traveling", "Historical"),
         effort_hours <= 6,
         effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10)

ne_land <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", "ne_land") |> 
  st_geometry()
ne_country_lines <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", "ne_country_lines") |> 
  st_geometry()
ne_state_lines <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", "ne_state_lines") |> 
  st_geometry()
study_region <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", "ne_states") |> 
  filter(state_code %in% c("US-FL")) |> 
  st_geometry()
checklists_sf <- totaldata |> 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |> 
  dplyr::select(species_observed)
# map
par(mar = c(0.25, 0.25, 4, 0.25))
# set up plot area
plot(st_geometry(checklists_sf), 
     main = "American Kestrel eBird Observations\nJune 1976-2022",
     col = NA, border = NA)
# contextual gis data
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)
plot(study_region, col = "#e6e6e6", border = NA, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
# not observed
plot(filter(checklists_sf, !species_observed),
     pch = 19, cex = 0.1, col = alpha("#111111", 0.25),
     add = TRUE)
# observed
plot(filter(checklists_sf, species_observed),
     pch = 19, cex = 0.3, col = alpha("#8894ae", 0.25),
     add = TRUE)
# legend
legend("bottomright", bty = "n",
       col = c("#111111", "#8894ae"),
       legend = c("eBird checklist", "Kestrel sighting"),
       pch = 19)
box()

# summarize data by hour long bins
breaks <- seq(0, 6)
labels <- breaks[-length(breaks)] + diff(breaks) / 2
checklists_duration <- totaldata |> 
  mutate(duration_bins = cut(effort_hours, 
                             breaks = breaks, 
                             labels = labels,
                             include.lowest = TRUE),
         duration_bins = as.numeric(as.character(duration_bins))) |> 
  group_by(duration_bins) |> 
  summarise(n_checklists = n(),
            n_detected = sum(species_observed),
            det_freq = mean(species_observed))

# histogram
g_duration_hist <- ggplot(checklists_duration) +
  aes(x = duration_bins, y = n_checklists) +
  geom_segment(aes(xend = duration_bins, y = 0, yend = n_checklists),
               color = "grey50") +
  geom_point() +
  scale_x_continuous(breaks = breaks) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Checklist duration [hours]",
       y = "# checklists",
       title = "Distribution of checklist durations")

# frequency of detection
g_duration_freq <- ggplot(checklists_duration |> filter(n_checklists > 100)) +
  aes(x = duration_bins, y = det_freq) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = breaks) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Checklist duration [hours]",
       y = "% checklists with detections",
       title = "Detection frequency")

# combine
grid.arrange(g_duration_hist, g_duration_freq)

# summarize data by 1 km bins
breaks <- seq(0, 10)
labels <- breaks[-length(breaks)] + diff(breaks) / 2
checklists_dist <- totaldata |> 
  mutate(dist_bins = cut(effort_distance_km, 
                         breaks = breaks, 
                         labels = labels,
                         include.lowest = TRUE),
         dist_bins = as.numeric(as.character(dist_bins))) |> 
  group_by(dist_bins) |> 
  summarise(n_checklists = n(),
            n_detected = sum(species_observed),
            det_freq = mean(species_observed))

# histogram
g_dist_hist <- ggplot(checklists_dist) +
  aes(x = dist_bins, y = n_checklists) +
  geom_segment(aes(xend = dist_bins, y = 0, yend = n_checklists),
               color = "grey50") +
  geom_point() +
  scale_x_continuous(breaks = breaks) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Distance travelled [km]",
       y = "# checklists",
       title = "Distribution of distance travelled")

# frequency of detection
g_dist_freq <- ggplot(checklists_dist |> filter(n_checklists > 100)) +
  aes(x = dist_bins, y = det_freq) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = breaks) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Distance travelled [km]",
       y = "% checklists with detections",
       title = "Detection frequency")

# combine
grid.arrange(g_dist_hist, g_dist_freq)
landcover2 = rast("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data-raw/landcover_mcd12q1_umd_us-ga_2014-2022.tif")
plot(as.factor(landcover[["NA_NALCMS_landcover_2020_30m"]]))

checklists_sf <- checklists |> 
  # identify unique location/year combinations
  distinct(locality_id, year_lc, latitude, longitude) |> 
  # generate a 10 km neighborhoods
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
buffers <- st_buffer(checklists_sf, dist = set_units(5, "km"))

  lsm <- NULL
  for (i in seq_len(nrow(buffers))) {
    buffer_i <- st_transform(buffers[i, ], crs = crs(landcover))
    year <- as.character(buffer_i$year_lc)
    
    # crop and mask landcover raster so all values outside buffer are missing
    lsm[[i]] <- crop(landcover[[year]], buffer_i) |> 
      mask(buffer_i) |> 
      # calcualte landscape metrics
      calculate_lsm(level = "class", metric = c("pland", "ed")) |> 
      # add variables to uniquely identify each point
      mutate(locality_id = buffer_i$locality_id, 
             year_lc = buffer_i$year_lc) |> 
      dplyr::select(locality_id, year_lc, class, metric, value)
  }
  lsm <- bind_rows(lsm)
  
  lsm_wide <- lsm |> 
    # fill missing classes with zeros
    complete(nesting(locality_id, year_lc),
             class = lc_classes$class,
             metric = c("ed", "pland"),
             fill = list(value = 0)) |> 
    # bring in more descriptive names
    inner_join(dplyr::select(lc_classes, class, label), by = "class") |> 
    # transform from long to wide format
    pivot_wider(values_from = value,
                names_from = c(class, label, metric),
                names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
                names_sort = TRUE) |> 
    arrange(locality_id, year_lc)
  
  
  elev_buffer <- exact_extract(elevation, buffers, fun = c("mean", "stdev"),
                               progress = FALSE) |> 
    # add variables to uniquely identify each point
    mutate(locality_id = buffers$locality_id, year_lc = buffers$year_lc) |> 
    dplyr::select(locality_id, year_lc, 
           elevation_mean = mean,
           elevation_sd = stdev)
  study_region <- read_sf("/Users/archith/Downloads/CBCSpreadsheets/ebird-best-practices-data/data/gis-data.gpkg", layer = "ne_states") |> 
    filter(state_code == c("US-FL")) |> 
    st_transform(crs = laea_crs)
  
  # create a raster template covering the region with 3 km resolution
  r <- rast(study_region, res = c(10000, 10000))
  
  # fill the raster with 1s inside the study region
  r <- rasterize(study_region, r, values = 1) |> 
    setNames("study_region")
  
  elevation <- rast("/Users/archith/Downloads/Landcover 2001-20/floridaelevation2.tif")
  
  # mean and standard deviation within each circular neighborhood
  elev_buffer <- exact_extract(elevation, buffers, fun = c("mean", "stdev"),
                               progress = FALSE) |> 
    # add variables to uniquely identify each point
    mutate(locality_id = buffers$locality_id, year_lc = buffers$year_lc) |> 
    dplyr::select(locality_id, year_lc, 
           elevation_mean = mean,
           elevation_sd = stdev)
  
  env_variables <- inner_join(elev_buffer, lsm_wide,
                              by = c("locality_id", "year_lc"))
  
  # attach and expand to checklists
  env_variables <- checklists |> 
    dplyr::select(checklist_id, locality_id, year_lc) |> 
    inner_join(env_variables, by = c("locality_id", "year_lc")) |> 
    dplyr::select(-locality_id, -year_lc)
  
  # save to csv, dropping any rows with missing variables
  write_csv(drop_na(env_variables), 
            "environmentalvarskestrels10km.csv", 
            na = "")
  buffers_pg <- as.data.frame(r, cells = TRUE, xy = TRUE) |> 
    dplyr::select(cell_id = cell, x, y) |> 
    st_as_sf(coords = c("x", "y"), crs = laea_crs, remove = FALSE) |> 
    st_transform(crs = 4326) |> 
    st_buffer(set_units(10, "km"))
  
  lsm_pg <- NULL
  for (i in seq_len(nrow(buffers_pg))) {
    buffer_i <- st_transform(buffers_pg[i, ], crs = crs(landcover))
    
    # crop and mask landcover raster so all values outside buffer are missing
    lsm_pg[[i]] <- crop(landcover[["florida2"]], buffer_i) |> 
      mask(buffer_i) |> 
      # calcualte landscape metrics
      calculate_lsm(level = "class", metric = c("pland", "ed")) |> 
      # add variable to uniquely identify each point
      mutate(cell_id = buffer_i$cell_id) |> 
      dplyr::select(cell_id, class, metric, value)
  }
  lsm_pg <- bind_rows(lsm_pg)
  
  # transform to wide format
  lsm_wide_pg <- lsm_pg |> 
    # fill missing classes with zeros
    complete(cell_id,
             class = lc_classes$class,
             metric = c("ed", "pland"),
             fill = list(value = 0)) |> 
    # bring in more descriptive names
    inner_join(dplyr::select(lc_classes, class, label), by = "class") |> 
    # transform from long to wide format
    pivot_wider(values_from = value,
                names_from = c(class, label, metric),
                names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
                names_sort = TRUE,
                values_fill = 0) |> 
    arrange(cell_id)
  elev_buffer_pg <- exact_extract(elevation, buffers_pg, 
                                  fun = c("mean", "stdev"),
                                  progress = FALSE) |> 
    # add variables to uniquely identify each point
    mutate(cell_id = buffers_pg$cell_id) |> 
    dplyr::select(cell_id, elevation_mean = mean, elevation_sd = stdev)
  # combine landcover and elevation
  env_variables_pg <- inner_join(elev_buffer_pg, lsm_wide_pg, by = "cell_id")
  
  # attach the xy coordinates of the cell centers
  env_variables_pg <- buffers_pg |> 
    st_drop_geometry() |> 
    dplyr::select(cell_id, x, y) |> 
    inner_join(env_variables_pg, by = "cell_id")
  
  # save as csv, dropping any rows with missing variables
  write_csv(drop_na(env_variables_pg),
            "/Users/archith/Downloads/environmental-variables_prediction-grid_us-fl10km.csv", 
            na = "")

totaldata$type <- if_else(runif(nrow(totaldata)) <= 0.8, "train", "test")
table(totaldata$type) / nrow(totaldata)
checklists <- totaldata |> 
  dplyr::select(checklist_id, observer_id, type,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         hours_of_day, 
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers)
write.csv(checklists, "/Users/archith/Downloads/kestrelfull.csv")
# summarize data by hourly bins
breaks <- seq(0, 24)
labels <- breaks[-length(breaks)] + diff(breaks) / 2
checklists_time <- checklists |> 
  mutate(hour_bins = cut(hours_of_day, 
                         breaks = breaks, 
                         labels = labels,
                         include.lowest = TRUE),
         hour_bins = as.numeric(as.character(hour_bins))) |> 
  group_by(hour_bins) |> 
  summarise(n_checklists = n(),
            n_detected = sum(species_observed),
            det_freq = mean(species_observed))

# histogram
g_tod_hist <- ggplot(checklists_time) +
  aes(x = hour_bins, y = n_checklists) +
  geom_segment(aes(xend = hour_bins, y = 0, yend = n_checklists),
               color = "grey50") +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Hours since midnight",
       y = "# checklists",
       title = "Distribution of observation start times")

# frequency of detection
g_tod_freq <- ggplot(checklists_time |> filter(n_checklists > 100)) +
  aes(x = hour_bins, y = det_freq) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Hours since midnight",
       y = "% checklists with detections",
       title = "Detection frequency")

# combine
grid.arrange(g_tod_hist, g_tod_freq)