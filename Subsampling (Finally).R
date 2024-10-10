env_vars = read.csv("environmentalvarskestrels30km.csv")
checklists_env = checky |> inner_join(env_vars, by = "checklist_id")
pred_grid = read.csv("/Users/archith/Downloads/environmental-variables_prediction-grid_us-fl.csv")
crs <- st_crs(r)
checklists_ss <- grid_sample_stratified(checklists,
                                        obs_column = "species_observed",
                                        sample_by = "type", res = c(10000, 10000, 30), 
                                        min_detection_probability = 0)

#checklists_ss30km = read.csv("10kmsample.csv")
write.csv(checklists_ss, "10kmsample.csv")
View(checklists_ss)
checklists_ss = drop_na(checklists_ss)
all_pts <- checkliststest |>  
  filter(type == "train") |> 
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) |>
  st_transform(crs = crs) |> 
  dplyr::select(species_observed)
ss_pts <- checklists_ss |>  
  filter(type == "train") |> 
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) |>
  st_transform(crs = crs) |> 
  dplyr::select(species_observed)
both_pts <- list(before_ss = all_pts, after_ss = ss_pts)

# map
p <- par(mfrow = c(1, 2))
for (i in seq_along(both_pts)) {
  par(mar = c(0.25, 0.25, 0.25, 0.25))
  # set up plot area
  plot(st_geometry(both_pts[[i]]), col = NA)
  # contextual gis data
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  plot(study_region, col = "#cccccc", border = NA, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  # ebird observations
  # not observed
  plot(st_geometry(both_pts[[i]]),
       pch = 19, cex = 0.1, col = alpha("#111111", 0.25),
       add = TRUE)
  # observed
  plot(filter(both_pts[[i]], species_observed) |> st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#8894ae", 0.5),
       add = TRUE)
  # legend
  legend("bottomright", bty = "n",
         col = c("#111111", "#8894ae"),
         legend = c("Non-detection", "Detection"),
         pch = 19)
  box()
  par(new = TRUE, mar = c(0, 0, 3, 0))
  if (names(both_pts)[i] == "before_ss") {
    title("Before subsampling")
  } else {
    title("After subsampling")
  }
}
par(p)

checklists_train <- checklists_ss30km |> 
  filter(type == "train") |> 
  # select only the columns to be used in the model
  dplyr::select(species_observed,
         year, day_of_year, hours_of_day,
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers, 
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"))
detection_freq <- mean(checklists$species_observed)

er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                   data = checklists_train,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE, 
                   sample.fraction = c(detection_freq, detection_freq))
# predicted encounter rate based on out of bag samples
er_pred <- er_model$predictions[, 2]
# observed detection, converted back from factor
det_obs <- as.integer(checklists_train$species_observed)
# construct a data frame to train the scam model
obs_pred <- data.frame(obs = det_obs, pred = er_pred)

# train calibration model
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                          gamma = 2,
                          data = obs_pred)
er_breaks <- seq(0, 1, by = 0.02)
mean_er <- obs_pred |>
  mutate(er_bin = cut(pred, breaks = er_breaks, include.lowest = TRUE)) |>
  group_by(er_bin) |>
  summarise(n_checklists = n(),
            pred = mean(pred), 
            obs = mean(obs),
            .groups = "drop")

# make predictions from the calibration model
calibration_curve <- data.frame(pred = er_breaks)
cal_pred <- predict(calibration_model, calibration_curve, type = "response")
calibration_curve$calibrated <- cal_pred

# compared binned mean encounter rates to calibration model
ggplot(calibration_curve) +
  aes(x = pred, y = calibrated) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(color = "blue") +
  geom_point(data = mean_er, 
             aes(x = pred, y = obs),
             size = 2, alpha = 0.6,
             show.legend = FALSE) +
  labs(x = "Estimated encounter rate",
       y = "Observed encounter rate",
       title = "Calibration model") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))


mcc_f1 <- mccf1(
  # observed detection/non-detection
  response = obs_pred$obs,
  # predicted encounter rate from random forest
  predictor = obs_pred$pred)

# identify best threshold
mcc_f1_summary <- summary(mcc_f1)
#>  mccf1_metric best_threshold
#>         0.399          0.508
threshold <- mcc_f1_summary$best_threshold[1]

# get the test set held out from training
checklists_test <- filter(checklists_ss, type == "test") |> 
  mutate(species_observed = as.integer(species_observed))

# predict to test data using random forest model
pred_er <- predict(er_model, data = checklists_test, type = "response")
# extract probability of detection
pred_er <- pred_er$predictions[, 2]
# convert predictions to binary (presence/absence) using the threshold
pred_binary <- as.integer(pred_er > threshold)
# calibrate
pred_calibrated <- predict(calibration_model, 
                           newdata = data.frame(pred = pred_er), 
                           type = "response") |> 
  as.numeric()
# constrain probabilities to 0-1
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1
# combine observations and estimates
obs_pred_test <- data.frame(id = seq_along(pred_calibrated),
                            # actual detection/non-detection
                            obs = as.integer(checklists_test$species_observed),
                            # binary detection/on-detection prediction
                            pred_binary = pred_binary,
                            # calibrated encounter rate
                            pred_calibrated = pred_calibrated)

# mean squared error (mse)
mse <- mean((obs_pred_test$obs - obs_pred_test$pred_calibrated)^2, na.rm = TRUE)

# precision-recall auc
em <- precrec::evalmod(scores = obs_pred_test$pred_binary, 
                       labels = obs_pred_test$obs)
pr_auc <- precrec::auc(em) |> 
  filter(curvetypes == "PRC") |> 
  pull(aucs)

# calculate metrics for binary prediction: sensitivity, specificity
pa_metrics <- obs_pred_test |> 
  dplyr::select(id, obs, pred_binary) |> 
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# mcc and f1
mcc_f1 <- calculate_mcc_f1(obs_pred_test$obs, obs_pred_test$pred_binary)

# combine ppms together
ppms <- data.frame(
  mse = mse,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  pr_auc = pr_auc,
  mcc = mcc_f1$mcc,
  f1 = mcc_f1$f1
)
knitr::kable(pivot_longer(ppms, everything()), digits = 3)

# extract predictor importance from the random forest model object
pred_imp <- er_model$variable.importance
pred_imp <- data.frame(predictor = names(pred_imp), 
                       importance = pred_imp) |> 
  arrange(desc(importance))
# plot importance of top 20 predictors
ggplot(head(pred_imp, 20)) + 
  aes(x = reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL, 
       y = "Predictor Importance (Gini Index)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5))

# function to calculate partial dependence for a single predictor
calculate_pd <- function(predictor, er_model, calibration_model,
                         data, x_res = 25, n = 1000) {
  # create prediction grid using quantiles
  x_grid <- quantile(data[[predictor]],
                     probs = seq(from = 0, to = 1, length = x_res),
                     na.rm = TRUE)
  # remove duplicates
  x_grid <- x_grid[!duplicated(signif(x_grid, 8))]
  x_grid <- unname(unique(x_grid))
  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  data <- data[sample(seq.int(nrow(data)), size = n, replace = FALSE), ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict encounter rate
  p <- predict(er_model, data = grid)
  
  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         .groups = "drop")
  
  # calibrate
  pd$encounter_rate <- predict(calibration_model, 
                               newdata = data.frame(pred = pd$encounter_rate), 
                               type = "response")
  pd$encounter_rate <- as.numeric(pd$encounter_rate)
  # constrain to 0-1
  pd$encounter_rate[pd$encounter_rate < 0] <- 0
  pd$encounter_rate[pd$encounter_rate > 1] <- 1
  
  return(pd)
}

#Now we’ll use this function to calculate partial dependence for the top 6 predictors.

# calculate partial dependence for each of the top 6 predictors
pd <- NULL
for (predictor in head(pred_imp$predictor)) {
  pd <- calculate_pd(predictor, 
                     er_model = er_model, 
                     calibration_model = calibration_model,
                     data = checklists_train) |> 
    bind_rows(pd)
}
head(pd)

# plot partial dependence
ggplot(pd) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor, levels = rev(unique(predictor))), 
             ncol = 2, scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

# find peak time of day from partial dependence
pd_time <- calculate_pd("hours_of_day",
                        er_model = er_model, 
                        calibration_model = calibration_model,
                        data = checklists_train) |> 
  dplyr::select(hours_of_day = x, encounter_rate)

# histogram
g_hist <- ggplot(checklists_train) +
  aes(x = hours_of_day) +
  geom_histogram(binwidth = 1, center = 0.5, color = "grey30",
                 fill = "grey50") +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Hours since midnight",
       y = "# checklists",
       title = "Distribution of observation start times")

# partial dependence plot
g_pd <- ggplot(pd_time) +
  aes(x = hours_of_day, y = encounter_rate) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(x = "Hours since midnight",
       y = "Encounter rate",
       title = "Observation start time partial dependence")

# combine
grid.arrange(g_hist, g_pd) #Peak time 8:05 AM

pred_grid_eff <- pred_grid |> 
  mutate(observation_date = ymd("2000-06-15"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         hours_of_day = t_peak,
         effort_distance_km = 2,
         effort_hours = 1,
         effort_speed_kmph = 2,
         number_observers = 1)

pred_er <- predict(er_model, data = pred_grid_eff, type = "response")
pred_er <- pred_er$predictions[, 2]
# define range-boundary
pred_binary <- as.integer(pred_er > threshold)
# apply calibration
pred_er_cal <- predict(calibration_model, 
                       data.frame(pred = pred_er), 
                       type = "response") |> 
  as.numeric()
# constrain to 0-1
pred_er_cal[pred_er_cal < 0] <- 0
pred_er_cal[pred_er_cal > 1] <- 1
# combine predictions with coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          in_range = pred_binary, 
                          encounter_rate = pred_er_cal)

r_pred <- predictions |> 
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs) |> 
  dplyr::select(in_range, encounter_rate) |> 
  # rasterize
  rasterize(r, field = c("in_range", "encounter_rate"))
print(r_pred)


par(mar = c(0.25, 0.25, 1.25, 0.25))
# set up plot area
plot(study_region, 
     main = "Wood Thrush Range (June 2023)",
     col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# convert binary prediction to categorical
r_range <- as.factor(r_pred[["in_range"]])
plot(r_range, col = c("#e6e6e6", "forestgreen"),
     maxpixels = ncell(r_range),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region, col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks
brks <- global(r_pred[["encounter_rate"]], fun = quantile, 
               probs = seq(0, 1, 0.1), na.rm = TRUE) |> 
  as.numeric() |> 
  unique()
# label the bottom, middle, and top value
lbls <- round(c(0, median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_pred[["encounter_rate"]], 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Kestrel Encounter Rate (MODIS Landcover)"
image.plot(zlim = c(0, 1), legend.only = TRUE, 
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.03, 0.06),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))