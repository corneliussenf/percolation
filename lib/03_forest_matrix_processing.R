
library(terra)
library(sf)
library(landscapemetrics)
library(tidyverse)

forest <- rast("data/forestlanduse_mask_EUmosaic_3035.tif") # !!! NEEDS TO BE DOWNLOADED AND MOSAICED FROM: https://doi.org/10.5281/zenodo.8389085
disturbance <- rast("data/latest_disturbance_eu_v211_3035.tif") # !!! NEEDS TO BE DOWNLOADED AND MOSAICED FROM: https://doi.org/10.5281/zenodo.8389085

disturbance_summary <- read_csv("data/disturbance_summary.csv")

disturbance_gaps <- classify(disturbance, matrix(c(NA, NA, 0, 0, 2023 - 30, 0, 2023 - 30, 2024, 1), ncol = 3, byrow = TRUE))

countries <- read_sf("data/gis/countries_europe_map.shp")
countries <- st_buffer(countries, 0)

grid <- read_sf("data/gis/grid_100km_epsg3035.shp")

samp <- 1:nrow(grid)

landscape_metrics <- vector("list", length(samp))

k <- 0

for (i in samp) {
  
  k <- k + 1
  
  print(paste0(round(k / length(samp) * 100, 0), "%"))
  
  grid_tmp <- grid[i, ]
  
  forest_tmp <- crop(forest, grid_tmp)
  disturbance_gaps_tmp <- crop(disturbance_gaps, grid_tmp)
  
  forest_matrix_tmp <- forest_tmp - disturbance_gaps_tmp
  
  forest_matrix_binary <- classify(forest_matrix_tmp, matrix(c(NA, NA, 0, NA, 1, 1), ncol = 2, byrow = TRUE))
  
  np <- lsm_l_np(forest_matrix_binary, directions = 4)
  ai <- lsm_l_ai(forest_matrix_binary, directions = 4)
  
  area_patch <- lsm_p_area(forest_matrix_binary, directions = 4)
  
  area_patch <- area_patch %>%
    summarize(
      area_mn = mean(value, na.rm = TRUE),
      area_md = median(value, na.rm = TRUE),
      area_max = max(value, na.rm = TRUE),
      area_q95 = quantile(value, 0.95, na.rm = TRUE),
      area_q99 = quantile(value, 0.99, na.rm = TRUE),
      area_sm = sum(value, na.rm = TRUE)
    ) %>%
    mutate(
      lpi_q99 = area_q99 / area_sm,
      lpi_max = area_max / area_sm
    )
  
  landscape_metrics[[k]] <- data.frame(
    np = np$value,
    ai = ai$value,
    area_patch
  )
  
}

dat <- landscape_metrics %>%
  bind_rows() %>%
  mutate(
    gridid = grid[samp, ]$gridid
  ) %>%
  pivot_longer(
    cols = 1:(ncol(.) - 1)
  ) %>%
  left_join(
    disturbance_summary
  )

write_csv(dat, "data/forest_matrix.csv")
