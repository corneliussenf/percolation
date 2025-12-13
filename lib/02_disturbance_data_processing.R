
library(tidyverse)
library(terra)
library(sf)
library(exactextractr)
library(landscapemetrics)
library(foreach)
library(doParallel)

# Countries/grids to select
countries <- read_sf("data/gis/countries_europe_map.shp")
countries <- st_buffer(countries, 0)
grid <- read_sf("data/gis/grid_100km_epsg3035.shp")

# Load disturbances and characterize change patterns annually ----

recovery <- 30 # years to recovery

distdat <- list.files(
  "data/mosaics_annual_disturb_v211", # !!! NEEDS TO BE DOWNLOADED AND MOSAICED FROM: https://doi.org/10.5281/zenodo.8389085
  pattern = ".tif$",
  full.names = TRUE
)

distdat_agent <- list.files(
  "data/mosaics_agents/v211/", # !!! NEEDS TO BE DOWNLOADED AND MOSAICED FROM: https://doi.org/10.5281/zenodo.8389085
  pattern = ".tif$",
  full.names = TRUE
)

years <- basename(distdat) %>%
  substring(., 1, 4)

n.cores <- parallel::detectCores() - 10

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

dat_full_summary_collector <- foreach (j = 1:length(unique(grid$gridid)), .packages = c("tidyverse", "terra", "sf", "landscapemetrics")) %dopar% {
  
  grid_sel <- grid[grid$gridid == unique(grid$gridid)[j], ]
  
  if (!file.exists(paste0("temp/dat_full_summary_", unique(grid_sel$gridid), ".csv"))) {
    
    distdat_tmp1 <- rast(distdat[1])
    grid_sel <- st_transform(grid_sel, st_crs(distdat_tmp1))
    distdat_tmp1 <- crop(distdat_tmp1, grid_sel)
    distdat_tmp1 <- classify(distdat_tmp1, matrix(c(NA, 0, 1, 1), ncol = 2, byrow = TRUE))
    
    out <- vector("list", length(distdat) - 1)
    out_natural <- vector("list", length(distdat) - 1)
    out_harvest <- vector("list", length(distdat) - 1)
    
    counter <- distdat_tmp1
    
    if (sum(values(distdat_tmp1)) > 0) {
      
      for (i in 2:length(distdat)) {
        
        distdat_tmp2 <- rast(distdat[i])
        distdat_agent_tmp2 <- rast(distdat_agent[i])
        
        year2 <- years[i]
        
        print(year2)
        
        distdat_tmp2 <- crop(distdat_tmp2, grid_sel)
        distdat_tmp2 <- classify(distdat_tmp2, matrix(c(NA, 0, 1, 1), ncol = 2, byrow = TRUE))
        
        distdat_agent_tmp2 <- crop(distdat_agent_tmp2, grid_sel)
        distdat_natural_tmp2 <- classify(distdat_agent_tmp2, matrix(c(NA, 0, 1, 1, 2, 1, 3, 0), ncol = 2, byrow = TRUE))
        distdat_harvest_tmp2 <- classify(distdat_agent_tmp2, matrix(c(NA, 0, 1, 0, 2, 0, 3, 1), ncol = 2, byrow = TRUE))
        
        distdat_tmp12 <- distdat_tmp1 + distdat_tmp2
        distdat_tmp12 <- classify(distdat_tmp12, matrix(c(0, NA, 1, 1, 2, 1), ncol = 2, byrow = TRUE)) 
        
        distdat_natural_tmp12 <- distdat_tmp1 + distdat_natural_tmp2
        distdat_natural_tmp12 <- classify(distdat_natural_tmp12, matrix(c(0, NA, 1, 1, 2, 1), ncol = 2, byrow = TRUE)) 
        
        distdat_harvest_tmp12 <- distdat_tmp1 + distdat_harvest_tmp2
        distdat_harvest_tmp12 <- classify(distdat_harvest_tmp12, matrix(c(0, NA, 1, 1, 2, 1), ncol = 2, byrow = TRUE)) 
        
        distdat_tmp12_patches <- get_patches(distdat_tmp12, directions = 4)
        distdat_natural_tmp12_patches <- get_patches(distdat_natural_tmp12, directions = 4)
        distdat_harvest_tmp12_patches <- get_patches(distdat_harvest_tmp12, directions = 4)
        
        dat <- rast(c(distdat_tmp12_patches$layer_1, distdat_tmp1, distdat_tmp2))
        dat <- as_tibble(dat)
        
        out[[i-1]] <- dat %>%
          set_names(
            c("patchid", "t0", "t1")
          ) %>%
          filter(
            !is.na(patchid)
          ) %>%
          group_by(
            patchid
          ) %>%
          summarize(
            area_ha_t0 = sum(t0) * 0.09,
            area_ha_t1 = sum(t1) * 0.09
          ) %>%
          ungroup() %>%
          mutate(
            expansion = ifelse(area_ha_t0 > 0 & area_ha_t1 > 0, 1, 0)
          ) %>%
          mutate(
            year = year2
          )
        
        dat_natural <- rast(c(distdat_natural_tmp12_patches$layer_1, distdat_tmp1, distdat_natural_tmp2))
        dat_natural <- as_tibble(dat_natural)
        
        out_natural[[i-1]] <- dat_natural %>%
          set_names(
            c("patchid", "t0", "t1")
          ) %>%
          filter(
            !is.na(patchid)
          ) %>%
          group_by(
            patchid
          ) %>%
          summarize(
            area_ha_t0 = sum(t0) * 0.09,
            area_ha_t1 = sum(t1) * 0.09
          ) %>%
          ungroup() %>%
          mutate(
            expansion = ifelse(area_ha_t0 > 0 & area_ha_t1 > 0, 1, 0)
          ) %>%
          mutate(
            year = year2
          )
        
        dat_harvest <- rast(c(distdat_harvest_tmp12_patches$layer_1, distdat_tmp1, distdat_harvest_tmp2))
        dat_harvest <- as_tibble(dat_harvest)
        
        out_harvest[[i-1]] <- dat_harvest %>%
          set_names(
            c("patchid", "t0", "t1")
          ) %>%
          filter(
            !is.na(patchid)
          ) %>%
          group_by(
            patchid
          ) %>%
          summarize(
            area_ha_t0 = sum(t0) * 0.09,
            area_ha_t1 = sum(t1) * 0.09
          ) %>%
          ungroup() %>%
          mutate(
            expansion = ifelse(area_ha_t0 > 0 & area_ha_t1 > 0, 1, 0)
          ) %>%
          mutate(
            year = year2
          )
        
        distdat_tmp1 <- distdat_tmp12
        distdat_tmp1[is.na(distdat_tmp1)] <- 0
        
        counter <- counter + distdat_tmp1
        distdat_tmp1[counter >= recovery] <- 0
        
      }
      
      dat_full <- out %>%
        bind_rows()
      
      dat_full_summary <- dat_full %>%
        filter(
          area_ha_t1 > 0
        ) %>%
        group_by(year, expansion) %>%
        summarize(
          n = length(unique(patchid)),
          area_ha_t0 = sum(area_ha_t0),
          area_ha_t1 = sum(area_ha_t1)
        ) %>%
        ungroup() %>%
        mutate(
          gridid = unique(grid_sel$gridid)
        )
      
      write_csv(dat_full_summary, paste0("temp/dat_full_summary_", unique(grid_sel$gridid), ".csv"))
      
      dat_natural_full <- out_natural %>%
        bind_rows()
      
      dat_natural_full_summary <- dat_natural_full %>%
        filter(
          area_ha_t1 > 0
        ) %>%
        group_by(year, expansion) %>%
        summarize(
          n = length(unique(patchid)),
          area_ha_t0 = sum(area_ha_t0),
          area_ha_t1 = sum(area_ha_t1)
        ) %>%
        ungroup() %>%
        mutate(
          gridid = unique(grid_sel$gridid)
        )
      
      write_csv(dat_natural_full_summary, paste0("temp/dat_natural_full_summary_", unique(grid_sel$gridid), ".csv"))
      
      dat_harvest_full <- out_harvest %>%
        bind_rows()
      
      dat_harvest_full_summary <- dat_harvest_full %>%
        filter(
          area_ha_t1 > 0
        ) %>%
        group_by(year, expansion) %>%
        summarize(
          n = length(unique(patchid)),
          area_ha_t0 = sum(area_ha_t0),
          area_ha_t1 = sum(area_ha_t1)
        ) %>%
        ungroup() %>%
        mutate(
          gridid = unique(grid_sel$gridid)
        )
      
      write_csv(dat_harvest_full_summary, paste0("temp/dat_harvest_full_summary_", unique(grid_sel$gridid), ".csv"))
      
      list(
        dat_full_summary,
        dat_natural_full_summary,
        dat_harvest_full_summary
      )
      
    } else {
      
      NULL
      
    }
    
  } else {
    
    NULL
    
  }
  
}

parallel::stopCluster(cl = my.cluster)

dat_full_summary_collector <- list.files(
  "temp/",
  pattern = glob2rx("dat_full_summary*"),
  full.names = TRUE
) %>%
  map(read_csv)

dat_full_summary <- dat_full_summary_collector %>%
  discard(~ nrow(.) == 0) %>%
  bind_rows()

write_csv(dat_full_summary, "data/spatial_pattern/dat_full_summary.csv")
dat_full_summary <- read_csv("data/spatial_pattern/dat_full_summary.csv")

dat_natural_full_summary_collector <- list.files(
  "temp/",
  pattern = glob2rx("dat_natural_full_summary*"),
  full.names = TRUE
) %>%
  map(read_csv)

dat_natural_full_summary <- dat_natural_full_summary_collector %>%
  discard(~ nrow(.) == 0) %>%
  bind_rows()

write_csv(dat_natural_full_summary, "data/spatial_pattern/dat_natural_full_summary.csv")
dat_natural_full_summary <- read_csv("data/spatial_pattern/dat_natural_full_summary.csv")

dat_harvest_full_summary_collector <- list.files(
  "temp/",
  pattern = glob2rx("dat_harvest_full_summary*"),
  full.names = TRUE
) %>%
  map(read_csv)

dat_harvest_full_summary <- dat_harvest_full_summary_collector %>%
  discard(~ nrow(.) == 0) %>%
  bind_rows()

write_csv(dat_harvest_full_summary, "data/spatial_pattern/dat_harvest_full_summary.csv")
dat_harvest_full_summary <- read_csv("data/spatial_pattern/dat_harvest_full_summary.csv")

disturbance_summary <- list.files(
  "data/aggregated_disturbances",
  pattern = "disturbance_100km",
  full.names = TRUE
) %>%
  map(
    read_csv
  ) %>%
  bind_rows() %>%
  group_by(
    gridid
  ) %>%
  summarise(
    dist_px_sum = mean(dist_px_sum)
  )

disturbance_summary <- read_csv("data/aggregated_disturbances/forest_100km.csv") %>%
  right_join(
    disturbance_summary
  ) %>%
  mutate(
    rate = dist_px_sum / forest_px_sum
  )

write_csv(disturbance_summary, "data/disturbance_summary_1985_2023.csv")
disturbance_summary <- read_csv("data/disturbance_summary_1985_2023.csv")

years_dat <- list.files(
  "data/aggregated_disturbances",
  pattern = "disturbance_100km"
) %>%
  substr(., 1, 4)

disturbance_summary_30yr_recovery <- list.files(
  "data/aggregated_disturbances",
  pattern = "disturbance_100km",
  full.names = TRUE
) %>%
  map(
    read_csv
  ) %>%
  set_names(
    years_dat
  ) %>%
  bind_rows(
    .id = "year"
  ) %>%
  filter(
    year > (1985 + 30)
  ) %>%
  group_by(
    gridid
  ) %>%
  summarise(
    dist_px_sum = mean(dist_px_sum)
  )

disturbance_summary_30yr_recovery <- read_csv("data/aggregated_disturbances/forest_100km.csv") %>%
  right_join(
    disturbance_summary_30yr_recovery
  ) %>%
  mutate(
    rate = dist_px_sum / forest_px_sum
  )

write_csv(disturbance_summary_30yr_recovery, "data/disturbance_summary_2015_2023.csv")
disturbance_summary_30yr_recovery <- read_csv("data/disturbance_summary_2015_2023.csv")
