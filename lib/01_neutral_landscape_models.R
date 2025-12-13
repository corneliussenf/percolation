
library(tidyverse)
library(terra)
library(NLMR)
library(foreach)
library(doParallel)
library(landscapemetrics)

extent <- 100
grain <- 30
dr <- seq(0.001, 0.05, length.out = 20)
p_forests <- c(0.1, 0.3, 0.5, 0.7, 0.9)
clustering <- c(0.1, 0.3, 0.5)
p_nodists <- c(0, 0.2, 0.4, 0.6)
r <- 30

n_year <- 100 + r

length_out <- length(dr) * length(p_forests) * length(clustering) * length(p_nodists)

reps <- 1:20

n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
clusterEvalQ(cluster, .libPaths("/home/csenf/R/x86_64-pc-linux-gnu-library/4.5"))
registerDoParallel(cluster)

dat <- foreach(i = reps, .packages = c("tidyverse", "terra", "NLMR", "landscapemetrics")) %dopar% {
  
  dr_collector <- c()
  dr_real_collector <- c()
  p_forest_collector <- c()
  y_collector <- c()
  p_nodist_collector <- c()
  p <- c()
  cluster_collector <- c()
  
  landscape_metrics <- vector("list", length_out)
  
  z <- 0
  
  for (p_nodist in p_nodists) {
    
    for (p_forest in p_forests) {
      
      for (clust in clustering) {
        
        for (k in 1:length(dr)) {
          
          print(k)
          
          z <- z + 1
          
          landscape <- NLMR::nlm_randomcluster(
            ncol = extent, 
            nrow = extent, 
            resolution = grain, 
            p = clust, 
            ai = c(1 - p_forest, p_forest)
          )
          landscape[landscape == 0] <- NA
          landscape[landscape == 1] <- 0
          cellids <- data.frame(
            cellnumer = 1:ncell(landscape),
            cellvalue = values(landscape)
          )
          
          pixels_to_set_nodist <- 1:ncell(landscape)
          pixels_to_set_nodist <- pixels_to_set_nodist[which(values(landscape) == 0)]
          pixels_to_set_nodist <- sample(pixels_to_set_nodist, floor(length(pixels_to_set_nodist) * p_nodist))
          
          landscape[pixels_to_set_nodist] <- -1
          
          # fun_lognormal <- function(mn, a = 0.2, b = 0.2, n) { # log-normal
          #   r <- rnorm(n)
          #   sd <- sqrt((exp(a + b * log(mn))))
          #   mmn <- log(mn) - 0.5 * log((sd/mn)^2 + 1)
          #   ssd <- sqrt(log((sd/mn)^2 + 1))
          #   rr <- r * ssd + mmn
          #   return(exp(rr))
          # }
          
          #dr_tmp <- fun_lognormal(mn = dr[k], n = 1)
          
          dr_tmp <- dr[k]
          
          n_samp_init <- floor(sum(values(landscape) == 0, na.rm = TRUE) * dr_tmp)
          
          pixels_to_sample_init <- cellids[!is.na(cellids$cellvalue), "cellnumer"]
          
          if (n_samp_init > length(pixels_to_sample_init)) {
            #n_samp <- length(pixels_to_sample)
            n_samp_init <- 0
          }
          
          landscape[sample(pixels_to_sample_init, n_samp_init)] <- 1
          
          landscape <- rast(landscape)
          
          for(i in 1:(n_year + r)) {
            
            dist <- terra::distance(landscape > 0, target = 0, exclude = NA)
            #dist <- terra::classify(dist, matrix(c(0, grain, 1, grain, 2*(extent*grain)^2, 0), ncol = 3, byrow = TRUE)) # Rook contiguity
            dist <- terra::classify(dist, matrix(c(0, sqrt(2 * grain^2), 1, sqrt(2 * grain^2), 2*(extent*grain)^2, 0), ncol = 3, byrow = TRUE)) # Queen contiguity
            
            if (is.nan(global(dist, sum, na.rm = TRUE)[[1]])) values(dist) <- 0
            
            pixels_to_sample <- 1:ncell(landscape)
            pixels_to_sample <- pixels_to_sample[which(values(landscape)[,1] == 0)]
            
            # dr_tmp <- fun_lognormal(mn = dr[k], n = 1)
            
            dr_tmp <- dr[k]
            
            n_samp <- floor(sum(!is.na(values(landscape)[, 1])) * dr_tmp)
            
            if (n_samp <= length(pixels_to_sample)) {
              
              dist_new <- sample(
                pixels_to_sample, 
                n_samp
              )
              
              landscape_updater <- landscape
              landscape_updater[landscape_updater == -1] <- 0
              landscape_updater[landscape_updater >= 1] <- 1
              
              landscape <- landscape + landscape_updater
              
              landscape[dist_new] <- 1
              
              y <- (1:(n_year + r))[i]
              
              dist_sel <- dist[dist_new][, 1]
              
              if (length(dist_new) > 0) {
                p <- c(p, mean(dist_sel))
              } else {
                p <- c(p, NA)
              }
              
            } else {
              
              #n_samp <- length(pixels_to_sample)
              n_samp <- 0
              
              dist_new <- sample(
                pixels_to_sample, 
                n_samp
              )
              
              landscape_updater <- landscape
              landscape_updater[landscape_updater == -1] <- 0
              landscape_updater[landscape_updater >= 1] <- 1
              
              landscape <- landscape + landscape_updater
              
              landscape[dist_new] <- 1
              
              y <- (1:(n_year + r))[i]
              
              dist_sel <- dist[dist_new][, 1]
              
              if (length(dist_new) > 0) {
                p <- c(p, mean(dist_sel))
              } else {
                p <- c(p, NA)
              }
              
            }
            
            dr_collector <- c(dr_collector, dr[k])
            dr_real_collector <- c(dr_real_collector, dr_tmp)
            p_forest_collector <- c(p_forest_collector, p_forest)
            y_collector <- c(y_collector, y)
            p_nodist_collector <- c(p_nodist_collector, p_nodist)
            cluster_collector <- c(cluster_collector, clust)
            
            landscape[landscape > r] <- 0
            
          }
          
          # Landscape analysis
          landscape_binary <- classify(landscape, matrix(c(NA, NA, NA, -2, 0, 0, 1, (n_year + r), 1), ncol = 3, byrow = TRUE))
          
          landscape_metrics[[z]] <- calculate_lsm(
            landscape_binary,
            level = "class",
            directions = 4
          )
          
          landscape_metrics[[z]]$dr <- dr[k]
          landscape_metrics[[z]]$p_forest <- p_forest
          landscape_metrics[[z]]$p_nodist <- p_nodist
          landscape_metrics[[z]]$cluster <- clust
          
        }
        
        
      }
      
    } 
    
  }
  
  dat <- data.frame(
    year = y_collector,
    p_forest = p_forest_collector,
    p = p,
    dr = dr_collector,
    dr_real = dr_real_collector,
    p_nodist = p_nodist_collector,
    cluster = cluster_collector
  )
  
  list(dat, landscape_metrics)
  
}

stopCluster(cl = cluster)

dat_simulations <- dat %>% 
  map(., ~.[[1]]) %>%
  bind_rows(
    .id = "iteration"
  )

write_csv(dat_simulations, "results/simulations.csv")
dat_simulations <- read_csv("results/simulations.csv")

dat_landscapemetrics <- dat %>%
  map(
    ., ~ .[[2]]
  ) %>%
  map(
    ., ~ bind_rows(.)
  ) %>%
  bind_rows(
    .id = "iteration"
  )

write_csv(dat_landscapemetrics, "results/simulations_landscapemetrics.csv")
dat_landscapemetrics <- read_csv("results/simulations_landscapemetrics.csv")

ggplot(
  data = dat_simulations %>%
    filter(
      year > r
    ) %>%
    sample_n(
      1000
    )
) +
  geom_point(
    aes(
      x = as.double(dr),
      y = as.double(dr_real)
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1
  )

A_star <- function(lambda, r, N) {
  (r * lambda * N) / (1 + r * lambda)
}

P_exact <- function(lambda, r, N) {
  A <- A_star(lambda, r, N)
  1 - (1 - 8/N)^A
}

ggplot(
  data = dat_simulations %>%
    filter(
      year > 30
    ) %>%
    group_by(
      iteration,
      p_forest,
      dr,
      p_nodist,
      cluster
    ) %>%
    summarize(
      p = mean(p, na.rm = TRUE)
    )) +
  geom_line(
    aes(
      x = dr * 100,
      y = p,
      group = interaction(p_nodist, iteration),
      col = factor(p_nodist)
    )
  ) +
  facet_grid(
    cluster ~ p_forest
  ) +
  labs(
    x = "Disturbance rate (% per yr.)",
    y = "Gap grwoth rate",
    col = "Proportion undisturbed",
    fill = "Proportion undisturbed"
  ) +
  geom_line(
    data = data.frame(
      dr0 = seq(0, 0.05, length.out = 100),
      p0 = P_exact(seq(0, 0.05, length.out = 100), r = 30, N = (1000000/30)^2)
    ),
    aes(
      x = dr0 * 100,
      y = p0
    ),
    linetype = "dashed"
  )

ggplot(
  data = dat_simulations %>%
    filter(
      year > 30
    ) %>%
    group_by(
      p_forest,
      dr,
      cluster,
      p_nodist
    ) %>%
    summarize(
      p = mean(p, na.rm = TRUE)
    ) %>%
    group_by(
      dr
    ) %>%
    summarise(
      p_mn = mean(p, na.rm = TRUE),
      p_sd = sd(p, na.rm = TRUE),
      p_md = median(p, na.rm = TRUE),
      p_lower = quantile(p, 0.25, na.rm = TRUE),
      p_upper = quantile(p, 0.75, na.rm = TRUE),
      p_llower = quantile(p, 0.025, na.rm = TRUE),
      p_uupper = quantile(p, 0.975, na.rm = TRUE)
    )
) +
  geom_ribbon(
    aes(
      x = dr * 100,
      ymin = p_llower,
      ymax = p_uupper
    ),
    alpha = 0.2
  ) +
  geom_ribbon(
    aes(
      x = dr * 100,
      ymin = p_lower,
      ymax = p_upper
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = dr * 100,
      y = p_mn
    )
  ) +
  labs(
    x = "Disturbance rate (% per yr.)",
    y = "Gap grwoth rate",
    col = "Proportion undisturbed",
    fill = "Proportion undisturbed"
  )

dat_simulations %>%
  group_by(
    dr,
    p_nodist,
    cluster
  ) %>%
  summarize() %>%
  ungroup()

ggplot(
  data = dat_landscapemetrics
) +
  geom_point(
    aes(
      x = dr,
      y = value,
      col = factor(class)
    )
  ) +
  facet_wrap(
    ~ metric,
    scales = "free"
  )

ggplot(
  data = dat_landscapemetrics %>%
    filter(
      class == 1 & dr < 0.03 & p_nodist == 0
    ) %>%
    filter(
      metric %in% c("area_mn")
    ) %>%
    group_by(
      dr,
      p_nodist,
      p_forest,
      cluster,
      metric
    ) %>%
    summarise(
      value = mean(value)
    ) %>%
    group_by(
      p_nodist,
      p_forest,
      cluster,
      metric
    ) %>%
    mutate(
      value_stan = (value / min(value))
    )
) +
  geom_point(
    aes(
      x = dr,
      y = value_stan,
      col = factor(cluster),
      shape = factor(p_forest)
    )
  ) +
  geom_smooth(
    aes(
      x = dr,
      y = value_stan,
      col = factor(cluster),
      group = interaction(cluster, p_forest)
    ),
    method = "gam",
    se = FALSE
  ) +
  scale_y_log10()

ggplot(
  data = dat_landscapemetrics %>%
    filter(
      class == 0
    ) %>%
    group_by(
      dr,
      p_nodist,
      p_forest,
      cluster,
      metric
    ) %>%
    summarize(
      value = mean(value)
    )
) +
  geom_point(
    aes(
      x = dr,
      y = value,
      col = factor(p_forest),
      shape = factor(p_nodist),
      group = factor(cluster)
    )
  ) +
  facet_wrap(
    ~ metric,
    scales = "free"
  )


