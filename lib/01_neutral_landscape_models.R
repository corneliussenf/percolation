
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

dat_simulations_plot <- dat_simulations %>%
  filter(
    iteration == 1
  ) %>%
  filter(
    p_nodist == 0,
    p_forest == 0.7,
    dr %in% sample(unique(dat_simulations$dr), 5)
  )

ggplot(
  data = dat_simulations_plot
  ) +
  geom_area(
    data = data.frame(
      x = c(0, 0, 30 ,30),
      y = c(0, 1, 1, 0)
    ),
    aes(
      x = x,
      y = y
    ),
    fill = "darkgrey",
    alpha = 0.75
  ) +
  geom_line(
    aes(
      x = year,
      y = p,
      col = dr * 100,
      group = dr
    )
  ) +
  labs(
    x = "Simulation year",
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = bquote("Disturbance rate ("*lambda*"; %"~yr^-1*")"),
    fill = NULL,
    linetype = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    strip.background = element_blank(),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_color_viridis_c() +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  geom_segment(
    data = dat_simulations_plot %>%
      filter(
        year > r
      ) %>%
      group_by(
        dr
      ) %>%
      summarise(
        pm = mean(p, na.rm = TRUE)
      ),
    aes(
      x = 30,
      xend = 160,
      y = pm,
      yend = pm,
      col = dr * 100
    ),
    alpha = 0.5
  )

ggsave("results/extended_figure04.pdf", width = 3.5, height = 3.5)

P_adj_analytical <- function(lambda, r, k = 8) {
  1 - (1 / (1 + r * lambda))^k
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
      p0 = P_adj_analytical(seq(0, 0.05, length.out = 100), r = 30)
    ),
    aes(
      x = dr0 * 100,
      y = p0
    ),
    linetype = "dashed"
  )



