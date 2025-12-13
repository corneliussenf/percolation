
library(tidyverse)
library(patchwork)

# Settings ----

recovery <- 30

# Plot theoretical solution ----

A_star <- function(lambda, r, N) {
  (r * lambda * N) / (1 + r * lambda)
}

analytical_solution_a_star <- expand.grid(
  dr0 = seq(0, 0.05, length.out = 100),
  r = c(10, 20, 30, 40)
  ) %>%
  mutate(
    disturbancearea = A_star(dr0, r = r, N = 1000) / 1000,
    forestarea = 1 - (A_star(dr0, r = r, N = 1000) / 1000)
  ) %>%
  mutate(
    p0 = P_adj_analytical(dr0, r = r, k = k),
    r_label = paste0(r, " years")
  )

p0a <- ggplot() +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate ("*lambda*"; %"~yr^-1*")"),
    y = bquote("Landscape composition (%)"),
    col = NULL,
    fill = NULL,
    linetype = NULL,
    title = bquote("Recovery period ("*r*"):")
  ) +
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
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_x_continuous(limits = c(0, 5.2)) +
  geom_line(
    data = analytical_solution_a_star,
    aes(
      x = dr0 * 100,
      y = disturbancearea * 100,
      col = "Open forests"
    )
  ) +
  geom_line(
    data = analytical_solution_a_star,
    aes(
      x = dr0 * 100,
      y = forestarea * 100,
      col = "Closed forest"
    )
  ) +
  scale_color_manual(
    values = c("#810f7c", "#006d2c")
  ) +
  facet_wrap(
    ~r_label,
    ncol = 4,
    scales = "free"
  )

ggsave(plot = p0a, "results/steady_state.pdf", width = 7.5, height = 2.5)

P_adj_analytical <- function(lambda, r, k = 8) {
  1 - (1 / (1 + r * lambda))^k
}

analytical_solution_plot <- expand.grid(
  dr0 = seq(0, 0.05, length.out = 100), 
  k = 1:8,
  r = c(10, 20, 30, 40)
  ) %>%
  mutate(
    p0 = P_adj_analytical(dr0, r = r, k = k),
    r_label = paste0(r, " years")
  )

p0b <- ggplot() +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate ("*lambda*"; %"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = bquote("Number of neighbors ("*k*")"),
    fill = NULL,
    linetype = NULL,
    title = bquote("Recovery period ("*r*"):")
  ) +
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
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  geom_line(
    data = analytical_solution_plot,
    aes(
      x = dr0 * 100,
      y = p0,
      col = k,
      group = k
    )
  ) +
  scale_color_viridis_c(
    breaks = c(1, 4, 7)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  facet_wrap(
    ~r_label,
    ncol = 4,
    scales = "free"
  )

ggsave(plot = p0b, "results/theoretical_solution_full.pdf", width = 7.5, height = 2.5)

p0a + labs(x = NULL) +
  p0b + labs(title = NULL) + theme(strip.text = element_blank()) +
  plot_layout(
    ncol = 1
    ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

ggsave("results/extended_figure01.pdf", width = 7.5, height = 5)

# Get data ----

dat_full_summary <- read_csv("data/dat_full_summary.csv")
dat_natural_full_summary <- read_csv("data/dat_natural_full_summary.csv")
dat_harvest_full_summary <- read_csv("data/dat_harvest_full_summary.csv")
disturbance_summary <- read_csv("data/disturbance_summary_1985_2023.csv")
disturbance_summary_change <- read_csv("data/disturbance_summary_change.csv")
forest_matrix <- read_csv("data/forest_matrix.csv")
simulations <- read_csv("results/simulations.csv")
country_intersection <- read_csv("data/gridcell_country.csv")

# Calculate gap growth probabilities and plot ----

## Calulcate gap growth probability

expansion_rates <- dat_full_summary %>%
  pivot_wider(
    names_from = "expansion",
    values_from = c("n", "area_ha_t0", "area_ha_t1")
  ) %>%
  mutate(
    expansion_rate_n = n_1 / (n_0 + n_1),
    expansion_rate_area = area_ha_t1_1 / (area_ha_t1_0 + area_ha_t1_1),
    total_dist = (area_ha_t1_0 + area_ha_t1_1)
  ) %>%
  filter(
    year > (1986 + recovery)
  ) %>%
  group_by(
    gridid
  ) %>%
  summarize(
    expansion_rate_n_mn = mean(expansion_rate_n, na.rm = TRUE),
    expansion_rate_area_mn = mean(expansion_rate_area, na.rm = TRUE),
    expansion_rate_n_se = sd(expansion_rate_n, na.rm = TRUE) / sqrt(sum(!is.na(expansion_rate_n))),
    expansion_rate_area_se = sd(expansion_rate_area, na.rm = TRUE) / sqrt(sum(!is.na(expansion_rate_area)))
  ) %>%
  left_join(
    disturbance_summary,
    by = "gridid"
  )

expansion_rates_natural <- dat_natural_full_summary %>%
  pivot_wider(
    names_from = "expansion",
    values_from = c("n", "area_ha_t0", "area_ha_t1")
  ) %>%
  mutate(
    expansion_rate_n = n_1 / (n_0 + n_1),
    expansion_rate_area = area_ha_t1_1 / (area_ha_t1_0 + area_ha_t1_1),
    total_dist_rate = (area_ha_t1_0 + area_ha_t1_1)
  ) %>%
  filter(
    year > (1986 + recovery)
  ) %>%
  group_by(
    gridid
  ) %>%
  summarize(
    expansion_rate_n_mn = mean(expansion_rate_n, na.rm = TRUE),
    expansion_rate_area_mn = mean(expansion_rate_area, na.rm = TRUE),
    expansion_rate_n_se = sd(expansion_rate_n, na.rm = TRUE) / sqrt(sum(!is.na(expansion_rate_n))),
    expansion_rate_area_se = sd(expansion_rate_area, na.rm = TRUE) / sqrt(sum(!is.na(expansion_rate_area)))
  ) %>%
  left_join(
    disturbance_summary,
    by = "gridid"
  )

expansion_rates_harvest <- dat_harvest_full_summary %>%
  pivot_wider(
    names_from = "expansion",
    values_from = c("n", "area_ha_t0", "area_ha_t1")
  ) %>%
  mutate(
    expansion_rate_n = n_1 / (n_0 + n_1),
    expansion_rate_area = area_ha_t1_1 / (area_ha_t1_0 + area_ha_t1_1),
    total_dist_rate = (area_ha_t1_0 + area_ha_t1_1)
  ) %>%
  filter(
    year > (1986 + recovery)
  ) %>%
  group_by(
    gridid
  ) %>%
  summarize(
    expansion_rate_n_mn = mean(expansion_rate_n, na.rm = TRUE),
    expansion_rate_area_mn = mean(expansion_rate_area, na.rm = TRUE),
    expansion_rate_n_se = sd(expansion_rate_n, na.rm = TRUE) / sqrt(sum(!is.na(expansion_rate_n))),
    expansion_rate_area_se = sd(expansion_rate_area, na.rm = TRUE) / sqrt(sum(!is.na(expansion_rate_area)))
  ) %>%
  left_join(
    disturbance_summary,
    by = "gridid"
  )

## Combine and save

expansion_rates_all <- list(
  expansion_rates,
  expansion_rates_natural,
  expansion_rates_harvest
) %>%
  set_names(
    c("total", "natural", "harvest")
  ) %>%
  bind_rows(
    .id = "subset"
  )

write_csv(
  expansion_rates_all,
  "results/gap_growth_probabilities.csv"
)

## Plot distribution

p1 <- ggplot() +
  geom_histogram(
    data = expansion_rates_natural,
    aes(
      x = expansion_rate_n_mn,
      y = ..density..,
      fill = "natural"
    ),
    alpha = 0.5
  ) +
  geom_histogram(
    data = expansion_rates_harvest,
    aes(
      x = expansion_rate_n_mn,
      y = ..density..,
      fill = "harvest"
    ),
    alpha = 0.5
  ) +
  theme_classic() +
  labs(
    x = bquote("Gap expansion probability ("*p[g]*")"),
    y = "Density",
    fill = NULL,
    title = "Empirical data"
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_fill_manual(
    values = c("#E4632D", "#8CBC68"),
    labels = c("Harvest", "Natural disturbances")
  )
 
## Model gap growth probabilities over disturbance rate

# exp_model <- function(data, p_name = "p", dr_name = "dr") {
#   data$p <- data[, p_name][[1]]
#   data$dr <- data[, dr_name][[1]]
#   s.0 <- 1
#   k.0 <- 15
#   start <- list(k = k.0, s = s.0)
#   model <- nls(p ~ s - s * exp(-k * dr), data = data, start = start, weights = forest_px_sum)
#   return(model)
# }

analytical_model <- function(data, p_name = "p", dr_name = "dr") {
  data$p <- data[, p_name][[1]]
  data$dr <- data[, dr_name][[1]]
  r <- 30
  k.0 <- 8
  start <- list(k = k.0)
  model <- nls(p ~ 1 - (1 / (1 + r * dr))^k, data = data, start = start)
  return(model)
}

# model_emp_natural <- exp_model(
#   expansion_rates_natural %>% na.omit(), 
#   p = "expansion_rate_n_mn", 
#   dr = "rate")

model_emp_natural <- analytical_model(
  expansion_rates_natural %>% na.omit(), 
  p = "expansion_rate_n_mn", 
  dr = "rate")

# model_emp_harvest <- exp_model(
#   expansion_rates_harvest %>% na.omit(), 
#   p = "expansion_rate_n_mn",
#   dr = "rate")

model_emp_harvest <- analytical_model(
  expansion_rates_harvest %>% na.omit(), 
  p = "expansion_rate_n_mn",
  dr = "rate")

k <- coef(model_emp_natural)["k"]
r <- 30
#s <- coef(model_emp_natural)["s"]
pred_emp_natural <- nlraa::predict2_nls(model_emp_natural, 
                                        newdata = data.frame(dr = seq(0, max(expansion_rates_all$rate, na.rm = TRUE), 0.001)),
                                        interval = "conf"
)
pred_emp_natural$dr <- seq(0, max(expansion_rates_all$rate, na.rm = TRUE), 0.001)

k <- coef(model_emp_harvest)["k"]
r <- 30
#s <- coef(model_emp_harvest)["s"]
pred_emp_harvest <- nlraa::predict2_nls(model_emp_harvest, 
                                        newdata = data.frame(dr = seq(0, max(expansion_rates_all$rate, na.rm = TRUE), 0.001)),
                                        interval = "conf"
)
pred_emp_harvest$dr <- seq(0, max(expansion_rates_all$rate, na.rm = TRUE), 0.001)

pred_emp <- list(
  pred_emp_natural,
  pred_emp_harvest
) %>%
  bind_rows()

simulations_summary <- simulations %>%
  filter(
    year > recovery
  ) %>%
  filter(
    p_nodist == 0
  ) %>%
  group_by(
    dr,
    cluster,
    p_forest
  ) %>%
  summarize(
    mn = mean(p, na.rm = TRUE),
    lower = quantile(p, 0.25, na.rm = TRUE),
    upper = quantile(p, 0.75, na.rm = TRUE)
  ) %>%
  ungroup()

simulations_summary <- list(
  simulations_summary,
  expand_grid(
    dr = 0,
    cluster = c(0.1, 0.3, 0.5), 
    p_forest = c(0.1, 0.3, 0.5, 0.7, 0.9),
    mn = 0, 
    lower = 0, 
    upper = 0
  )
  ) %>%
  bind_rows()

simulations_summary_p_nodist <- simulations %>%
  filter(
    year > recovery
  ) %>%
  filter(
    cluster == 0.3 & p_forest == 0.3
  ) %>%
  group_by(
    dr,
    iteration,
    p_nodist
  ) %>%
  summarize(
    p_mn = mean(p, na.rm = TRUE)
  ) %>%
  ungroup()

simulations_summary_p_nodist <- list(
  simulations_summary_p_nodist,
  expand_grid(
    dr = 0,
    p_nodist = c(0, 0.2, 0.4, 0.6),
    p_mn = 0
  )
  ) %>%
  bind_rows()

simulations_summary_p_forest <- simulations %>%
  filter(
    year > recovery
  ) %>%
  filter(
    cluster == 0.3 & p_nodist == 0
  ) %>%
  group_by(
    dr,
    iteration,
    p_forest
  ) %>%
  summarize(
    p_mn = mean(p, na.rm = TRUE)
  ) %>%
  ungroup()

simulations_summary_p_forest <- list(
  simulations_summary_p_forest,
  expand_grid(
    dr = 0,
    p_forest = c(0.1, 0.3, 0.5, 0.7, 0.9),
    p_mn = 0
  )
) %>%
  bind_rows()

simulations_summary_clustering <- simulations %>%
  filter(
    year > recovery
  ) %>%
  filter(
    p_forest == 0.3 & p_nodist == 0
  ) %>%
  group_by(
    dr,
    iteration,
    cluster
  ) %>%
  summarize(
    p_mn = mean(p, na.rm = TRUE)
  ) %>%
  ungroup()

simulations_summary_clustering <- list(
  simulations_summary_clustering,
  expand_grid(
    dr = 0,
    cluster = c(0.1, 0.3, 0.5),
    p_mn = 0
  )
) %>%
  bind_rows()

P_adj_analytical <- function(lambda, r, k = 8) {
  1 - (1 / (1 + r * lambda))^k
}

analytical_solution <- expand.grid(
  dr0 = seq(0, 0.05, length.out = 100), 
  k = 1:8
  ) %>%
  mutate(
    p0 = P_adj_analytical(dr0, r = 30, k = k)
  )

p2b <- ggplot() +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = bquote("Number of neighbors ("*k*")"),
    fill = NULL,
    linetype = NULL,
    title = "Theoretical solution"
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      col = k,
      group = k
    )
  ) +
  scale_color_viridis_c(
    breaks = c(1, 4, 7, 10)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

p3 <- ggplot() +
  geom_ribbon(
    data = simulations_summary,
    aes(
      x = dr * 100,
      ymin = lower,
      ymax = upper,
      group = interaction(cluster, p_forest)
    ),
    fill = "#33BBEE",
    alpha = 0.1
  ) +
  geom_line(
    data = simulations_summary,
    aes(
      x = dr * 100,
      y = mn,
      group = interaction(cluster, p_forest),
      col = "Simulations"
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability"),
    col = NULL,
    fill = NULL,
    linetype = NULL,
    title = "Simulations"
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      group = factor(k),
      col = "Analytical solution"
    ),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("#BBBBBB", "#33BBEE")
  )

p4a <- ggplot() +
  geom_line(
    data = simulations_summary_p_nodist,
    aes(
      x = dr * 100,
      y = p_mn,
      col = factor(p_nodist * 100),
      group = interaction(p_nodist, iteration)
    ),
    alpha = 0.4
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = "Excluded (%)",
    fill = NULL,
    linetype = NULL,
    title = "  "
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      group = factor(k),
      col = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_color_manual(
    values = c("#b3cde3", "#8c96c6", "#8856a7", "#810f7c", "#BBBBBB")
  )

p4b <- ggplot() +
  geom_line(
    data = simulations_summary_p_forest,
    aes(
      x = dr * 100,
      y = p_mn,
      col = factor(p_forest * 100),
      group = interaction(p_forest, iteration)
    ),
    alpha = 0.4
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = "Forest cover (%)",
    fill = NULL,
    linetype = NULL,
    title = "  "
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  scale_linetype_manual(values = c(2, 1)) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      group = factor(k),
      col = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824", "#BBBBBB")
  )

p4c <- ggplot() +
  geom_line(
    data = simulations_summary_clustering,
    aes(
      x = dr * 100,
      y = p_mn,
      col = factor(cluster),
      group = interaction(cluster, iteration)
    ),
    alpha = 0.4
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = "Clustering",
    fill = NULL,
    linetype = NULL,
    title = "  "
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.25, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  scale_linetype_manual(values = c(2, 1)) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      group = factor(k),
      col = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("#7fcdbb", "#41b6c4", "#1d91c0", "#BBBBBB")
  )

p4a + p4b + p4c

p5 <- ggplot() +
  geom_point(
    data = expansion_rates_harvest,
    aes(
      x = rate * 100,
      y = expansion_rate_n_mn,
      size = forest_px_sum,
      col = "Observed"
    ) 
  ) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      group = factor(k),
      col = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = pred_emp_harvest,
    aes(
      x = dr * 100,
      ymin = Q2.5,
      ymax = Q97.5
    ),
    alpha = 0.5
  ) +
  geom_line(
    data = pred_emp_harvest,
    aes(
      x = dr * 100,
      y = Estimate
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = NULL,
    fill = NULL,
    linetype = NULL,
    title = " "
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.5, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  guides(size = "none") +
  scale_color_manual(values = c("#E4632D", "#BBBBBB")) +
  scale_linetype_manual(values = c(2, 1))

p6 <- ggplot() +
  geom_point(
    data = expansion_rates_natural,
    aes(
      x = rate * 100,
      y = expansion_rate_n_mn,
      # col = "Natural disturbance",
      size = forest_px_sum,
      col = "Observed"
    )
  ) +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0,
      group = factor(k),
      col = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = pred_emp_natural,
    aes(
      x = dr * 100,
      ymin = Q2.5,
      ymax = Q97.5
    ),
    alpha = 0.5
  ) +
  geom_line(
    data = pred_emp_natural,
    aes(
      x = dr * 100,
      y = Estimate
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap expansion probability ("*p[g]*")"),
    col = NULL,
    fill = NULL,
    linetype = NULL,
    title = " "
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing = unit(-0.5, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 5.2)) +
  geom_text(
    data = analytical_solution  %>%
      filter(
        dr0 == max(dr0)
      ) %>%
      filter(
        k == 8
      ),
    aes(
      x = dr0 * 100,
      y = p0 * 100,
      group = factor(k),
      label = paste0("k = ", k)
    ),
    size = 1.5,
    hjust = -0.1
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  scale_color_manual(values = c("#8CBC68", "#BBBBBB")) +
  guides(size = "none") +
  scale_linetype_manual(values = c(2, 1))

p2b + 
  p4b + labs(title = "Simulations") +
  p4a +
  p1 + coord_flip() + theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
  p6 + p5 +
  plot_layout(
    ncol = 3
  ) +
  plot_annotation(
    tag_levels = c("a"),
    tag_suffix = ")"
  )

ggsave("results/figure01.pdf", width = 7.5, height = 7.5 / 3 * 2)

# By country

ggplot() +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0 * 100,
      group = factor(k),
      color = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  geom_point(
    data = expansion_rates_harvest %>%
      left_join(
        country_intersection
      ) %>%
      filter(
        !is.na(name0)
      ) %>%
      mutate(
        name0 = case_when(
          name0 == "U.K. of Great Britain and Northern Ireland" ~ "United Kingdom",
          name0 == "Moldova, Republic of" ~ "Moldova",
          name0 == "The former Yugoslav Republic of Macedonia" ~ "Northern Macedonia",
          name0 == "Czech Republic" ~ "Czechia",
          TRUE ~ name0
        )
      ),
    aes(
      x = rate * 100,
      y = expansion_rate_n_mn * 100,
      col = "Harvest"
    )
  ) +
  geom_point(
    data = expansion_rates_natural %>%
      left_join(
        country_intersection
      ) %>%
      filter(
        !is.na(name0)
      ) %>%
      mutate(
        name0 = case_when(
          name0 == "U.K. of Great Britain and Northern Ireland" ~ "United Kingdom",
          name0 == "Moldova, Republic of" ~ "Moldova",
          name0 == "The former Yugoslav Republic of Macedonia" ~ "Northern Macedonia",
          name0 == "Czech Republic" ~ "Czechia",
          TRUE ~ name0
        )
      ),
    aes(
      x = rate * 100,
      y = expansion_rate_n_mn * 100,
      col = "Natural disturbances"
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap growth probability ("*p[g]*")"),
    col = NULL
  ) +
  scale_color_manual(
    values = c("#E4632D", "#8CBC68", "#BBBBBB")
  ) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 6)
  ) +
  scale_y_continuous(limits = c(0, 1.08) * 100, 
                     breaks = seq(0, 1, 0.2) * 100) +
  scale_x_continuous(limits = c(0, 5.2)) +
  facet_wrap(
    ~name0,
    scales = "free"
  )

ggsave("results/extended_figure02.pdf", width = 8.5, height = 8.5)

ggplot() +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0 * 100,
      group = factor(k),
      color = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  geom_point(
    data = expansion_rates_harvest %>%
      left_join(
        country_intersection
      ) %>%
      filter(
        !is.na(name0)
      ),
    aes(
      x = rate * 100,
      y = expansion_rate_n_mn * 100,
      col = "Harvest"
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap growth probability"),
    col = NULL
  ) +
  scale_color_manual(
    values = c("#E4632D", "#BBBBBB")
  ) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 6)
  ) +
  scale_y_continuous(limits = c(0, 1.08) * 100, 
                     breaks = seq(0, 1, 0.2) * 100) +
  scale_x_continuous(limits = c(0, 5.2)) +
  facet_wrap(
    ~name0
  )

ggsave("results/gap_growth_rate_harvest_by_country.pdf", width = 7.5, height = 7.5)

ggplot() +
  geom_line(
    data = analytical_solution,
    aes(
      x = dr0 * 100,
      y = p0 * 100,
      group = factor(k),
      color = "Theoretical solution"
    ),
    linetype = "dashed"
  ) +
  geom_point(
    data = expansion_rates_natural %>%
      left_join(
        country_intersection
      ) %>%
      filter(
        !is.na(name0)
      ),
    aes(
      x = rate * 100,
      y = expansion_rate_n_mn * 100,
      col = "Natural disturbance"
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = bquote("Gap growth probability"),
    col = NULL
  ) +
  scale_color_manual(
    values = c("#8CBC68", "#BBBBBB"),
    labels = c("Natural disturbance")
  ) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 6)
  ) +
  scale_y_continuous(limits = c(0, 1.08) * 100, 
                     breaks = seq(0, 1, 0.2) * 100) +
  scale_x_continuous(limits = c(0, 5.2)) +
  facet_wrap(~name0) +
  scale_linetype_manual(values = c(2, 1))

ggsave("results/gap_growth_rate_natural_by_country.pdf", width = 7.5, height = 7.5)

# Difference in gap size new versus expanding ---- 

exp_model <- function(data, y_name, dr_name = "rate") {
  data$y <- data[, y_name][[1]]
  data$rate <- data[, dr_name][[1]]
  start <- list(a = min(data$y), b = 250)
  model <- nls(y ~ I(a * exp(b * rate)), 
               data = data, 
               start = start, 
               weights = forest_px_sum
  )
  return(model)
}

logit_model <- function(data, y_name, dr_name = "rate") {
  data$y <- data[, y_name][[1]]
  data$rate <- data[, dr_name][[1]]
  #start <- list(a = min(data$y), b = 250, k = 0.5)
  model <- nls(y ~ SSlogis(rate, a, b, c), 
               data = data, 
               #start = start, 
               weights = forest_px_sum
  )
  return(model)
}

# Harvest disturbances

dat_harvest_summary_subset_plot <- dat_harvest_full_summary %>% 
  left_join(
    disturbance_summary
  ) %>%
  mutate(
    forest_share = forest_px_sum * 9e-04 / 100^2
  ) %>%
  filter(
    forest_share >= 0.1 
  ) %>%
  filter(
    rate > 0
  ) %>%
  group_by(
    expansion,
    year, 
    forest_px_sum,
    gridid,
    rate
  ) %>%
  summarize(
    n = sum(n),
    area_ha_t0 = sum(area_ha_t0),
    area_ha_t1 = sum(area_ha_t1)
  ) %>%
  mutate(
    patch_size_mn = (area_ha_t0 + area_ha_t1) / n
  ) %>%
  group_by(
    expansion,
    gridid, 
    forest_px_sum
  ) %>%
  summarize(
    rate = unique(rate),
    patch_size_mn = mean(patch_size_mn)
  ) %>%
  ungroup() %>%
  left_join(
    expansion_rates_harvest
  )

lin.model <- lm(
  patch_size_mn ~ rate,
  data = dat_harvest_summary_subset_plot %>%
    filter(
      expansion == 1
    ),
  weights = forest_px_sum
)

exp.model <- exp_model(
  dat_harvest_summary_subset_plot %>%
    filter(
      expansion == 1
    ), 
  y_name = "patch_size_mn"
)

logit.model <- logit_model(
  dat_harvest_summary_subset_plot %>%
    filter(
      expansion == 1
    ), 
  y_name = "patch_size_mn")

AIC(lin.model, exp.model, logit.model)

broom::tidy(exp.model)
broom::glance(exp.model)

broom::tidy(logit.model)
broom::glance(logit.model)

broom::tidy(lin.model)
broom::glance(lin.model)

saturation_point_harvest <- coef(logit.model)["a"]

final_model <- logit.model

model.df.harvest <- expand_grid(
  rate = seq(min(dat_harvest_summary_subset_plot$rate), 0.05, length.out = 100)
)

model.df.harvest <- cbind(
  rate = model.df.harvest$rate,
  nlraa::predict2_nls(final_model, 
                      newdata = model.df.harvest,
                      interval = "conf"
  )
)

# Natural disturbances

dat_natural_summary_subset_plot <- dat_natural_full_summary %>% 
  left_join(
    disturbance_summary
  ) %>%
  mutate(
    forest_share = forest_px_sum * 9e-04 / 100^2
  ) %>%
  filter(
    forest_share >= 0.1 
  ) %>%
  filter(
    rate > 0
  ) %>%
  group_by(
    expansion,
    year, 
    forest_px_sum,
    gridid,
    rate
  ) %>%
  summarize(
    n = sum(n),
    area_ha_t0 = sum(area_ha_t0),
    area_ha_t1 = sum(area_ha_t1)
  ) %>%
  mutate(
    patch_size_mn = (area_ha_t0 + area_ha_t1) / n
  ) %>%
  group_by(
    expansion,
    gridid, 
    forest_px_sum
  ) %>%
  summarize(
    rate = unique(rate),
    patch_size_mn = mean(patch_size_mn)
  ) %>%
  ungroup() %>%
  left_join(
    expansion_rates_harvest
  )

lin.model <- lm(
  patch_size_mn ~ rate,
  data = dat_natural_summary_subset_plot %>%
    filter(
      expansion == 1
    ),
  weights = forest_px_sum
)

exp.model <- exp_model(
  dat_natural_summary_subset_plot %>%
    filter(
      expansion == 1
    ), 
  y_name = "patch_size_mn")

logit.model <- logit_model(
  dat_natural_summary_subset_plot %>%
    filter(
      expansion == 1
    ), 
  y_name = "patch_size_mn")

AIC(lin.model, exp.model, logit.model)

broom::tidy(exp.model)
broom::glance(exp.model)

broom::tidy(logit.model)
broom::glance(logit.model)

broom::tidy(lin.model)
broom::glance(lin.model)

final_model <- exp.model

a <- coef(final_model)["a"]
b <- coef(final_model)["b"]
model.df.natural <- expand_grid(
  rate = seq(min(dat_natural_summary_subset_plot$rate), 0.05, length.out = 100)
)

model.df.natural <- cbind(
  rate = model.df.natural$rate,
  nlraa::predict2_nls(final_model, 
                      newdata = model.df.natural,
                      interval = "conf"
  )
)

# Joint plot

p5 <- ggplot() +
  geom_point(
    data = dat_harvest_summary_subset_plot %>% filter(expansion == 1),
    aes(
      x = rate * 100,
      y = patch_size_mn,
      size = forest_px_sum,
      shape = "Expanding gap"
    ),
    col = "#E4632D",
    #shape = 19
  ) +
  geom_point(
    data = dat_harvest_summary_subset_plot %>% filter(expansion == 0),
    aes(
      x = rate * 100,
      y = patch_size_mn,
      size = forest_px_sum,
      shape = "New gap"
    ),
    col = "#E4632D",
    #shape = 21
  ) +
  geom_ribbon(
    data = model.df.harvest,
    aes(
      x = rate * 100,
      ymin = Q2.5,
      ymax = Q97.5
      #group = factor(expansion)
    ),
    fill = "#000000",
    alpha = 0.25
  ) +
  geom_line(
    data = model.df.harvest,
    aes(
      x = rate * 100,
      y = Estimate,
      #group = factor(expansion)
    ),
    col = "#000000"
  ) +
  scale_y_log10(
    limits = c(0.5, 4800)
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = "Average gap size (ha)",
    col = NULL,
    title = "Harvest",
    shape = NULL
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm")
  ) +
  guides(
    size = "none",
    shape = guide_legend(override.aes = list(col = "#000000"))
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  scale_shape_manual(
    values = c(19, 21)
  ) +
  geom_hline(
    yintercept = saturation_point_harvest,
    linetype = "dashed"
  )

p6 <- ggplot() +
  geom_point(
    data = dat_natural_summary_subset_plot %>% filter(expansion == 1),
    aes(
      x = rate * 100,
      y = patch_size_mn,
      size = forest_px_sum,
      shape = "Expanding gap"
    ),
    col = "#8CBC68",
    #shape = 19
  ) +
  geom_point(
    data = dat_natural_summary_subset_plot %>% filter(expansion == 0),
    aes(
      x = rate * 100,
      y = patch_size_mn,
      size = forest_px_sum,
      shape = "New gap"
    ),
    col = "#8CBC68",
    #shape = 21
  ) +
  geom_ribbon(
    data = model.df.natural,
    aes(
      x = rate * 100,
      ymin = Q2.5,
      ymax = Q97.5
      #group = factor(expansion)
    ),
    fill = "#000000",
    alpha = 0.25
  ) +
  geom_line(
    data = model.df.natural,
    aes(
      x = rate * 100,
      y = Estimate,
      # group = factor(expansion)
    ),
    col = "#000000"
  ) +
  scale_y_log10(
    limits = c(1, 4800)
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = "Average gap size (ha)",
    col = NULL,
    title = "Natural disturbances",
    shape = NULL
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm")
  ) +
  guides(
    size = "none",
    shape = guide_legend(override.aes = list(col = "#000000"))
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  scale_shape_manual(
    values = c(19, 21)
  )

p7 <- ggplot() +
  geom_histogram(
    data = dat_natural_summary_subset_plot,
    aes(
      x = patch_size_mn,
      fill = "Natural",
      alpha = ifelse(expansion == 1, "Expanding gaps", "New gaps")
    ),
    binwidth = 0.11,
    position = "identity"
  ) +
  geom_text(
    data = dat_natural_summary_subset_plot %>% group_by(expansion) %>% summarise(med = median(patch_size_mn)),
    aes(
      x = med,
      y = c(125, 102),
      label = paste(round(med, 0), "ha")
    ),
    size = 2,
    inherit.aes = FALSE,
    check_overlap = TRUE
  ) +
  geom_histogram(
    data = dat_harvest_summary_subset_plot,
    aes(
      x = patch_size_mn,
      fill = "Harvest",
      alpha = ifelse(expansion == 1, "Growing gaps", "New gaps")
    ),
    binwidth = 0.11,
    position = "identity"
  ) +
  geom_text(
    data = dat_harvest_summary_subset_plot %>% group_by(expansion) %>% summarise(med = median(patch_size_mn)),
    aes(
      x = med,
      y = c(215, 110),
      label = paste(round(med, 0), "ha")
    ),
    size = 2,
    inherit.aes = FALSE,
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = c("Natural" = "#8CBC68", "Harvest" = "#E4632D")
  ) +
  scale_alpha_manual(
    values = c("Expanding gaps" = 1, "New gaps" = 0.4)
  ) +
  scale_x_log10() +
  theme_classic() +
  labs(
    x = "Average gap size (ha)",
    y = "Count",
    title = " ",
    fill = NULL,
    alpha = NULL
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(1, 1, 1, 1)
  )

p7 + 
  p5 +
  p6 +
  plot_layout(
    ncol = 3
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

ggsave("results/figure02.pdf", width = 7.3, height = 7.5 / 3)

# Forest matrix analyses ----

plotdat <- forest_matrix %>%
  mutate(
    forest_area_open_closed = cut(
      forest_px_sum, 
      c(0, 0.4 * 100^2 / 9e-04, max(forest_px_sum)),
      include.lowest = TRUE,
      labels = c("Open", "Closed")
    ),
    forest_area_cut = cut(
      forest_px_sum,
      # c(c(0, 0.2, 0.4, 0.6, 0.8) * 100^2 / 9e-04, max(forest_px_sum)),
      c(c(0, 0.2, 0.6) * 100^2 / 9e-04, max(forest_px_sum)),
      include.lowest = TRUE,
      labels = c("0-20%", "20-60%", ">60%"),
      # labels = c("0-20%", "20-40%", "40-60%", "60-80%", ">80%")
    ),
    forest_share = forest_px_sum * 9e-04 / 100^2
  )

plotdat_nona <- plotdat %>% 
  filter(forest_share >= 0.1) %>%
  filter(rate > 0)

# Size
fit_area.log1 <- lm(log(value) ~ (rate) * forest_area_cut, 
                    data = plotdat_nona %>% filter(name %in% c("area_mn")),
                    weights = forest_px_sum)
fit_area.log2 <- lm(log(value) ~ log(rate) * forest_area_cut, 
                    data = plotdat_nona %>% filter(name %in% c("area_mn")),
                    weights = forest_px_sum)

AIC(fit_area.log1, fit_area.log2)

fit_area.log.df <- expand_grid(
  rate = seq(0.002, 0.05, length.out = 100),
  name = "area",
  #forest_area_cut = "20-60%"
  forest_area_cut = unique(plotdat_nona$forest_area_cut)
  #forest_px_sum = mean(plotdat_nona$forest_px_sum)
)
fit_area.log.df$pred = predict(fit_area.log1, newdata = fit_area.log.df, se.fit = TRUE, interval = "prediction")$fit[, 1]
fit_area.log.df$pred_lower = predict(fit_area.log1, newdata = fit_area.log.df, se.fit = TRUE, interval = "confidence")$fit[, 2]
fit_area.log.df$pred_upper = predict(fit_area.log1, newdata = fit_area.log.df, se.fit = TRUE, interval = "confidence")$fit[, 3]

# Number of patches
fit_np.log1 <- lm(log(value / 100^2) ~ (rate) * forest_area_cut, 
                  data = plotdat_nona %>% filter(name %in% c("np")),
                  weights = forest_px_sum)
fit_np.log2 <- lm(log(value / 100^2) ~ log(rate) * forest_area_cut, 
                  data = plotdat_nona %>% filter(name %in% c("np")),
                  weights = forest_px_sum)

AIC(fit_np.log1, fit_np.log2)

fit_np.log.df <- expand_grid(
  rate = seq(0.002, 0.05, length.out = 100),
  name = "np",
  #forest_area_cut = "20-60%"
  forest_area_cut = unique(plotdat_nona$forest_area_cut)
  #forest_px_sum = mean(plotdat_nona$forest_px_sum)
)
fit_np.log.df$pred = predict(fit_np.log1, newdata = fit_np.log.df, se.fit = TRUE, interval = "confidence")$fit[, 1]
fit_np.log.df$pred_lower = predict(fit_np.log1, newdata = fit_np.log.df, se.fit = TRUE, interval = "confidence")$fit[, 2]
fit_np.log.df$pred_upper = predict(fit_np.log1, newdata = fit_np.log.df, se.fit = TRUE, interval = "confidence")$fit[, 3]

# Aggregation
fit_ai.log1 <- lm(log(value) ~ (rate) * forest_area_cut, 
                  data = plotdat_nona %>% filter(name %in% c("ai")),
                  weights = forest_px_sum)
fit_ai.log2 <- lm(log(value) ~ log(rate) * forest_area_cut, 
                  data = plotdat_nona %>% filter(name %in% c("ai")),
                  weights = forest_px_sum)

AIC(fit_ai.log1, fit_ai.log2)

fit_ai.log.df <- expand_grid(
  rate = seq(0.002, 0.05, length.out = 100),
  name = "ai",
  forest_area_cut = unique(plotdat_nona$forest_area_cut)
  #forest_px_sum = mean(plotdat_nona$forest_px_sum)
)
fit_ai.log.df$pred = predict(fit_ai.log1, newdata = fit_np.log.df, se.fit = TRUE, interval = "confidence")$fit[, 1]
fit_ai.log.df$pred_lower = predict(fit_ai.log1, newdata = fit_np.log.df, se.fit = TRUE, interval = "confidence")$fit[, 2]
fit_ai.log.df$pred_upper = predict(fit_ai.log1, newdata = fit_np.log.df, se.fit = TRUE, interval = "confidence")$fit[, 3]

# Plots
p_area <- ggplot() +
  geom_point(
    data = plotdat_nona  %>%
      filter(
        name == "area_mn"
      ),
    aes(
      x = rate * 100,
      y = value,
      size = forest_px_sum,
      col = forest_area_cut
    )
  ) +
  geom_ribbon(
    data = fit_area.log.df,
    aes(
      x = (rate) * 100,
      ymin = exp(pred_lower),
      ymax = exp(pred_upper),
      group = forest_area_cut
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = fit_area.log.df,
    aes(
      x = (rate) * 100,
      y = exp(pred),
      group = forest_area_cut
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = "Average patch size (ha)",
    #y = "Hectare",
    col = "Forest cover"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#41ae76", "#005824")
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  ylim(0, 45) +
  guides(
    size = "none"
  ) +
  facet_wrap(
    ~forest_area_cut,
    scales = "free"
  )

p_area_clean <- ggplot() +
  geom_point(
    data = plotdat_nona  %>%
      filter(
        name == "area_mn"
      ),
    aes(
      x = rate * 100,
      y = value,
      size = forest_px_sum,
      col = forest_area_cut
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = "Average patch size (ha)",
    #y = "Hectare",
    col = "Forest cover"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#41ae76", "#005824")
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  ylim(0, 45) +
  guides(
    size = "none"
  )

p_np <- ggplot() +
  geom_point(
    data = plotdat_nona  %>%
      filter(
        name == "np"
      ),
    aes(
      x = rate * 100,
      y = value / 100^2,
      size = forest_px_sum,
      col = forest_area_cut
    )
  ) +
  geom_ribbon(
    data = fit_np.log.df,
    aes(
      x = (rate) * 100,
      ymin = exp(pred_lower),
      ymax = exp(pred_upper),
      group = forest_area_cut
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = fit_np.log.df,
    aes(
      x = (rate) * 100,
      y = exp(pred),
      group = forest_area_cut
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    #title = "Number of patches",
    y = bquote("Number of patches (100"~km^-2*")"),
    col = "Forest cover"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#41ae76", "#005824")
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  ylim(0, 18) +
  guides(
    size = "none"
  ) +
  facet_wrap(
    ~forest_area_cut,
    scales = "free"
  )

p_np_clean <- ggplot() +
  geom_point(
    data = plotdat_nona  %>%
      filter(
        name == "np"
      ),
    aes(
      x = rate * 100,
      y = value / 100^2,
      size = forest_px_sum,
      col = forest_area_cut
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    #title = "Number of patches",
    y = bquote("Number of patches (100"~km^-2*")"),
    col = "Forest cover"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#41ae76", "#005824")
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  ylim(0, 18) +
  guides(
    size = "none"
  )

p_ai <- ggplot() +
  geom_point(
    data = plotdat_nona  %>%
      filter(
        name == "ai"
      ),
    aes(
      x = rate * 100,
      y = value,
      size = forest_px_sum,
      col = forest_area_cut
    )
  ) +
  geom_ribbon(
    data = fit_ai.log.df,
    aes(
      x = (rate) * 100,
      ymin = exp(pred_lower),
      ymax = exp(pred_upper),
      group = forest_area_cut
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = fit_ai.log.df,
    aes(
      x = (rate) * 100,
      y = exp(pred),
      group = forest_area_cut
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = "Aggregation index (%)",
    #y = "%",
    col = "Forest cover"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#41ae76", "#005824")
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  ylim(50, 100) +
  guides(
    size = "none"
  ) +
  facet_wrap(
    ~forest_area_cut,
    scales = "free"
  )

p_ai_clean <- ggplot() +
  geom_point(
    data = plotdat_nona  %>%
      filter(
        name == "ai"
      ),
    aes(
      x = rate * 100,
      y = value,
      size = forest_px_sum,
      col = forest_area_cut
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Disturbance rate (%"~yr^-1*")"),
    y = "Aggregation index (%)",
    #y = "%",
    col = "Forest cover"
  ) +
  scale_color_manual(
    values = c("#99d8c9", "#41ae76", "#005824")
  ) +
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_size_continuous(
    range = c(0.01, 2)
  ) +
  ylim(50, 100) +
  guides(
    size = "none"
  )

p_np + theme(legend.position = "none") +
  p_area + theme(legend.position = "none") +
  p_ai +
  plot_layout(
    ncol = 1
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

p_area_clean + theme(legend.position = "none") +
  p_np_clean + theme(legend.position = "none") +
  p_ai_clean +
  plot_layout(
    ncol = 3
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

ggsave("results/figure03.pdf", width = 7.3 / 3 * 3, height = 7.5 / 3)

dat_core_np_area <- plotdat_nona  %>%
  filter(
    name == "area"
  ) %>%
  select(
    gridid,
    area = value,
    forest_px_sum,
    forest_area_cut,
    rate
  ) %>%
  left_join(
    plotdat_nona  %>%
      filter(
        name == "np"
      ) %>%
      select(
        gridid,
        np = value,
        forest_px_sum,
        forest_area_cut,
        rate
      )
  ) %>%
  left_join(
    plotdat_nona  %>%
      filter(
        name == "ai"
      ) %>%
      select(
        gridid,
        ai = value,
        forest_px_sum,
        forest_area_cut,
        rate
      )
  )

plotdat %>%
  mutate(
    forest_area_cut = cut(
      forest_px_sum,
      c(c(0, 0.2, 0.6) * 100^2 / 9e-04, max(forest_px_sum)),
      include.lowest = TRUE,
      labels = c("0-40%", "20-60%", ">60%")
    )
  ) %>%
  group_by(
    name,
    forest_area_cut,
    rate_cutoff = rate > 0.015
  ) %>%
  summarize(
    n = n(),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    mn = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n)
  ) %>%
  filter(
    name %in% c("area", "np")
  ) %>%
  View()

pl

# Maps ----

library(sf)

grid <- read_sf("data/gis/grid_100km_epsg3035.shp")
countries <- read_sf("data/gis/countries_europe_map.shp")
countries <- countries %>% st_transform(st_crs(grid))

centerpoints <- grid %>%
  st_centroid()

country_intersection <- st_intersection(centerpoints, countries)

write_csv(country_intersection, "data/gridcell_country.csv")

country_intersection <- country_intersection %>%
  st_drop_geometry() %>%
  group_by(
    gridid
  ) %>%
  summarise(
    country = unique(name0)[1]
  )

countries_map <- countries %>% 
  filter(name0 %in% country_intersection$country) %>%
  st_crop(grid)

grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "total"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = factor(ifelse(rate * 100 < 1.5, ifelse(rate * 100 > 1, 1, 0), 2)),
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Disturbance rate (%"~yr.^-1~")")
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  )

map1a <- grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "total"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = rate * 100,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c() +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Disturbance rate (%"~yr.^-1~")")
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

map1b_natural <- grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "natural"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = expansion_rate_n_mn * 100,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c(
    breaks = c(0.2, 0.4, 0.6, 0.8, 1) * 100,
    limits = c(0.2, 1) * 100
  ) +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Gap growth probability"),
    subtitle = "Natural disturbances"
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

map1b_harvest <- grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "harvest"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = expansion_rate_n_mn * 100,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c(
    breaks = c(0.2, 0.4, 0.6, 0.8, 1) * 100,
    limits = c(0.2, 1) * 100
  ) +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Gap growth probability"),
    subtitle = "Harvest"
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

map1c <- grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "total"
      ) %>%
      mutate(
        gap_growth_analytical = P_exp_approx(lambda = rate, r = 30, k = 8)
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = gap_growth_analytical * 100,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c(
    breaks = c(0.2, 0.4, 0.6, 0.8, 1) * 100,
    limits = c(0.2, 1) * 100
  ) +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Analytical solution (%"~yr.^-1~")")
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

map1d_natural <- grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "natural"
      ) %>%
      mutate(
        gap_growth_analytical = P_exp_approx(lambda = rate, r = 30, k = 8)
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = gap_growth_analytical - expansion_rate_n_mn,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_gradient2(
    limits = c(-0.35, 0.35)
  ) +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Difference analytical - empirical"),
    subtitle = "Natural disturbances"
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

map1d_harvest <- grid %>%
  left_join(
    expansion_rates_all %>%
      filter(
        subset == "harvest"
      ) %>%
      mutate(
        gap_growth_analytical = P_exp_approx(lambda = rate, r = 30, k = 8)
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = gap_growth_analytical - expansion_rate_n_mn,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_gradient2(
    limits = c(-0.35, 0.35)
  ) +
  scale_size_continuous(
    range = c(0.01, 0.8), 
    guide = "none"
  )  +
  theme_minimal() +
  labs(
    col = NULL,
    title = bquote("Difference analytical - empirical"),
    subtitle = "Harvest"
  ) +
  theme(
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 7),
    axis.text = element_blank(),
    legend.position = c(-0.05, 1.05),
    legend.justification = c(0, 1),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5
    )
  )

map1a + map1b_harvest + map1b_natural +
  map1c + map1d_harvest + map1d_natural +
  plot_layout(
    ncol = 3
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

map2 <- grid %>%
  left_join(
    plotdat_nona %>%
      filter(
        name == "area"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = (value),
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(0.05, 0.75), guide = "none")  +
  theme_minimal() +
  labs(
    col = bquote("Average patch size (ha)")
  ) +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    legend.position = "bottom",
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top"
    )
  )

map3 <- grid %>%
  left_join(
    plotdat_nona %>%
      filter(
        name == "np"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = (value) / 100^2,
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(0.05, 0.75), guide = "none")  +
  theme_minimal() +
  labs(
    col = bquote("Number of patches (100"~km^-2*")")
  ) +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    legend.position = "bottom",
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top"
    )
  )

map4 <- grid %>%
  left_join(
    plotdat_nona %>%
      filter(
        name == "ai"
      ),
    by = "gridid"
  ) %>%
  st_centroid() %>%
  ggplot(
    data = .
  ) +
  geom_sf(
    data = countries_map,
    fill = "grey"
  ) +
  geom_sf(
    aes(
      col = (value),
      size = forest_px_sum * 9e-04 / 100^2
    )
  ) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(0.05, 0.75), guide = "none")  +
  theme_minimal() +
  labs(
    col = bquote("Aggregation index (%)")
  ) +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    legend.position = "bottom",
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ) +
  guides(
    col = guide_colorbar(
      direction = "horizontal",
      title.position = "top"
    )
  )

map2 + map3 + map4 + plot_layout(ncol = 3) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )
