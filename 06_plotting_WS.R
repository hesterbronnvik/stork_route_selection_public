### Code to generate plots for Bronnvik et al., 2024
### Hester Bronnvik
### 2023-08-29
### hbronnvik@ab.mpg.de

library(gridExtra)
library(grid)
library(cowplot)
library(tidyverse)

# the processed data used to look at relationships within the data
a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/a_data_2023-09-05.rds")

# the model results and predictions used to visualize results
pre_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_pre_2023-09-05.rds")
post_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_post_2023-09-05.rds")

pre_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_pre_2023-09-05.rds")
post_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_post_2023-09-05.rds")

# a universal palette (cool to warm)
# cerulean, Munsell, verdigris, Tiffany, light orange, melon, salmon, bittersweet
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# universal facet labels
fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")

### Plots for the models

# Figure 1
# the coefficients of the models for spring and fall
pre_graph <- confint(pre_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!grepl("id", Factor)) %>% 
  mutate(Variable = c("Conspecific density", "Uplift", "Spring migrations", #"Step length", "Turning angle", 
                      "Conspecifics X Uplift", "Conspecifics X Migrations", "Uplift X Migrations", "Conspecifics X Uplift X Migrations"))
colnames(pre_graph)[2:3] <- c("Lower", "Upper") 

pre_graph$p <- as.data.frame(summary(pre_mod)[[6]][1])[,4]
pre_graph$significance <- ifelse(pre_graph$p < 0.001, "***", 
                                 ifelse(pre_graph$p > 0.001 & pre_graph$p < 0.01, "**",
                                        ifelse(pre_graph$p > 0.01 & pre_graph$p < 0.05, "*",
                                               "")))
pre_graph <- pre_graph %>% 
  mutate(Variable = factor(Variable, levels = c("Conspecifics X Uplift X Migrations", "Uplift X Migrations",
                                                "Conspecifics X Uplift", "Conspecifics X Migrations", "Spring migrations", 
                                                "Uplift", "Conspecific density"))) %>% 
  arrange(Variable)

post_graph <- confint(post_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!grepl("id", Factor)) %>% 
  mutate(Variable = c("Conspecific density", "Uplift", "Fall migrations", #"Step length", "Turning angle", 
                      "Conspecifics X Uplift", "Conspecifics X Migrations", "Uplift X Migrations", "Conspecifics X Uplift X Migrations"))
colnames(post_graph)[2:3] <- c("Lower", "Upper") 

post_graph$p <- as.data.frame(summary(post_mod)[[6]][1])[,4]
post_graph$significance <- ifelse(post_graph$p < 0.001, "***", 
                                  ifelse(post_graph$p > 0.001 & post_graph$p < 0.01, "**",
                                         ifelse(post_graph$p > 0.01 & post_graph$p < 0.05, "*",
                                                "")))
post_graph <- post_graph %>% 
  mutate(Variable = factor(Variable, levels = c("Conspecifics X Uplift X Migrations", "Uplift X Migrations",
                                                "Conspecifics X Uplift", "Conspecifics X Migrations", "Fall migrations", 
                                                "Uplift", "Conspecific density"))) %>% 
  arrange(Variable)

post_coefs <- ggplot(post_graph %>% filter(Variable != "Fall migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Fall") +
  theme_classic() +
  ggplot2::annotate(geom="text", x=16, y=c(1:6), label=post_graph$significance[post_graph$Variable != "Fall migrations"], color="black", size = 6)  +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.5, 'cm')) 
pre_coefs <- ggplot(pre_graph %>% filter(Variable != "Spring migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Spring") +
  theme_classic() +
  ggplot2::annotate(geom="text", x=5, y=c(1:6), label=pre_graph$significance[pre_graph$Variable != "Spring migrations"], color="black", size = 6)  +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.5, 'cm'))

coefs_plot <-  plot_grid(post_coefs, pre_coefs, labels = c("A", "B"), nrow = 2,
                         align = 'v', axis = 'l')
coefs_plot

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/coefs_glmm_05-09.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
coefs_plot
dev.off()

# Figure 2
# the model predictions for the 3-term interaction
my_labels <- c("Low", "Mean", "High")
inter_plots <- function(data, x_id, y_id, xlab, ylab, xmin, xmax){
  ggplot(data, aes_string(x_id, y_id, fill = "probs")) +
    geom_tile(color = "white", lwd = 0, linetype = 1) +
    scale_x_continuous(expand = c(0, 0), breaks = c(-0.75, 0, 1, 2), labels = scales::number_format(accuracy = 1)) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = c(min(data$sqrt_ud), mean(data$sqrt_ud), max(data$sqrt_ud)), 
                       labels = my_labels) +
    labs(x = "", y = "", fill = "Selection \nprobability") +
    scale_fill_gradientn("Selection \nprobability", colors = colfunc(135)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 15),
          axis.line = element_line(linewidth = 1.2),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "none",
          # aspect.ratio = 1,
          axis.title.x=element_text(color = "white"),
          axis.text.x=element_text(color = "white"),
          axis.ticks.x=element_blank(),
          plot.margin=unit(c(-1,1,-1,0), "cm"))
}

spring_plots <- lapply(1:8, function(x){
  p <- inter_plots(pre_preds %>% filter(interaction == "uplift_migration_ud" & migrations == x),
                   "w_star", "sqrt_ud",
                   "Uplift (m/s)", "Conspecific density", 
                   1, 9) + labs(title = " ") + theme(axis.ticks.y=element_blank(),
                                                     axis.text.y=element_blank(),
                                                     plot.margin=unit(c(-1,4,-1,-1.5), "cm"))
  yloc <- pre_preds %>% 
    filter(interaction == "uplift_migration_ud") %>% 
    dplyr::select(sqrt_ud) %>% 
    deframe() %>% 
    mean()
  xloc <- pre_preds %>% 
    filter(interaction == "uplift_migration_ud") %>% 
    dplyr::select(w_star) %>% 
    deframe() %>% 
    max() + .75
  
  if(x == 1){p <- p + 
    labs(title = "Spring")}
  if(x == 8){p <- p + theme(axis.title.x=element_text(color = "black"),
                            axis.text.x=element_text(color = "black"),
                            axis.ticks.x=element_blank())}
  p
})

fall_plots <- lapply(1:9, function(x){
  p <- inter_plots(post_preds %>% filter(interaction == "uplift_migration_ud" & migrations == x),
                   "w_star", "sqrt_ud",
                   "Uplift (m/s)", "Conspecific density", 
                   1, 9) + labs(title = " ")
  yloc <- post_preds %>% 
    filter(interaction == "uplift_migration_ud") %>% 
    dplyr::select(sqrt_ud) %>% 
    deframe() %>% 
    mean()
  xloc <- post_preds %>% 
    filter(interaction == "uplift_migration_ud") %>% 
    dplyr::select(w_star) %>% 
    deframe() %>% 
    max() + .75
  
  if(x == 1){p <- p + 
    labs(title = "Fall")}
  if(x == 9){p <- p + theme(axis.title.x=element_text(color = "black"),
                            axis.text.x=element_text(color = "black"),
                            axis.ticks.x=element_blank())}
  p
})

yleft <- textGrob("Conspecific density", rot = 90, gp = gpar(fontsize = 25))
xbottom <- textGrob("Uplift (m/s)", gp = gpar(fontsize = 25))

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/pcol2.png",
    height = 23.4, width = 16.6/2, units = "in", res = 500)
gridExtra::grid.arrange(fall_plots[[1]], spring_plots[[1]], fall_plots[[2]], spring_plots[[2]], fall_plots[[3]], spring_plots[[3]], 
                        fall_plots[[4]], spring_plots[[4]], fall_plots[[5]], spring_plots[[5]], fall_plots[[6]], spring_plots[[6]], fall_plots[[7]], spring_plots[[7]], 
                        fall_plots[[8]], spring_plots[[8]], fall_plots[[9]], 
                        ncol = 2, nrow = 9, left = yleft, bottom = xbottom, 
                        vp = viewport(width=0.9, height=0.9))
dev.off()

# Supplementary Figure S1
# the model predictions covering the interaction of age and uplift in the spring
inter_plots <- function(data, x_id, y_id, xlab, ylab, xmin, xmax){
  ggplot(data, aes_string(x_id, y_id, fill = "probs")) +
    geom_tile(color = "white", lwd = 0, linetype = 1) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-0.75, 0, 1, 2), labels = scales::number_format(accuracy = 1)) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1:8)) +
    labs(x = xlab, y = ylab, fill = "Selection \nprobability") +
    scale_fill_gradientn("Selection \nprobability", colors = colfunc(135)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 15),
          axis.line = element_line(linewidth = 1.2),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          aspect.ratio = 1,
          legend.key.size = unit(0.75, 'cm'))
}
w_age_s <- inter_plots(pre_preds %>% filter(interaction == "uplift_migration"),
                       "migrations", "w_star",
                       "Number of spring migrations", "Uplift (m/s)", 
                       1, 9)
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/uplift_age.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
w_age_s
dev.off()

### Plots for data exploration
a_data$datestamp <- a_data$timestamp
year(a_data$datestamp) <- 2024

# Figure 3
# look at the distributions of variances
k_data <- a_data %>% 
  group_by(stratum) %>% 
  mutate(obs_diff = w_star[which(used == 1)]-mean(w_star[which(used == 0)]),
         strat_var_w = var(w_star),
         strat_var_s = var(ud_pdf)) %>% 
  slice(1) %>% 
  ungroup()
w_dens <- ggplot(k_data, aes(log(strat_var_w), fill = as.factor(migrations))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = colfunc(9)) +
  labs(x = "log per-hour variance in uplift", y = "Density", fill = "Migrations") +
  theme_classic()  +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm')) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
s_dens <- ggplot(k_data, aes(log(strat_var_s), fill = as.factor(migrations))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = colfunc(9)) +
  labs(x = "log per-hour variance in conspecific density", y = "Density", fill = "Migrations") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(0.75, 'cm')) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
legend <- get_legend(
  # create some space to the left of the legend
  s_dens + theme(legend.box.margin = margin(0, 0, 0, 12))
)
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/stratum_densities.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggpubr::ggarrange(s_dens + theme(legend.position="none"), w_dens + theme(legend.position="none"), 
                  ncol = 1, common.legend = T, legend.grob = legend, legend = "right")
dev.off()

# Supplementary Figure S2
# look at the number of recordings of migration on each day for each age group
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/migration_timing_sp23.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
a_data %>% 
  mutate(yd = date(alignment)) %>% 
  group_by(yd, migrations) %>% 
  mutate(count = n()) %>% 
  slice(1) %>%
  ungroup() %>% 
  arrange(migrations) %>% 
  ggplot(aes(x = yd, y = count, group = as.factor(migrations), color = as.factor(migrations))) +
  geom_segment(aes(alpha = 0.7, x=yd, xend=yd, y=0, yend=count), linewidth = 1.5) +
  geom_point(size = 3) +
  labs(x = "Day", y = "Observations", color = "Migration") +
  scale_color_manual(values = c("#0081A7", "#0098B0", "#00AFB9", "#7FC4B8", "#F7A58F", "#F27E71", "#F07167", "#ED5145", "#EE5E53")) +
  guides(alpha = "none") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%m-%d") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75, 'cm'))
dev.off()

# Supplementary Figure S3
# look at uplift distributions per age group
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/seasonal_uplifts.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(a_data %>% mutate(season = ifelse(season == "pre", "Spring", "Fall")), aes(as.factor(migrations), w_star, fill = season)) +
  geom_boxplot(lwd = 1.5, aes(alpha = forcats::fct_rev(as.factor(migrations)))) +
  scale_fill_manual(values = c("#FE6D5D", "#0096CC"))  +
  guides(alpha = "none") +
  labs(x = "Migrations", y = "Uplift (m/s)", fill = "Season") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75, 'cm'))
dev.off()