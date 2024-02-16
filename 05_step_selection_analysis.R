### Code for running step selection analyses
### Hester Bronnvik
### 2023-07-13
### hbronnvik@ab.mpg.de

library(corrr)
library(performance)
library(glmmTMB)
library(tidyverse)

# the annotated tracks
a_data <- readRDS(file = "/home/hbronnvik/Documents/storkSSFs/annotations/HR_030923.rds")

# the total number of migrations attempted and completed by each animal in each season 
# this file was generated in wip_seg.R before any of the migrations were filtered based on length or success
meta <- readRDS("/home/hbronnvik/Documents/storkSSFs/metadata_migration_tracks_2023-08-30.rds")

# only the metadata we need
meta <- meta %>% 
  rename(track = trackID) %>% 
  dplyr::select(individual.id, track, journey_number, total_journeys)

# add the number of journeys
a_data <- a_data %>% 
  arrange(timestamp) %>% 
  left_join(meta)

# scale the predictors
a_data <- a_data %>% 
  rename(ud_pdf = UD_PDF,
         migrations = journey_number) %>% 
  mutate_at(c("migrations", "w_star","ud_pdf", "step_length", "turning_angle", "blh"),
            list(z = ~(scale(.)))) %>% 
  # remove the geopotential heights
  dplyr::select(-colnames(.)[15:28])

# square root transform the social influence proxy for normality
a_data <- a_data %>% 
  mutate(sqrt_ud = sqrt(ud_pdf),
         sqrt_ud_z = scale(sqrt_ud))

# save out the data 
# saveRDS(a_data, file = paste0("/home/hbronnvik/Documents/storkSSFs/a_data_", Sys.Date(), ".rds"))

a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/a_data_2023-09-05.rds") %>%  # or most recent version
  mutate(season = ifelse(grepl("fall", track), "post", "pre"))

# the number of journeys that passed all the filtering
a_data %>% 
  filter(season == "pre") %>% 
  group_by(track) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(migrations) %>% 
  table()
a_data %>% 
  filter(season == "post") %>% 
  group_by(track) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(migrations) %>% 
  table()

# look at correlation
a_data %>% 
  dplyr::select(c(ud_pdf, w_star, step_length, turning_angle, migrations)) %>% 
  correlate() # w* and conspecific density = 0.0169

# run Kolmogorov-Smirnov tests to see how the environmental variables are distributed in each migration
k_data <- a_data %>% 
  group_by(stratum) %>% 
  mutate(strat_var_w = var(w_star),
         strat_var_s = var(ud_pdf)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(group = paste(season, migrations, sep = "_"))

k_data_ls <- split(k_data, k_data$group)
ks_stats <- lapply(k_data_ls, function(k){
  variable <- lapply(c("strat_var_w", "strat_var_s"), function(v){
    # generalze the name for the group of interest
    vk <- k %>% 
      rename("variable" = all_of(v))
    # take all the other data in the right season but other migrations
    au <- k_data %>% 
      rename("variable" = all_of(v)) %>% 
      filter(group != unique(vk$group) & season == unique(vk$season))
    output <- data.frame(season = unique(vk$season),
                         migration = unique(vk$migrations),
                         variable = v,
                         t_stat = round(ks.test(vk$variable, au$variable)[[1]], 3),
                         p_val = format(round(ks.test(vk$variable, au$variable)[[2]], 4), scientific = F))
    return(output)
  }) %>% reduce(rbind)
  return(variable)
}) %>% reduce(rbind) %>% 
  arrange(variable)

# STEP 1: run the model ------------------------------------------------------------------ 
#this is based on Muff et al:
#https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1

a_d1 <- a_data %>% 
  filter(season == "post") %>% 
  mutate(stratum_ID = as.factor(stratum),
         startum_ID = as.numeric(stratum_ID),
         individual.id = as.numeric(individual.id))

TMB_struc <- glmmTMB(used ~ -1 + sqrt_ud_z*w_star_z*migrations_z + step_length_z + turning_angle_z + 
                       (1|stratum_ID) + 
                       (0 + sqrt_ud_z | individual.id) + 
                       (0 + w_star_z | individual.id) + 
                       (0 + migrations_z | individual.id), 
                     family = poisson,
                     data = a_d1, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA,1:3))), # 3 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3),0,0,0))) #add a 0 for each random slope. in this case, 2

start_time <- Sys.time()
TMB_M <- glmmTMB:::fitTMB(TMB_struc)
Sys.time() - start_time # Time difference of 2.992602 mins
summary(TMB_M) # fall AIC: 428810.6

TMB_struc <- glmmTMB(used ~ -1 + sqrt_ud_z*w_star_z*migrations_z + step_length_z + turning_angle_z + 
                       (1|stratum_ID) + 
                       (0 + sqrt_ud_z | individual.id) + 
                       (0 + w_star_z | individual.id) + 
                       (0 + migrations_z | individual.id), 
                     family = poisson, ziformula = ~1,
                     data = a_d1, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA,1:3))), # 3 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3),0,0,0))) #add a 0 for each random slope. in this case, 2

start_time <- Sys.time()
TMB_M_zi <- glmmTMB:::fitTMB(TMB_struc)
Sys.time() - start_time 
summary(TMB_M_zi) # fall AIC: 428319.6

# the AIC is slightly better without considering zero-inflation, so we proceed with that model

lapply(split(a_data, a_data$season), function(data){
  season <- unique(data$season)
  # prep the data
  prep_d <- data %>% 
    mutate(stratum_ID = as.factor(stratum),
           stratum_ID = as.numeric(stratum_ID),
           individual.id = as.numeric(individual.id))
  # define the structure of the model
  TMB_struc <- glmmTMB(used ~ -1 + sqrt_ud_z*w_star_z*migrations_z + #step_length_z + turning_angle_z + 
                         (1|stratum_ID) + 
                         (0 + sqrt_ud_z | individual.id) + 
                         (0 + w_star_z | individual.id) + 
                         (0 + migrations_z | individual.id), 
                       family = poisson,
                       data = prep_d, doFit = FALSE,
                       #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                       map = list(theta = factor(c(NA,1:3))), # 3 is the n of random slopes
                       #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                       start = list(theta = c(log(1e3),0,0,0))) #add a 0 for each random slope. in this case, 2
  # fit the model
  TMB_M <- glmmTMB:::fitTMB(TMB_struc)
  # save the model
  saveRDS(TMB_M, file = paste0("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_", season, "_", Sys.Date(), ".rds"))
  
  # the second part demands that we create values and have the model predict use of them
  
  ## make a grid of the length of the data to fill in with predicted values (one for each predictor)
  #to make sure the predictions cover the parameter space, create a dataset with all possible combinations. one per interaction term. merge later on
  grd_up <- expand.grid(x = (1:max(prep_d$migrations)),
                        y = seq(from = quantile(prep_d$w_star, 0.025, na.rm = T), to = quantile(prep_d$w_star, 0.975, na.rm = T), length.out = 15)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
    rename(migrations = x,
           w_star = y) %>% 
    mutate(sqrt_ud = mean(prep_d$sqrt_ud, na.rm = T), #set other variables to their mean
           # migrations = mean(prep_d$journey_number),
           step_length = mean(prep_d$step_length, na.rm = T),
           turning_angle = mean(prep_d$turning_angle, na.rm = T),
           interaction = "uplift_migration")
  
  grd_soc <- expand.grid(x = (1:max(prep_d$migrations)),
                         y = seq(from = min(prep_d$sqrt_ud, na.rm = T), to = max(prep_d$sqrt_ud, na.rm = T), length.out = 15)) %>% # n = 135
    rename(migrations = x,
           sqrt_ud = y) %>% 
    mutate(w_star = mean(prep_d$w_star), #set other variables to their mean
           # migrations = mean(prep_d$journey_number),
           step_length = mean(prep_d$step_length, na.rm = T),
           turning_angle = mean(prep_d$turning_angle, na.rm = T),
           interaction = "ud_migration")
  
  grd_soc_up <- expand.grid(x = seq(from = min(prep_d$w_star, na.rm = T), to = max(prep_d$w_star, na.rm = T), length.out = max(prep_d$migrations)),
                            y = seq(from = min(prep_d$sqrt_ud, na.rm = T), to = max(prep_d$sqrt_ud, na.rm = T), length.out = 15)) %>% # n = 135
    rename(w_star = x,
           sqrt_ud = y) %>% 
    mutate(#set other variables to their mean
      migrations = mean(prep_d$migrations),
      step_length = mean(prep_d$step_length, na.rm = T),
      turning_angle = mean(prep_d$turning_angle, na.rm = T),
      interaction = "ud_up")
  
  grd_3 <- expand.grid(x = (1:max(prep_d$migrations)),
                       y = seq(from = quantile(prep_d$w_star, 0.025, na.rm = T), to = quantile(prep_d$w_star, 0.975, na.rm = T), length.out = 25),
                       z = seq(from = quantile(prep_d$sqrt_ud, 0.025, na.rm = T), to = quantile(prep_d$sqrt_ud, 0.975, na.rm = T), length.out = 25)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
    rename(migrations = x,
           w_star = y,
           sqrt_ud = z) %>% 
    mutate(step_length = mean(prep_d$step_length, na.rm = T),
           turning_angle = mean(prep_d$turning_angle, na.rm = T),
           interaction = "uplift_migration_ud")
  
  grd_all <- bind_rows(grd_up, grd_soc, grd_soc_up, grd_3) 
  
  set.seed(770)
  n <- nrow(grd_all)
  
  new_data_only <- prep_d %>%
    group_by(stratum_ID) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
    # only keep the columns that I need
    dplyr::select(c("stratum_ID", "individual.id")) %>% 
    bind_cols(grd_all) %>% 
    #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
    mutate(w_star_z = (w_star - mean(a_data$w_star))/(sd(a_data$w_star)),
           sqrt_ud_z = (sqrt_ud - mean(a_data$sqrt_ud))/(sd(a_data$sqrt_ud)),
           migrations_z = (migrations - mean(a_data$migrations))/(sd(a_data$migrations)),
           step_length_z = (step_length - mean(a_data$step_length, na.rm = T))/(sd(a_data$step_length, na.rm = T)), 
           turning_angle_z = (turning_angle - mean(a_data$turning_angle, na.rm = T))/(sd(a_data$turning_angle, na.rm = T)))
  
  new_data <- prep_d %>% 
    mutate(interaction = "OG_data") %>% 
    dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
    bind_rows(new_data_only) %>% 
    #accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
    #(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."
    mutate(stratum_ID = NA)#,
  # individual.id = NA)
  
  # now that we have the values to predict, run the model on them
  preds <- predict(TMB_M, newdata = new_data, type = "link")
  
  preds_pr <- new_data %>% 
    mutate(preds = preds) %>% 
    rowwise() %>% 
    mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob
  
  inter_preds <- preds_pr %>% 
    filter(interaction != "OG_data") 
  
  saveRDS(inter_preds, file = paste0("/home/hbronnvik/Documents/storkSSFs/glmm_preds_", season, "_", Sys.Date(), ".rds"))
})

pre_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_pre_2023-09-05.rds")
post_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_post_2023-09-05.rds")

pre_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_pre_2023-09-05.rds")
post_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_post_2023-09-05.rds")

# model validation 
diagnose(pre_mod)
# "large negative and/or positive components in binomial or Poisson
# conditional parameters suggest (quasi-)complete separation"
# "Complete or quasi-complete separation occurs when there is a combination of 
# regressors in the model whose value can perfectly predict one or both outcomes." 
# (https://www.zeileis.org/news/biasreduction/)

# extract coefficient estimates and confidence intervals
confint(pre_mod)

# extract individual-specific random effects:
# ranef(pre_mod)[[1]]$individual.id

#calculate the RMSE
performance_rmse(pre_mod) # 0.139612

# the same for fall
confint(post_mod)
performance_rmse(post_mod) # 0.1381841
