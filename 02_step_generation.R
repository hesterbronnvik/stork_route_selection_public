### Step generation for step selection functions using only the complete routes
### Adapted from E. Nourani 2021 https://github.com/mahle68/global_seascape_public/blob/main/step_generation.R
### Hester Bronnvik
### hbronnvik@ab.mpg.de
### 2023-02-08

library(move)
library(CircStats)
library(circular)
library(fitdistrplus)
library(lubridate)
library(tidyverse)
#import required functions

NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}

# set criteria for the steps
hr <- 60 #minutes; determine the sub-sampling interval
tolerance <- 15 #minutes; tolerance for sub-sampling
n_alt <- 100 #number of alternative steps. (create more than are strictly needed to remove some later)
meters_proj <- CRS("+proj=moll +ellps=WGS84")#Mollweide projection (in meters) for accurate calculation of length
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# the location data during migration
full_ml <- readRDS("/home/hbronnvik/Documents/storkSSFs/migration_locations_40km70km30daySpeed_2023-09-27.rds")

# only migratory flight locations
# filter ground speeds and daily distances (speed is scalar and a few missing data can affect it dramatically)
full_ml <- full_ml %>% 
  filter(daily_dist > 40*1000 & ground_speed_15 > 2)
# arrange to the satisfaction of the move object
full_ml <- full_ml %>% 
  arrange(individual.id, timestamp) %>% 
  as.data.frame()
# make a move object
mv <- move(x = full_ml$location.long, y = full_ml$location.lat, time = full_ml$timestamp, data =full_ml, 
           proj = CRS("+proj=longlat +datum=WGS84 +no_defs"), animal = full_ml$trackID)

## Burst the tracks
(b <- Sys.time()) 
sp_obj_ls <- lapply(split(mv), function(track){ #for each track,
  print(track@idData$trackID)

  #--STEP 1: thin the data 
  track_th <- track %>%
    thinTrackTime(interval = as.difftime(hr, units = 'mins'),
                  tolerance = as.difftime(tolerance, units = 'mins')) #the un-selected bursts are the large gaps between the selected ones
  
  #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst. Longer gaps will divide the bursts) 
  track_th$selected <- c(as.character(track_th@burstId),NA) #assign "selected" as a variable
  track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define a value for first row
  
  if(nrow(track_th@data) == 1){
    track_th@data$burst_id <- track_th$burst_id
  } else {for(i in 2:nrow(track_th@data)){
    #if(i == nrow(track_th@data)){
    #  track_th@data$burst_id[i] <- NA #why?!
    #} else
    if(track_th@data[i-1,"selected"] == "selected"){
      track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"]
    } else {
      track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"] + 1
    }
  }
  }
  
  #convert back to a move object (from move burst)
  track_th <- as(track_th,"Move")
  
  #--STEP 3: calculate step lengths and turning angles 
  #sl_ and ta_ calculations should be done for each burst.
  burst_ls <- split(track_th, track_th$burst_id)
  burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with fewer than 3 observations
  
  burst_ls <- lapply(burst_ls, function(burst){
    burst$step_length <- c(distance(burst),NA) 
    burst$turning_angle <- c(NA,turnAngleGc(burst),NA)
    burst
  })
  
  #put burst_ls into one dataframe
  bursted_sp <- do.call(rbind, burst_ls)
  
  #reassign values
  if(length(bursted_sp) >= 1){
    bursted_sp$track <- track@idData$trackID
    bursted_sp$individual.id <- track@idData$individual.id
  }
  #bursted_sp$track<-track@idData$seg_id 
  
  return(bursted_sp)
  
}) %>% 
  Filter(function(x) length(x) > 1, .) #remove segments with no observation

if(length(sp_obj_ls) > 0){
  #--STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  burst_df <- sp_obj_ls %>%  
    reduce(rbind) %>% 
    as.data.frame() %>% 
    dplyr::select(-c("coords.x1","coords.x2"))
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
  mu <- mean.circular(rad(burst_df$turning_angle[complete.cases(burst_df$turning_angle)]))
  kappa <- est.kappa(rad(burst_df$turning_angle[complete.cases(burst_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!!
  sl <- burst_df$step_length[complete.cases(burst_df$step_length) & burst_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  # plot turning angle and step length distributions
  png(paste0("/home/hbronnvik/Documents/storkSSFs/figures/", hr, "_", tolerance, ".png"),
      width = 11.5, height = 8, units = "in", res = 500)
  par(mfrow=c(1,2))
  hist(sl, freq=F, main="", xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  hist(rad(burst_df$turning_angle[complete.cases(burst_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  dev.off()

  #diagnostic plots for step length distribution
  png(paste0("/home/hbronnvik/Documents/storkSSFs/figures/", hr, "_", tolerance, "_diag.png"),
      width = 11.5, height = 8, units = "in", res = 500)
  plot(fit.gamma1)
  dev.off()
  
  #--STEP 5: produce alternative steps
  used_av_track <- lapply(sp_obj_ls, function(track){ #for each trackment
    
    lapply(split(track,track$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
      if(kappa != Inf){
      
      
        lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
          
          current_point<- burst[this_point,]
          previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
          used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
          
          #calculate bearing of previous point
          #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
          prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                            previous_point@coords[,1], current_point@coords[,1])
          
          current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
          
          #randomly generate n alternative points
          rnd <- data.frame(turning_angle = as.vector(rvonmises(n = n_alt, mu = mu, kappa = kappa)), #randomly generate n step lengths and turning angles
                            step_length = rgamma(n = n_alt, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000) %>% 
            #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
            mutate(lon = current_point_m@coords[,1] + step_length*cos(turning_angle),
                   lat = current_point_m@coords[,2] + step_length*sin(turning_angle))
          
          
          #convert back to lat-lon proj
          rnd_sp <- rnd
          coordinates(rnd_sp) <- ~lon+lat
          proj4string(rnd_sp) <- meters_proj
          rnd_sp <- spTransform(rnd_sp, wgs)
          
          #put used and available points together
          df <- used_point@data %>%  
            mutate(x = location.long,
                   y = location.lat) %>% 
            slice(rep(row_number(), n_alt+1)) %>% #paste each row n_alt times for the used and alternative steps
            mutate(location.long = c(head(x,1),rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
                   location.lat = c(head(y,1),rnd_sp@coords[,2]),
                   turning_angle = c(head(turning_angle,1),deg(rnd_sp$turning_angle)),
                   step_length = c(head(step_length,1),rnd_sp$step_length),
                   used = c(1,rep(0,n_alt)))  %>%
            dplyr::select(-c("x","y")) %>% 
            rowwise() %>% 
            mutate(heading = NCEP.loxodrome.na(lat1 = current_point@coords[,2], lat2 = location.lat, lon1 = current_point@coords[,1], lon2 = location.long)) %>% 
            as.data.frame()
          
          df
        
        }) %>% 
        reduce(rbind)
      }else{
          df <- burst@data %>% 
            mutate(used = NA,
                   heading = NA)
        }
      
    }) %>% 
      reduce(rbind)
    
  }) %>% 
    reduce(rbind)
}
  

Sys.time() - b # Time difference of 1.044653 hours

#create one dataframe with movebank specifications
used_av_df <- used_av_track %>% 
  dplyr::select( c("timestamp", "location.lat", "location.long", "selected", "individual.id",  "burst_id", "step_length", "turning_angle", "track", "step_id", "used", "heading")) %>% 
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>%
  mutate(stratum = paste(track, burst_id, step_id, sep = "_")) %>% #create unique stratum id
  as.data.frame()

# take only the steps created over land (this is why we made more than we need)

# call in a buffered raster of land outlines + 30km
mm <- raster::raster("/home/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
spatMM <- terra::rast(mm)
poi_vals <- terra::extract(spatMM, terra::vect(used_av_df, geom = c("location.long", "location.lat")))[,2]
unique(poi_vals)
used_av_df$poi_vals <- poi_vals
used_av_df$land <- ifelse(is.na(poi_vals), "water", "land")

# there are no actual, observed points classified as being over water
table(is.na(used_av_df$poi_vals), used_av_df$used)

used_av_df %>% 
  filter(is.na(poi_vals)) %>% 
  group_by(stratum) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  dplyr::select(n, stratum, used, step_length) %>% 
  arrange(desc(n))

# remove the "available" locations that we generated over the sea and pare down the data
# so that all strata are the same size (not strictly needed for Bayesian approaches)
used_ad <- used_av_df %>% 
  filter(used == 1)
avail_ad <- used_av_df %>% 
  filter(used == 0 & land == "land") %>% 
  group_by(stratum) %>% 
  slice_sample(n = 49) %>% 
  ungroup()

used_av_df <- used_ad %>% 
  rbind(avail_ad) %>% 
  arrange(individual.id, stratum)

used_av_df %>% 
  group_by(stratum) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  dplyr::select(n) %>% 
  unique()

# saveRDS(used_av_df, file = "/home/hbronnvik/Documents/storkSSFs/burst_data/used_av_df_230612.rds")

#rename lat and lon columns
colnames(used_av_df)[c(2,3)] <- c("location-lat", "location-long")

# add a basis for splitting the df into small enough chunks to run on Movebank
used_av_df <- used_av_df %>% 
  mutate(splitter = c(rep("A", times = 1000000), rep("B", times = 1000000), rep("C", times = n()-2000000)))

# save to files for Movebank annotation
lapply(split(used_av_df, used_av_df$splitter), function(x){
  df <- x %>% 
    dplyr::select(-splitter)
  write.csv(df, file = paste0("/home/hbronnvik/Documents/storkSSFs/burst_data/", unique(x$splitter), "_", Sys.Date(), ".csv"))
})
