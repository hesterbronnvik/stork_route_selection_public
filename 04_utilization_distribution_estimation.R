### Code for estimating the utilization distributions of white storks migrating on the western flyway
### Hester Bronnvik
### 2023-07-13
### hbronnvik@ab.mpg.de

library(move)
library(ctmm)
library(raster)
library(terra)
library(data.table)
library(lubridate)
library(tidyverse)

# read in the saved data for the 397 birds for all locations during migration
locs <- readRDS("/home/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_2023-08-30.rds")

# the total number of migrations attempted and completed by each animal in each season 
# this file was generated in wip_seg.R before any of the migrations were filtered based on length or success
meta <- readRDS("/home/hbronnvik/Documents/storkSSFs/metadata_migration_tracks_2023-08-30.rds")

# align these locations in time
locs$datestamp <- locs$timestamp
year(locs$datestamp) <- 2024
second(locs$datestamp) <- 00
locs$datestamp <- round_date(locs$datestamp, "hour")

# read in the saved data for the 158 birds that completed migrations
# these have already been burst (see 02_step_generation.R) and annotated with uplift (see 03_uplift_annotation.R)
a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/full_annotated_310823.rds")

# align these locations in time
a_data$alignment <- round_date(a_data$timestamp, unit = "hour")
year(a_data$alignment) <- 2024
second(a_data$alignment) <- 00

# get the outlines of land to remove the possibility of information over water
# tmax <- raster::getData('worldclim', var = "tmax", res = 10)
# mask <- crop(raster(tmax, 1), raster::extent(min(first_locs$location.long), max(first_locs$location.long), min(first_locs$location.lat), max(first_locs$location.lat)))
# raster::plot(mask)
# # the Straits of Gibraltar are at most 60km across
# mm <- buffer(mask, 30000)
# raster::plot(mm)
# writeRaster(mm, "C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")

# call in a buffered raster of land outlines + 30km
mm <- raster::raster("/home/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
outlines <- raster::rasterToPolygons(mm, dissolve=TRUE)

# now, take one of the 158 animals at a time
# take it out of the 367 birds that migrated, and use the rest to build one AKDE per hour
# run the one bird across those AKDEs, annotating its locations with the AKDE for that hour
# lather, rinse, repeat 
start_time <- Sys.time()
HR <- lapply(unique(a_data$individual.id), function(x){
  print(x)
  # isolate the locations of a single, focal bird
  ind <- a_data %>% 
    filter(individual.id == x)
  # take all the locations that are not from that bird
  soc <- locs %>% 
    filter(individual.id != x) %>% 
    mutate(individual.id = as.character(datestamp),
           individual.local.identifier = as.character(datestamp)) %>% 
    group_by(individual.id) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    # only use hours with more than 25 locations
    filter(count >= 25)%>% 
    filter(datestamp %in% unique(a_data$alignment))
  # compress into a leap year
  year(soc$timestamp) <- 2024
  
  # take all of the migratory locations within the hours of migrations by the focal bird
  soc <- soc %>% 
    filter(datestamp %in% unique(ind$alignment))
  # make a list of telemetry objects
  soc_tele <- suppressWarnings(suppressMessages(as.telemetry(soc, projection = "ESRI:54009")))
  # make a list of model fits
  fits <- lapply(soc_tele, function(x){
    tryCatch({
      guess <- ctmm.guess(x, interactive=FALSE)
      # fit the models
      fit <- ctmm.select(x, guess)
      return(fit)
    }, error = function(e){print(geterrmessage(), quote = F)})
  }) %>% suppressWarnings() # warning: Duplicate timestamps require an error model.
  # (we ignore this warning because we are deliberately using only one timestamp and a population not an individual)
  # estimate utilization distributions of all (the non-focal) storks in each hour 
  # cropped to the buffered map, and on a grid of set extent and 30km resolution (approximately matching the weather data)
  kdes <- akde(soc_tele, fits, SP = outlines, grid = list(dr = 30000, extent = extent(-2468787, 3006683, 0, 6876759)))
  # go through the KDEs to annotate the tracks
  akde_data <- lapply(kdes, function(y){
    # the hour covered by this KDE
    kde_hour <- y@info$identity
    # the locations in that hour
    tracks <- ind %>% 
      filter(as.character(alignment) == kde_hour)
    if(nrow(tracks) > 0){
      # convert to a spatial object
      coordinates(tracks) <- ~ long + lat
      proj4string(tracks) <- CRS("EPSG:4326")
      tracks <- spTransform(tracks, "ESRI:54009")
      # the KDE as a raster ("PDF" gives the average probability density per cell)
      ud <- raster(y, DF = "PDF")
      # extract the UD value
      vals <- raster::extract(ud, tracks)
      # append
      df <- ind %>% 
        filter(as.character(alignment) == kde_hour) %>% 
        mutate(UD_PDF = vals)
      return(df)
    }
  }) %>% 
    discard(is.null)
  akde_data <- rbindlist(akde_data)
  print(paste0("Annotated data for individual ", which(unique(a_data$individual.id) == x), " of ", 
               length(unique(a_data$individual.id)), "."))
  return(akde_data)
}) %>% rbindlist()
Sys.time() - start_time # Time difference of 18.4626 hours

a_data <- a_data %>% 
  full_join(HR)
# saveRDS(HR, file = "/home/hbronnvik/Documents/storkSSFs/annotations/HR_only_030923.rds")
# saveRDS(a_data, file = "/home/hbronnvik/Documents/storkSSFs/annotations/HR_030923.rds")