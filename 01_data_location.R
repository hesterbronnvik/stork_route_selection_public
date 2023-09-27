### Code to determine which birds to include, clean their data, and define migration.
### Hester Brønnvik
### hbronnvik@ab.mpg.de
### 29.09.2022

library(lubridate)
library(geosphere)
library(move)
library(stringr)
library(tidyverse)
library(data.table)
library(mapview)

# required information
d_thresh <- 40000 # meters
w_thresh <- 6 # weeks
s_thresh <- 30 # days
l_thresh <- 3 # degrees latitude
g_thresh <- 7 # days (allowable gap in transmission)
# access to the Movebank data
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633, 
             1562253659, 9648615, 2204484313, 10449318, 2106043019)
# visual estimate of the death dates of individuals in these studies
visDODs <- readRDS("/home/hbronnvik/Documents/storkSSFs/visual_estimated_deaths_2023-09-05.rds")

# determine the identities of the nestling birds (remove any care center adults)
info <- lapply(studies, function(x){
  print(x)
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653)
  if("animal_life_stage" %in% colnames(birds)){
    chicks <- birds %>% 
      filter(grepl("juv|chick|nestling", animal_life_stage, ignore.case = T) & grepl("release", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier)
    juv <- birds %>% 
      filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier)
    chicks <- rbind(chicks, juv)
  }else{
    chicks <- birds %>% 
      filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier)
  }
  return(chicks)
}) %>% reduce(rbind)

# remove an empty tag
info <- info %>% 
  filter(animal_local_identifier != "Rheindelta 4 9692 E0709")
# download data
lapply(1:nrow(info), function(x){
  print(x)
  df <- getMovebankLocationData(info$study_id[x], sensorID = 653, 
                                animalName = info$animal_id[x], login = loginStored)
  saveRDS(df, file = paste0("/home/hbronnvik/Documents/storkSSFs/full_data/", info$animal_id[x], "_30082023.rds"))
})

# select the nestling data files
full_files <- list.files("/home/hbronnvik/Documents/storkSSFs/full_data/", pattern = "_30082023.rds", full.names = T)

files_ids <- lapply(full_files, function(x){sub("data//", "", str_split(x, "_")[[1]][2])}) %>% unlist()

# remove errors
start_time <- Sys.time()
clean_locations <- lapply(full_files, function(x){
  # load the data
  ind <- readRDS(x) 
  # clean the data
  locs_df <- ind %>% 
    drop_na(location.long) %>% 
    mutate(index = row_number())
  # remove duplicated locations because they prevent accurate calculations of distance and speed
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),] %>% 
    filter(is.na(height.above.ellipsoid))
  
  locs_df <- locs_df %>% 
    filter(!index %in% doubles$index) 
  
  # warn if a duplicated timestamp contains information other than location (not usually the case)
  if(nrow(locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp),
                                                       "timestamp"]),]) > 0){print("Duplicates containing HAE and DOP values exist.", quote = F)}
  
  # remake doubles in the event of duplicates that hold values
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),]%>% 
    mutate(event = round_date(timestamp, "minute"))
  
  if(nrow(doubles) > 0){
    # if there is more than one instance of duplicates with information
    check <- lapply(unique(round_date(doubles$timestamp, "minute")), function(q){
      # take the duplicates within one minute
      dd <- doubles %>% 
        filter(event == q)
      # determine the last location before the duplicates (point of interest)
      poi <- locs_df %>% 
        filter(timestamp < dd$timestamp[1]) %>% 
        slice(n())
      # find the distance from each duplicate to the poi, even for true points this may increase
      cc <- lapply(1:nrow(dd), function(p){
        d <- distVincentyEllipsoid(c(poi$location.long, poi$location.lat), c(dd$location.long[p], dd$location.lat[p]))
        return(d)
      }) %>% unlist()
      # calculate the ground speeds from each duplicate to the poi
      dd <- dd %>% 
        mutate(dist_from_unique = cc, 
               time_since_unique = as.numeric(difftime(dd$timestamp[1], poi$timestamp, units = "sec")),
               speed_after_unique = dist_from_unique/time_since_unique) %>% 
        filter(speed_after_unique > s_thresh)
      return(dd)
    }) %>% reduce(rbind)
    # filter out the ground speeds higher than reasonable
    locs_df <- locs_df %>% 
      filter(!index %in% check$index)
  }
  
  # down-sample to 15 minutes
  locs_df <- locs_df %>%
    mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
    filter(td >= 300) %>%
    mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>%
    group_by(seq15) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-seq15, -td) %>% 
    # calculate ground speeds
    mutate(distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
           timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           ground_speed_15 = distance/timediff) %>% 
    # remove absurd speeds to clean outliers
    filter(ground_speed_15 < 50)
  
  # add daily metrics
  ind <- locs_df %>% 
    filter(!is.na(location.lat)) %>% 
    mutate(date = date(timestamp)) %>%
    group_by(date) %>% 
    # the Haversine distance between first and last locations of the day
    mutate(daily_dist = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
           # the rhumbline bearing between first and last locations of the day
           daily_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
           compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward")) %>% 
    ungroup() %>% 
    arrange(timestamp)

  print(ind$individual.local.identifier[1])
  return(ind)
}) 
Sys.time() - start_time

# saveRDS(clean_locations, file = paste0("/home/hbronnvik/Documents/storkSSFs/clean_locations_", Sys.Date(), ".rds"))
clean_locations <- readRDS("/home/hbronnvik/Documents/storkSSFs/clean_locations_2023-08-30.rds")

# isolate long-distance movement and label it
ld_locs <- lapply(clean_locations, function(x){tryCatch({
  print(x$individual.id[1])
  df <- x %>% 
    dplyr::select(individual.id, timestamp, location.lat, location.long, daily_dist, daily_direction, ground_speed_15) %>% 
    filter(daily_dist > d_thresh) %>%
    arrange(timestamp) %>%
    mutate(timeLag = as.numeric(timestamp - lag(timestamp), units = "secs"),
           # insert an enormous time lag for the first location to remove the NA
           timeLag = ifelse(is.na(timeLag), 1e5, timeLag),
           newBurst = ifelse(round(timeLag) <= weeks(1), F, T),
           newCluster = ifelse(round(timeLag) <= weeks(6), F, T),
           # take the cumulative sum to act as a unique ID for each burst
           cumu_check_for_burst = cumsum(newBurst),
           cumu_check_for_clust = cumsum(newCluster)) %>% 
    group_by(cumu_check_for_burst) %>% 
    # for the bursts, calculate the time difference between the last and first locations
    mutate(burstLength = length(unique(date(timestamp))),
           # add an ID to each burst
           burstID = as.character(cur_group_id())) %>%
    # remove the bursts that do not meet user-set criteria
    # filter(burstLength > MinBurstLength) %>% 
    ungroup() %>% 
    group_by(cumu_check_for_clust) %>% 
    # for the bursts, calculate the time difference between the last and first locations
    mutate(clusterLength = length(unique(date(timestamp))),
           # add an ID to each burst
           clusterID = as.character(cur_group_id())) %>% 
    ungroup() %>% 
    # clean up the sorting columns
    dplyr::select(-"newBurst", -"cumu_check_for_burst", -"newCluster", -"cumu_check_for_clust")
  
 
  if(nrow(df) > 0){
    class_df <- df %>%
      group_by(burstID) %>% 
      mutate(burst_angle = ifelse(location.lat[1] > location.lat[n()], "south", "north"),
             burst_season = ifelse(month(timestamp[n()]) %in% c(8:11), "late", 
                                   ifelse(month(timestamp[n()]) %in% c(1:6), "early", "change")),
             burst_class = ifelse(burst_season == "early", "spring", 
                                  ifelse(burst_season == "late", "fall",
                                         ifelse(burst_angle == "north" & burst_season == "change", "spring", 
                                                ifelse(burst_angle == "south" & burst_season == "change", "fall", NA)))),
             burst_year = ifelse(month(timestamp[n()]) == 12 & burst_class == "spring", year(years(1) + timestamp[n()]), year(timestamp[n()])),
             trackID = paste(unique(individual.id), burst_class, burst_year, sep = "_")) %>% 
      ungroup()
    
    # remove ragged edges of tracks the bird survived
    # identify ragged edges as bursts that do not have high speed days in them
    ragged <- class_df %>% 
      group_by(burstID) %>% 
      mutate(high_speed = max(daily_dist > 70000)) %>%
      ungroup() %>% 
      filter(high_speed == F)
    # did the bird die?
    if(unique(x$individual.id) %in% visDODs$individual_id){
      dod <- visDODs %>% 
        filter(individual_id == unique(x$individual.id)) %>% 
        dplyr::select(dod) %>% 
        deframe()
    }else{
      dod <- NA
    }
    if(nrow(ragged) > 0){
      if(!is.na(dod)){
        # did the bird die in a ragged edge?
        cutoff <- date(dod) %in% unique(date(ragged$timestamp))
        if(cutoff == T){
          # how ragged?
          messy <- lapply(2:length(class_df$burstID), function(z){
            temp <- class_df %>%
              filter(burstID == z)
            temp_prev <- class_df %>%
              filter(burstID == z-1)
            temp$prev_long <- !unique(temp_prev$burstID) %in% unique(ragged$burstID)
            temp
          }) %>% reduce(rbind)
          # these are preceded by high speed bursts, but are low speed
          messy <- messy %>% 
            filter(prev_long == T & burstID %in% unique(ragged$burstID))
          # did the bird die in a ragged edge adjacent to a high speed burst?
          keep <- date(dod) %in% unique(date(messy$timestamp))
          if(keep == T){
            edge <- unique(class_df$burstID[date(class_df$timestamp) == date(dod)])
            # discard ragged edges that do not contain a death
            class_df <- class_df %>% 
              filter(!burstID %in% ragged$burstID | burstID == edge)
          }
        }else{
          class_df <- class_df %>% 
            filter(!burstID %in% ragged$burstID)
        }
      }else{
        class_df <- class_df %>% 
          filter(!burstID %in% ragged$burstID)
      }
    }
    
    if(nrow(class_df) > 0){
      tracks <- class_df %>%
        group_by(trackID) %>% 
        slice(1, n())
      tracks <- lapply(split(tracks, tracks$trackID), function(c){
        ID <- unique(c$trackID)
        track <- x %>% 
          filter(between(timestamp, c$timestamp[1], c$timestamp[2])) %>% 
          mutate(trackID = ID)
        track
      }) %>% rbindlist()
    
      suppressMessages(x <- x %>% left_join(tracks))
      return(x)
    }
  }
}, error = function(e){
  print(geterrmessage())
})
})
# remove empty items (birds that do not migrate)
migration_locations <- ld_locs[lapply(ld_locs, length) > 1]
# compress the data into a single frame
migration_locations <- rbindlist(migration_locations, fill = T)
# remove all the data that are not classified as migration
migration_locations <- migration_locations %>% 
  filter(!is.na(trackID))
# get the track IDs from the migratory birds
rs_ids <- migration_locations %>% 
  group_by(individual.id, trackID) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, trackID) %>% 
  ungroup()
# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
  mutate(season = ifelse(grepl("fall", trackID), "fall", "spring")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()
# save these numbers before filtering migrations
# saveRDS(meta, file = paste0("/home/hbronnvik/Documents/storkSSFs/metadata_migration_tracks_", Sys.Date(),".rds"))

# add the number of journeys
migration_locations <- migration_locations %>% 
  left_join(meta)

# find any animals that slipped through the cracks for visual estimates of DOD
lost <- migration_locations %>% 
  group_by(individual.id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(individual.id) %>% 
  filter(!individual.id %in% visDODs$individual_id)

# define birds that took eastern routes as ones that are ever east of 16.5 longitude (East Germany)
eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 16.5])
# remove the eastern birds
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% eastern_birds) %>% 
  group_by(trackID) %>% 
  mutate(track_displacement = abs(location.lat[1] - location.lat[n()])) %>% 
  ungroup() 

# write out the file to use for social information estimates
# this includes all western migrants regardless of success or the length or 
# transmission rates of the tracks, but holds only locations that are during migration
# 397 individuals, 929 tracks
# saveRDS(migration_locations, file = paste0("/home/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_", Sys.Date(),".rds"))

# now filter for the route selection analysis:
migration_locations <- migration_locations %>%
  # use only the southwest German birds that travel more than 3 degrees N/S
  filter(study.id %in% c(173641633,1176017658,21231406,76367850,212096177,24442409) & track_displacement > l_thresh)

# the 2023 fall migrations were not complete yet, but might have started, when the data were retrieved
# remove these
migration_locations <- migration_locations %>% 
  filter(!grepl("fall_2023", trackID))

# define the tracks of birds for which data transmission errors prevent determining the track target
# some are also birds with only one or two locations per day, which precludes the possibility
# of accurate calculations of daily displacement and therefore of migration speeds, which should have been
# removed in the data cleaning step but were missed
safe <- migration_locations %>% 
  filter(trackID == "80439771_fall_2015")
ind_gaps <- c(1578604607, 23460528, 80436240, 1576790450, 80439771, 504440981)
track_gaps <- c("219392149_spring_2016", "293986322_spring_2018", "1176046609_spring_2022",
                "909029794_fall_2014", "1576782003_spring_2022", "173659608_spring_2018",
                "78031713_fall_2022", "23463463_spring_2018", "78031713_fall_2022", "80438400_spring_2019",
                "909029800_fall_2016", "1576782003_spring_2022", "23463463_fall_2017")

# 288 individuals, 690 tracks
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% ind_gaps & !trackID %in% track_gaps) %>% 
  # add the one useful track from individual 1178289602
  rbind(safe)

# find the final data addition for the remaining birds
info <- lapply(studies, function(x){
  df <- getMovebankAnimals(x, loginStored) %>% 
    filter(number_of_events > 0) %>% 
    mutate(timestamp_end = as.POSIXct(sub("\\.000", "", timestamp_end), tz = "UTC"),
           study_id = x) %>% 
    filter(sensor_type_id == 653)
}) %>% 
  reduce(rbind) %>% 
  filter(individual_id %in% migration_locations$individual.id)

# find the success/failure of each migration track
migration_locations <- lapply(split(migration_locations, migration_locations$trackID), function(x){
  print(x$trackID[1])
  # take the migratory route
  track <- x %>% 
    arrange(timestamp) 
  # if the last GPS time point is greater than the time the data were downloaded, use the download date
  t_end <- ifelse(info$timestamp_end[info$individual_id == unique(x$individual.id)] > "2023-08-30", 
                  "2023-08-30 00:00:00", info$timestamp_end[info$individual_id == unique(x$individual.id)])
  # take either the confirmed DOD or the last transmitted time stamp
  loss_time <- visDODs %>% 
    filter(individual_id == unique(x$individual.id) | local_identifier == unique(x$individual.local.identifier)) %>% 
    mutate(dod = as.POSIXct(ifelse(is.na(dod), t_end, dod), tz = "UTC", origin = "1970-01-01")) %>% 
    dplyr::select(dod) %>% 
    deframe()
  # compare the DOD to the end of the migration
  loss <- max(track$timestamp) > loss_time - days(s_thresh)
  # add a column containing the outcome of the migratory track
  x <- x %>% 
    mutate(track_status = ifelse(loss == T & trackID == unique(track$trackID), "incomplete", "complete"))
  return(x)
}) %>% reduce(rbind)

# finally, save out the complete, long-distance migrations of the 158 birds from southwest Germany
saveRDS(migration_locations, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40km70km30daySpeed_", Sys.Date(),".rds"))
