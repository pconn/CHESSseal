### output CHESS data to .csv files for archiving to axiom repository / calculate stuff for medadata


library("ChukchiPolarBear")
library(sf)
data(Chukchi_PB_data)

n_t = length(Grid_list)
n_s = nrow(Grid_list[[1]])

Cols_keep = c("cell","area_km","dist_land","land_cover","sea_ice","fast_ice","rsf")
Grid_df = Grid_list[[1]][,Cols_keep]
Grid_df$jday = 98

Centroids = st_coordinates(st_transform(st_centroid(Grid_list[[1]]),4326))
Grid_df$latitude = Centroids[,"Y"]
Grid_df$longitude = Centroids[,"X"]
Grid_df = st_drop_geometry(Grid_df)

#need to add snow melt data from CHESS paper
load('c:/users/paul.conn/git/CHESS/CHESS_data_Mar2021.RData')

for(iday in 2:n_t){
  Tmp_df = Grid_list[[iday]][,Cols_keep]
  Tmp_df$latitude=Centroids[,"Y"]
  Tmp_df$longitude = Centroids[,"X"]
  Tmp_df$jday = iday+97
  Grid_df = rbind(Grid_df,st_drop_geometry(Tmp_df))  
}
Grid_df$snow_depth = CHESS_data$X_s[,"snow_depth"]
colnames(Grid_df)[1]="cell"
colnames(Grid_df)[2]="area_km_sq"

write.csv(Grid_df,file='CHESS_grid_axiom.csv',row.names=FALSE)

#bounding box 


### Russian polar bear data

#note that the process used in coming up with C_i_ru for the analysis was slightly different
#than here but includes the same animals.  report everything in Alaska time so we're working with
#the same "day"

pb_data <- read.csv("pb_data.csv")
pb_data$dt <- apply(pb_data[,c(2,3)],1,paste,collapse=' ')
d1<-strptime(pb_data$dt, format="%m/%d/%Y %H:%M:%S", tz = "Asia/Anadyr")
pb_data$correct_dt <- as.POSIXct(d1, tz="Asia/Anadyr")  
attributes(pb_data$correct_dt)$tzone <- "US/Alaska" 
pb_data$jday = floor(julian(pb_data[,"correct_dt"], origin = as.POSIXct("2016-01-01", tz = "US/Alaska")))

#assign grid cell
obs_sf = st_transform(st_as_sf(pb_data,coords=c("Lon","Lat"),crs=4326,agr="constant"),st_crs(Grid_list[[1]]))
Inter = st_intersects(obs_sf,Grid_list[[1]],sparse=F)
pb_data$cell = 0
for(irow in 1:nrow(pb_data)){
  pb_data$cell[irow]=as.numeric(Grid_list[[1]][which(Inter[irow,]==1),"cell"])[1]
}
pb_data$jday = as.numeric(pb_data$jday)

Cols_keep = c("cell","jday","Adults","Cubs","Dist_m")
pb_data_axiom = pb_data[,Cols_keep]

write.csv(pb_data_axiom,file="pb_encounters_russia.csv")


### US polar bear data
CellDay_on = c("457 35","1023 40","415 7")
CellDay_off = c("1229 15","1309 9","1341 9","1342 9")
Unique_ct = paste(Cell_lookup_table[Area_table[,"Grid_ID"]],Area_table[,"Day"])
Which = rep(0,7)
for(i in 1:3)Which[i]=which(Unique_ct==CellDay_on[i])
for(i in 1:4)Which[i+3]=which(Unique_ct==CellDay_off[i])

pb_us_encounters <- data.frame("cell"=Grid_list[[1]]$cell[Which],
                               "IR_detect"=c(1,1,1,0,0,0,0),
                               "jday" = 97+c(35,40,7,15,9,9,9),
                               "Adults" = rep(1,7),
                               "Cubs" = c(2,rep(0,6)))
write.csv(pb_us_encounters,file="pb_encounters_us.csv")                               


#### u.s. seals, polar bear tracks, and thinning proportions

load("CHESS_us_data_TMB_Mar2021.RData")
load("tracks_filtered.Rda")

Effort_us = data.frame("cell"=Grid_list[[1]]$cell[(Data_chess_us$S_i+1) %% n_s]) #recall S_i index starts at 0 for TMB
n_i = nrow(Effort_us)
Effort_us$jday = 96 + ceiling((Data_chess_us$S_i+1)/n_s)  
Effort_us$area_photo_km2 = as.numeric(Data_chess_us$P_i*st_area(Grid_list[[1]][1,])/1000000)
Effort_us = cbind(Effort_us,data.frame(Data_chess_us$C_i))
colnames(Effort_us)[4:6]=c("bd_count","rd_count","unk_count")
Effort_us$bd_thin = plogis(Data_chess_us$thin_mu_logit[1:n_i])
Effort_us$rd_thin = plogis(Data_chess_us$thin_mu_logit[(n_i+1):(2*n_i)])

# #calculate number of photographs, number with tracks for each unique time & grid cell surveyed

tracks_df = st_drop_geometry(tracks)
jday = floor(julian(tracks_df[,"correct_dt"], origin = as.POSIXct("2016-01-01", tz = "US/Alaska")))
tracks$jday = jday
tracks$day = as.numeric(jday - min(jday))+1
tracks$cell_time = paste(Cell_lookup_table[tracks_df$GridID],tracks$day)
tracks$I_track = as.numeric(tracks_df[,"detect_bear_track"]=="Y")
tracks_df = st_drop_geometry(tracks)

Inter = st_intersects(tracks,Grid_list[[1]],sparse=F)
tracks_df$cell = 0
for(irow in 1:nrow(tracks)){
  tracks_df$cell[irow]=as.numeric(Grid_list[[1]][which(Inter[irow,]==1),"cell"])[1]
}
tracks_df$jday = as.numeric(tracks_df$jday)
tracks_df$cell_jday = apply(tracks_df[,c("cell","jday")],1,paste,collapse=' ')
ID_effort = apply(Effort_us[,c("cell","jday")],1,paste,collapse=' ')
N_i = T_i = rep(0,n_i)
for(irow in 1:nrow(Effort_us)){
  Cur_which = which(tracks_df$cell_jday==ID_effort[irow])
  if(length(Cur_which)>0){
    N_i[irow] = nrow(tracks_df[Cur_which,])
    T_i[irow] = sum(tracks_df[Cur_which,"I_track"])
  }
}

Effort_us$n_photos = N_i
Effort_us$n_photos_pb_tracks = T_i

write.csv(Effort_us,file="effort_us.csv")

#load('CHESS_data_Mar2021.RData')




### Russian CHESS effort
load('CHESS_seal_data_no_area_filter.RData')

Grid_df = st_drop_geometry(Grid_list[[1]])

Effort_rus = data.frame("cell"=Grid_list[[1]]$cell[(CHESS_data$S_rus_i+1) %% n_s])
n_i = nrow(Effort_rus)
Effort_rus$jday = 96 + ceiling((CHESS_data$S_rus_i+1)/n_s)  
Effort_rus$effective_area_surveyed_km2 = as.numeric(CHESS_data$P_rus_i*st_area(Grid_list[[1]][1,])/1000000)
Effort_rus = cbind(Effort_rus,data.frame(CHESS_data$C_rus))
colnames(Effort_rus)[4:6]=c("bd_count","rd_count","unk_count")
Effort_rus$no_photo_count = CHESS_data$Nophoto_i
Effort_rus$bd_thin = plogis(CHESS_data$Thin_mu_rus_i[1:n_i])
Effort_rus$rd_thin = plogis(CHESS_data$Thin_mu_rus_i[n_i + c(1:n_i)])
sum(Effort_rus$cell != pb_ru$cell) #order is same
Effort_rus$n_tracks = pb_ru$N_tracks

write.csv(Effort_rus,file="effort_russia.csv")