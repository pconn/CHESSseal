###  Attach environmental covariates to "big" (US + Russia) CHESS grid

install_pkg <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# Install libraries ----------------------------------------------
install_pkg("RPostgreSQL")
install_pkg("sf")
install_pkg("devtools")
#install_version("crawl",version="1.4-1")
#install_pkg("crawl")
install_pkg("dplyr")
install_pkg("sp")
install_pkg("rgeos")
install_pkg("automap")
install_pkg("rpostgis")
install_pkg("raster")
install_pkg("fasterize")
install_pkg("rgdal")

#function to drop the geometry of an st object (easier to work w/ the data.frame)
st_drop_geometry <- function(x) {
  if(inherits(x,"sf")) {
    x <- st_set_geometry(x, NULL)
    class(x) <- 'data.frame'
  }
  return(x)
}

# Run code -------------------------------------------------------
# Extract data from DB ------------------------------------------------------------------
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))
grid_sf_CHESS <- sf::st_read(con, 
                       query = "SELECT * FROM base.geo_analysis_grid_no_polarbear", 
                       geometry_column = "geom")
grid_sf_CHESS = grid_sf_CHESS[grid_sf_CHESS$cell>20000 & grid_sf_CHESS$cell<25000,]




effort <- sf::st_read(con, 
                         query = "SELECT DISTINCT flight_num, camera_loc, image_name, correct_dt, effort_type, gga_alt, geom 
                         FROM surv_chess.tbl_effort_raw WHERE effort_type = \'ON\'", 
                         geometry_column = "geom")
effort$gga_alt <- effort$gga_alt * 3.28084
effort <- st_transform(effort,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

bears <- sf::st_read(con,
                        query = "SELECT *
                        FROM surv_chess.geo_polar_bear
                        WHERE effort_type = \'ON\'
                        AND detection_skeyes = \'Y\'
                        AND hotspot_type <> \'Duplicate\'",
                       geometry_column = "geom")
bears<- st_transform(bears,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


bears_off <- sf::st_read(con,
                            query = "SELECT *
                            FROM surv_chess.geo_polar_bear
                            WHERE (effort_type = \'BEAR\'
                            OR (detection_skeyes = \'N\' AND effort_type = \'ON\'))
                            AND hotspot_type <> \'Duplicate\'",
                            geometry_column = "geom")
bears_off = bears_off[c(1,11,12,13),]  #take out cubs (their mother was detected on effort) and bears at the point (polar_bear_id = 7-18 were at bone site)
bears_off<- st_transform(bears_off,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


tracks <- sf::st_read(con,
                         query = "SELECT image_name, correct_dt, CASE WHEN detect_bear_track IS NULL THEN \'N\' ELSE detect_bear_track END, geom 
                         FROM surv_chess.tbl_unfilt_detect_beartrk 
                         RIGHT JOIN surv_chess.tbl_unfilt_detect_images 
                         USING (unfilt_image) 
                         INNER JOIN surv_chess.tbl_effort_raw
                         ON unfilt_dt = effort_dt",
                         geometry_column = "geom")
tracks<- st_transform(tracks,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
#limit tracks on 'on effort'
effort_df = st_drop_geometry(effort)
tracks_df = st_drop_geometry(tracks)
Which_on =  tracks_df[,"image_name"] %in% effort_df[,"image_name"]
tracks<- tracks[Which_on,]

fastice_dates <- c("20160407", "20160414", "20160421", "20160428", "20160512", "20160519", "20160526")
for (i in 1:length(fastice_dates)){
  assign(paste("fastice_", fastice_dates[i], sep = ""),
         sf::st_read(con,
                        query = paste("SELECT * FROM surv_chess.geo_fastice_", fastice_dates[i], sep = ""),
                        geometry_column = "geom"))
}



effort_df <- st_drop_geometry(effort)
Ecoords <- st_coordinates(effort)
effort_df = cbind(effort_df,Ecoords)
tracks_df <- st_drop_geometry(tracks)

Cell_lookup <- st_intersects(st_centroid(fastice_20160407),grid_sf_CHESS,sparse=FALSE)  
sum(diag(Cell_lookup)) #Looks like ordering of grid_sf and fast_ice are the same despite different cell #s

#NB: 4/07 and 4/14 had issues and Cynthia and Erin R. did not recommend them for use
grid_sf_CHESS$fi_20160421 =  fastice_20160421$fi_pro
grid_sf_CHESS$fi_20160428 = fastice_20160428$fi_pro
grid_sf_CHESS$fi_20160512 = fastice_20160512$fi_pro
grid_sf_CHESS$fi_20160519 = fastice_20160519$fi_pro
grid_sf_CHESS$fi_20160526 = fastice_20160526$fi_pro



#attach distance from mainland, distance from land
AK = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ak_dcw.shp')
AK = st_transform(AK,st_crs(grid_sf_CHESS))
Russia = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/russia_dcw.shp')
Russia = st_transform(Russia,st_crs(grid_sf_CHESS))

Area_AK = st_area(AK)
AK_mainland=AK[which(Area_AK==max(Area_AK)),2]
Area_RU = st_area(Russia)
Rus_mainland=Russia[which(Area_RU==max(Area_RU)),2]

Dist1 = as.numeric(st_distance(grid_sf_CHESS,AK_mainland))
Dist2 = as.numeric(st_distance(grid_sf_CHESS,Rus_mainland))
grid_sf_CHESS$dist_mainland = rep(0,length(Dist1))
for(i in 1:length(Dist1))grid_sf_CHESS$dist_mainland[i]=min(Dist1[i],Dist2[i])

Dist_AK=apply(st_distance(AK,grid_sf_CHESS,byid=TRUE),2,'min')
Dist_Rus=apply(st_distance(Russia,grid_sf_CHESS,byid=TRUE),2,'min')
Dist_land=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
for(i in 1:length(Dist_AK))grid_sf_CHESS$dist_land[i]=min(Dist_AK[i],Dist_Rus[i])

#calculate proportion area covered by land and remove cells with >99% of area on land
Land=st_union(st_union(AK),st_union(Russia))
n_cells=nrow(grid_sf_CHESS)
Land.area=rep(0,n_cells)
I.intersect = st_intersects(grid_sf_CHESS,Land)
for(icell in 1:n_cells){
  if(length(I.intersect[[icell]])>0)Land.area[icell]=st_area(st_intersection(grid_sf_CHESS[icell,],Land))
}
grid_sf_CHESS$land_cover=Land.area/628381060
grid_sf_CHESS=grid_sf_CHESS[-which(grid_sf_CHESS$land_cover>0.99),]

#attach depth
depth <- sf::st_read(con,
                         query = "SELECT *
                         FROM environ.geo_etopo1_sdr_bins", 
                         geometry_column = "geom")
depth = st_transform(depth,st_crs(grid_sf_CHESS))
depth = depth[-which(depth$gridcode>0),]  #take out land
#take out features that don't intersect w/ chess grid
chess_union = st_union(grid_sf_CHESS)
I.grid = st_intersects(depth,chess_union,sparse=FALSE)
depth=depth[I.grid>0,]
depth$mean = apply(cbind(depth$gridcode,depth$meandepth),1,'mean')  #using midpoint between min and max depth for each polygon

grid_sf_CHESS$depth = NA
for(i in 1:nrow(grid_sf_CHESS)){
  Pts = st_sample(grid_sf_CHESS[i,],400,type='regular')  #lay down 400 regular points in each cell
  N.cells = colSums(st_intersects(Pts,depth,sparse=FALSE))
  grid_sf_CHESS$depth[i]=sum(N.cells %*% depth$mean)/sum(N.cells)
} 
grid_sf_CHESS$depth[which(is.na(grid_sf_CHESS$depth))]=0  #in case some cells 400 points were all on land
#note depth truncated at 600m (deeper gets assigned 600m)

### attach snow melt, snow depth covariates
load('c:/users/paul.conn/git/haulout/ringed/snow_melt_grid_Feb2021.Rda')  #compiled w/ "attach_snow_melt_grid.R
grid_sf_CHESS$melt_MDSDA = grid_sf$MDSDA_2016[match(grid_sf_CHESS$cell,grid_sf$cell)]
grid_sf_CHESS$melt_Markus = grid_sf$Markus_2016[match(grid_sf_CHESS$cell,grid_sf$cell)]
grid_sf_CHESS$snow_depth = grid_sf$snow_2016[match(grid_sf_CHESS$cell,grid_sf$cell)]

#attach temperature
load('c:/users/paul.conn/git/haulout/ringed/Temps_cell_year.RDa')
grid_sf_CHESS$spring_temp =  Cell_temp_yr[match(grid_sf_CHESS$cell,grid_sf$cell),12]  #Cell temp yr columns go 2005-2019

# #old code to attach Nikita Platanov's MDSDA melt date covariate
# melt_shape = readOGR(dsn = 'c:/users/paul.conn/git/haulout/data/MDSDA/2016-12-05-MDSDA/mdm2016a.shp')
# melt_sf = st_as_sf(melt_shape)
# #krige values with low coverage
# melt_coords = st_coordinates(st_centroid(melt_sf))
# Which_miss = which(is.na(melt_sf$median))
# pred_loc=SpatialPointsDataFrame(melt_coords[Which_miss,],data.frame(median=matrix(rep(0,length(Which_miss)),ncol=1)),proj4string=CRS(st_crs(grid_sf_CHESS)))
# melt_points=SpatialPointsDataFrame(melt_coords[-Which_miss,],data.frame(median=matrix(melt_sf$median[-Which_miss],ncol=1)),proj4string=CRS(grid_proj))
# krige_out=autoKrige(median~1,input_data=melt_points,new_data=pred_loc)$krige_output[["var1.pred"]] 
# melt_sf$median[Which_miss]=krige_out
# new_grid_large$MDSDA_median_jday = NA
# #match CHESS grid to fuller grid from Nikita
# Which_entry = rep(0,nrow(new_grid_large))
# Distances = st_distance(st_centroid(new_grid_large),st_centroid(melt_sf))
# for(i in 1:nrow(new_grid_large)){
#   new_grid_large$median[i]=melt_sf$median[which(Distances[i,]==min(Distances[i,]))]
# }


#old code to attach time averaged RSF
#rsf = read.delim("./polar bear info to NMFS 30aug17/CS_pbtoNMFS_30aug17.csv",sep="\t")
#rsf_proj = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
#grid_trans = st_transform(new_grid_large,rsf_proj)
#coords_grid = st_coordinates(st_centroid(grid_trans))
#Dists = matrix(0,nrow(coords_grid),nrow(rsf))
#rsf_grid = rep(0,nrow(coords_grid))
#for(igrid in 1:nrow(coords_grid)){
#  Dists[igrid,] = (coords_grid[igrid,1]-rsf[,"x"])^2 + (coords_grid[igrid,2]-rsf[,"y"])^2
#  rsf_grid[igrid] = rsf[which(Dists[igrid,]==min(Dists[igrid,])),"rsf2016"]
#}
#new_grid_large[,"rsf"]=rsf_grid
#plot(new_grid_large[,"rsf"],main="original rsf")
#fill missing cells in by kriging
# library(automap)
# I.missing=is.na(rsf_grid)
# Which.missing = which(I.missing == 1)
# n.na=sum(I.missing)
# na.dm=data.frame(matrix(rep(0,n.na),ncol=1))
# colnames(na.dm)="rsf"
# rsf2=data.frame(matrix(rsf_grid[-Which.missing],ncol=1))
# colnames(rsf2)="rsf"
# coords = coords_grid[-Which.missing,]
# coords.pred = coords_grid[Which.missing,]
# rownames(coords)=c(1:nrow(coords))
# rownames(coords.pred)=c(1:nrow(coords.pred))
# rsf_points2=SpatialPointsDataFrame(coords,rsf2,proj4string=CRS(rsf_proj))
# pred_loc=SpatialPointsDataFrame(coords.pred,na.dm,proj4string=CRS(rsf_proj))
# krige_out=autoKrige(rsf~1,input_data=rsf_points2,new_data=pred_loc)$krige_output[["var1.pred"]] 
# new_grid_large[Which.missing,"rsf"]=krige_out
# plot(new_grid_large[,"rsf"],main="kriged rsf")


### get and attach sea ice data
qry <- "WITH a AS
(SELECT id, cell, ST_Transform(ST_Centroid(geom), 3338) as centroid
FROM base.geo_analysis_grid)
SELECT a.*, fdate, ST_Value(rast, centroid)
FROM a, environ.tbl_sic_cdr_conc
WHERE fdate >= '2016-04-07' AND fdate <= '2016-05-31'"

#ice <- sf::st_read(con,query=qry,sf_column_name="centroid")
ice <- sf::st_read(con,query=qry,geometry_column="centroid")

sea_ice_df = st_drop_geometry(ice)
Unique_days = unique(sea_ice_df[,"fdate"])
grid_df = st_drop_geometry(grid_sf_CHESS)
Grid_list = vector("list",55)
Ice = rep(0,nrow(grid_df))
Problem_cells = Problem_cells2 = rep(0,nrow(grid_df))  #try to figure out which cells often have problems with ice<0, etc.
for(iday in c(1:27,29:55)){
  Cur_ice = sea_ice_df[which(sea_ice_df[,"fdate"]==Unique_days[iday]),]
  for(icell in 1:nrow(grid_df))Ice[icell]=Cur_ice[which(Cur_ice[,"cell"]==grid_df[icell,"cell"]),"st_value"]
  
  if(sum(Ice<0)>0)Problem_cells[which(Ice<0)]=1
}
for(iday in 1:55){
  jday = iday+97  #April 7 = jday 98
  Grid_list[[iday]]=grid_sf_CHESS
  Grid_list[[iday]][,"dmelt_Markus"]= jday-Grid_list[[iday]]$melt_Markus
  Grid_list[[iday]][,"dmelt_MDSDA"]= jday-Grid_list[[iday]]$melt_MDSDA
  Cur_ice = sea_ice_df[which(sea_ice_df[,"fdate"]==Unique_days[iday]),]
  for(icell in 1:nrow(grid_df))Ice[icell]=Cur_ice[which(Cur_ice[,"cell"]==grid_df[icell,"cell"]),"st_value"]
  
  Grid_list[[iday]][,"sea_ice"]= Ice
  if(iday<=19)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160421"]
  if(iday>19 & iday<=30)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160428"]
  if(iday>30 & iday<=40)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160512"]
  if(iday>40 & iday<=47)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160519"]
  if(iday>47)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160526"]
}

Grid_list[[28]]=Grid_list[[27]]
#now, if a cell always has land cover > 0.25, set = NA
grid_df=st_drop_geometry(Grid_list[[1]])
Problem_cells2[which(grid_df[,"land_cover"]>0.25)]=1
# Problem_cells2[which(grid_df[,"sea_ice"]<0.2)]=1
# for(iday in 2:55){
#   grid_df=st_drop_geometry(Grid_list[[iday]])
#   for(icell in 1:nrow(grid_df)){
#     if(Problem_cells2[icell]==1 & grid_df[icell,"sea_ice"]>0.2)Problem_cells2[icell]=0
#   }
# }
Which_problem = which(Problem_cells==1 | Problem_cells2==1)
#now krige cells that either have >25% land or ever have ice<0.0
for(iday in 1:55){
  Grid_list[[iday]][Which_problem,"sea_ice"]=NA
  ice_points=as(st_centroid(Grid_list[[iday]]),'Spatial')
  fit_points = ice_points[-Which_problem,]
  pred_points = ice_points[Which_problem,]
  krige_out=autoKrige(sea_ice~1,input_data=fit_points,new_data=pred_points)$krige_output[["var1.pred"]] 
  if(sum(krige_out<0)>0){
    krige_out[which(krige_out<0)]=0
  }
  Grid_list[[iday]][Which_problem,"sea_ice"]=krige_out
}

#OLD: attach difference in julian day from median melt (MDSDA)
# start of grid = April 7 = jday 97
# for(iday in 1:55){
#   Grid_list[[iday]][,"melt_jday_diff"]=iday+96-grid_sf_CHESS$median
# }



oldwd=getwd()

setwd('c:/users/paul.conn/git/chess/RSF_daily')
library(raster)
library(rgdal)

Fnames = rep("",24+31)  #filenames
day1 = 7
for(it in 1:3){
  Fnames[it]=paste0('pred040',it+6,'2016.tif')
  day1 = day1+1
}
for(it in 4:24){
  Fnames[it]=paste0('pred04',it+6,'2016.tif')
  day1 = day1+1
}
day1 = 1
for(it in 1:9){
  Fnames[24+it]=paste0('pred050',it,'2016.tif')
  day1 = day1+1
}
for(it in 10:31){
  Fnames[24+it]=paste0('pred05',it,'2016.tif')
  day1 = day1+1
}

it=1
rsf_raster = raster(Fnames[it])
rsf_SPDF=as(rsf_raster, "SpatialPixelsDataFrame")
rsf_sf = st_as_sf(rsf_SPDF)
rsf_sf= st_transform(rsf_sf,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
rsf_coords = st_coordinates(st_centroid(rsf_sf))
#match w/ chess grid
Which_entry = rep(0,nrow(grid_sf_CHESS))
Distances = st_distance(st_centroid(grid_sf_CHESS),st_centroid(rsf_sf))
for(i in 1:nrow(grid_sf_CHESS)){
  Which_entry[i]=which(Distances[i,]==min(Distances[i,]))
  Grid_list[[1]]$rsf[i]=as.numeric(rsf_sf[Which_entry[i],])[1]
}
for(it in 2:55){
  rsf_raster = raster(Fnames[it])
  rsf_SPDF=as(rsf_raster, "SpatialPixelsDataFrame")
  rsf_sf = st_as_sf(rsf_SPDF)
  rsf_sf= st_transform(rsf_sf,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  Distances = st_distance(st_centroid(grid_sf_CHESS),st_centroid(rsf_sf))
  for(i in 1:nrow(grid_sf_CHESS)){
    Which_entry[i]=which(Distances[i,]==min(Distances[i,]))
    Grid_list[[it]]$rsf[i]=as.numeric(rsf_sf[Which_entry[i],])[1]
  }
  
}
setwd(oldwd)






save(Grid_list,file="Chess_grid_list_all_Feb2021.RDa")
RPostgreSQL::dbDisconnect(con)
rm(con, install_pkg, fastice_dates,ice,sea_ice_df,Grid_list)

