#######  
#
#    R script to calculate survey effort for CHESS (taking out turns, etc.)
#
#    Paul Conn, 11 Jan 2018
#
#######


#### (1) get data from database using code from S. Hardy
# CHESS: Export polar bear and associated effort data for Paul Conn's power analysis
# S. Hardy, 08NOV2017

# Create functions -----------------------------------------------
# Function to install packages needed
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
install_pkg("crawl")
install_pkg("dplyr")
install_pkg("sp")
install_pkg("rgeos")
install_pkg("automap")
install_pkg("rpostgis")
install_pkg("raster")

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
effort <- sf::st_read_db(con, 
                         query = "SELECT DISTINCT flight_num, camera_loc, image_name, correct_dt, effort_type, gga_alt, geom 
                         FROM surv_chess.tbl_effort_raw WHERE effort_type = \'ON\'", 
                         geom_column = "geom")
effort$gga_alt <- effort$gga_alt * 3.28084
effort <- st_transform(effort,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


# bears <- sf::st_read_db(con,
#                          query = "SELECT flight_num, camera_loc, image_name, correct_dt, species_id, e.geom
#                          FROM surv_chess.geo_polar_bear
#                          INNER JOIN surv_chess.tbl_process p
#                          ON e.effort_dt = p.process_dt_c
#                          WHERE effort_type = \'ON\'
#                          AND species_id = \'Polar Bear\'
#                          AND hotspot_type <> \'Duplicate\'",
#                          geom_column = "geom")

bears <- sf::st_read_db(con,
                       query = "SELECT *
                       FROM surv_chess.geo_polar_bear
                       WHERE effort_type = \'ON\'
                       AND detection_skeyes = \'Y\'
                       AND hotspot_type <> \'Duplicate\'",
                       geom_column = "geom")
bears<- st_transform(bears,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


bears_off <- sf::st_read_db(con,
                            query = "SELECT *
                            FROM surv_chess.geo_polar_bear
                            WHERE (effort_type = \'BEAR\'
                            OR (detection_skeyes = \'N\' AND effort_type = \'ON\'))
                            AND hotspot_type <> \'Duplicate\'",
                            geom_column = "geom")
bears_off = bears_off[c(1,11,12,13),]  #take out cubs (their mother was detected on effort) and bears at the point (polar_bear_id = 7-18 were at bone site)
bears_off<- st_transform(bears_off,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


tracks <- sf::st_read_db(con,
                         query = "SELECT image_name, correct_dt, CASE WHEN detect_bear_track IS NULL THEN \'N\' ELSE detect_bear_track END, geom 
                         FROM surv_chess.tbl_unfilt_detect_rate_beartrk 
                         RIGHT JOIN surv_chess.tbl_unfilt_detect_rate_images 
                         USING (unfilt_image) 
                         INNER JOIN surv_chess.tbl_effort_raw
                         ON unfilt_dt = effort_dt",
                         geom_column = "geom")
tracks<- st_transform(tracks,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
#limit tracks on 'on effort'
effort_df = st_drop_geometry(effort)
tracks_df = st_drop_geometry(tracks)
Which_on =  tracks_df[,"image_name"] %in% effort_df[,"image_name"]
tracks<- tracks[Which_on,]



#qry <- "select * from
#        pep.environ.tbl_sic_cdr_conc_2013"
#ice <- sf::st_read_db(con,query=qry,geom_column="geom")
        #where
        #se_mu_x < 50000
        #and se_mu_y < 50000
        #and percent_dry IS NOT NULL"
#library(raster)
#library(rgdal)

#dsn="PG:dbname='plots' host=localhost user='test' password='test' port=5432 schema='gisdata' table='map' mode=2"

#ras <- readGDAL(dsn) # Get your file as SpatialGridDataFrame


# Cynthia and Erin do no recommend using 20160407 or 20160414 because of questionable sea ice assignments in the raw data from these dates
fastice_dates <- c("20160407", "20160414", "20160421", "20160428", "20160512", "20160519", "20160526")
for (i in 1:length(fastice_dates)){
  assign(paste("fastice_", fastice_dates[i], sep = ""),
         sf::st_read_db(con,
                        query = paste("SELECT * FROM surv_chess.geo_fastice_", fastice_dates[i], sep = ""),
                        geom_column = "geom"))
}



#####  2) construct US chukchi analysis grid from 'full' chukchi grid



effort_df <- st_drop_geometry(effort)
Ecoords <- st_coordinates(effort)
effort_df = cbind(effort_df,Ecoords)

tracks_df <- st_drop_geometry(tracks)

#load grid from Chukchi power analysis for use in restricting new grid to areas east of the 
#EEZ
load('c:/users/paul.conn/git/chukchipower/Grid_chukchi.RDa')
old_grid <- st_as_sf(Grid_chukchi$Grid) 
rm(Grid_chukchi)
old_grid <- st_union(old_grid,by_feature=FALSE)
#reproject
old_grid <- st_transform(old_grid,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

new_grid_large <- fastice_20160407 
#new_grid_large[,"fi_20160407"] = fastice_20160407[,"fi_pro"]  #not recommended for use
#new_grid_large[,"fi_20160414"] = fastice_20160414[,"fi_pro"]
new_grid_large[,"fi_20160421"] = fastice_20160421[,"fi_pro"]
new_grid_large[,"fi_20160428"] = fastice_20160428[,"fi_pro"]
new_grid_large[,"fi_20160512"] = fastice_20160512[,"fi_pro"]
new_grid_large[,"fi_20160519"] = fastice_20160519[,"fi_pro"]
new_grid_large[,"fi_20160526"] = fastice_20160526[,"fi_pro"]

Intersects <- st_intersects(new_grid_large,old_grid)
I.intersect = rep(0,length(Intersects))
for(i in 1:length(Intersects))I.intersect[i]=length(Intersects[[i]])
new_grid_US <- new_grid_large[which(I.intersect==1),]

### take out cells that appear to have less than 50% of area on side of US EEZ
AK = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ak_dcw.shp')
AK = st_transform(AK,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
Russia = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/russia_dcw.shp')
Russia = st_transform(Russia,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
EEZ = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/alaska_EEZ_line_edit.shp')
EEZ = st_transform(EEZ,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
plot(new_grid_US[,1],col='white')
plot(EEZ[,1],add=TRUE)
plot(new_grid_US[c(1,2,11,19,27,37,48,61,89,105,122,141,161,162,182,183,203,222,240,256,257,271,272,286,287,300,301,313,325,336,337,348,349,360,370,380,381,387,393,394,400,401,410,420,431,432,447,448,464,480,495,510,511,517,518,524,525),1],col='blue',add=TRUE)

lat_long = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ne_10m_graticules_1/ne_10m_graticules_1.shp')
lat_long = st_transform(lat_long,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
plot(lat_long[516,1],add=TRUE)

new_grid_US=new_grid_US[-c(1,2,11,19,27,37,48,61,89,105,122,141,161,162,182,183,203,222,240,256,257,271,272,286,287,300,301,313,325,336,337,348,349,360,370,380,381,387,393,394,400,401,410,420,431,432,447,448,464,480,495,510,511,517,518,524,525),]

#calculate and attach distance from land
#Attach distance to mainland for each cell (distance from grid cell centroid to landscape polygon)
Area_AK = st_area(AK)
AK_mainland=AK[which(Area_AK==max(Area_AK)),2]
Area_RU = st_area(Russia)
Rus_mainland=Russia[which(Area_RU==max(Area_RU)),2]

Dist1 = as.numeric(st_distance(new_grid_US,AK_mainland))
Dist2 = as.numeric(st_distance(new_grid_US,Rus_mainland))
new_grid_US$dist_mainland = rep(0,length(Dist1))
for(i in 1:length(Dist1))new_grid_US$dist_mainland[i]=min(Dist1[i],Dist2[i])

#attach grid centroids
Centroids=st_coordinates(st_centroid(new_grid_US))
new_grid_US$easting=Centroids[,"X"]
new_grid_US$northing=Centroids[,"Y"]

#calculate proportion area covered by land and remove cells with >99% of area on land
Land=st_union(st_union(AK),st_union(Russia))
n_cells=nrow(new_grid_US)
Land.area=rep(0,n_cells)
I.intersect = st_intersects(new_grid_US,Land)
for(icell in 1:n_cells){
  if(length(I.intersect[[icell]])>0)Land.area[icell]=st_area(st_intersection(new_grid_US[icell,],Land))
}
new_grid_US$land_cover=Land.area/628381060
new_grid_US=new_grid_US[-which(new_grid_US$land_cover>0.99),]

#attach polar bear rsf to grid
rsf = read.delim("./polar bear info to NMFS 30aug17/CS_pbtoNMFS_30aug17.csv",sep="\t")

#reproject US chess grid to RSF scale 
rsf_proj = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
grid_trans = st_transform(new_grid_US,rsf_proj)
coords_grid = st_coordinates(st_centroid(grid_trans))
Dists = matrix(0,nrow(coords_grid),nrow(rsf))
rsf_grid = rep(0,nrow(coords_grid))
for(igrid in 1:nrow(coords_grid)){
  Dists[igrid,] = (coords_grid[igrid,1]-rsf[,"x"])^2 + (coords_grid[igrid,2]-rsf[,"y"])^2
  rsf_grid[igrid] = rsf[which(Dists[igrid,]==min(Dists[igrid,])),"rsf2016"]
}

new_grid_US[,"rsf"]=rsf_grid
plot(new_grid_US[,"rsf"],main="original rsf")

#fill missing cells in by kriging
library(automap)
I.missing=is.na(rsf_grid)
Which.missing = which(I.missing == 1)
n.na=sum(I.missing)
na.dm=data.frame(matrix(rep(0,n.na),ncol=1))
colnames(na.dm)="rsf"
rsf2=data.frame(matrix(rsf_grid[-Which.missing],ncol=1))
colnames(rsf2)="rsf"
coords = coords_grid[-Which.missing,]
coords.pred = coords_grid[Which.missing,]
rownames(coords)=c(1:nrow(coords))
rownames(coords.pred)=c(1:nrow(coords.pred))
rsf_points2=SpatialPointsDataFrame(coords,rsf2,proj4string=CRS(rsf_proj))
pred_loc=SpatialPointsDataFrame(coords.pred,na.dm,proj4string=CRS(rsf_proj))
krige_out=autoKrige(rsf~1,input_data=rsf_points2,new_data=pred_loc)$krige_output[["var1.pred"]] 
new_grid_US[Which.missing,"rsf"]=krige_out
plot(new_grid_US[,"rsf"],main="kriged rsf")

### get and attach sea ice data
qry <- "WITH a AS
(SELECT id, cell, ST_Transform(ST_Centroid(geom), 3338) as centroid
  FROM base.geo_analysis_grid)
SELECT a.*, fdate, ST_Value(rast, centroid)
FROM a, environ.tbl_sic_cdr_conc
WHERE fdate >= '2016-04-07' AND fdate <= '2016-05-31'"

ice <- sf::st_read_db(con,query=qry,geom_column="centroid")
sea_ice_df = st_drop_geometry(ice)
Unique_days = unique(sea_ice_df[,"fdate"])
grid_df = st_drop_geometry(new_grid_US)
Grid_list = vector("list",55)
Ice = rep(0,nrow(grid_df))
Problem_cells = Problem_cells2 = rep(0,nrow(grid_df))  #try to figure out which cells often have problems with ice<0, etc.
for(iday in c(1:27,29:55)){
  Cur_ice = sea_ice_df[which(sea_ice_df[,"fdate"]==Unique_days[iday]),]
  for(icell in 1:nrow(grid_df))Ice[icell]=Cur_ice[which(Cur_ice[,"cell"]==grid_df[icell,"cell"]),"st_value"]
  
  if(sum(Ice<0)>0)Problem_cells[which(Ice<0)]=1
}
for(iday in 1:55){
  Grid_list[[iday]]=new_grid_US
  Cur_ice = sea_ice_df[which(sea_ice_df[,"fdate"]==Unique_days[iday]),]
  for(icell in 1:nrow(grid_df))Ice[icell]=Cur_ice[which(Cur_ice[,"cell"]==grid_df[icell,"cell"]),"st_value"]
  
  Grid_list[[iday]][,"sea_ice"]= Ice
  if(iday<19)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160421"]
  if(iday>19 & iday<30)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160428"]
  if(iday>30 & iday<40)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160512"]
  if(iday>40 & iday<47)Grid_list[[iday]][,"fast_ice"]= grid_df[,"fi_20160519"]
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



save(Grid_list,file="US_chess_grid_list.RDa")
RPostgreSQL::dbDisconnect(con)
rm(con, install_pkg, fastice_dates,ice,sea_ice_df,Grid_list)

###### 3) intersect effort with new chukchi grid

#remove images on land or outside grid  - this step takes a long time!!!!
Grid_union = st_union(new_grid_US)
Intersect.land = st_intersects(effort,Land)
Intersect.grid = st_intersects(effort,Grid_union)

IL_vec = unlist(lapply(Intersect.land,"length"))
Which_IL = which(IL_vec>0)

IG_vec = unlist(lapply(Intersect.grid,"length"))
Which_OutGrid = which(IG_vec==0)

Which_out = unique(c(Which_IL,Which_OutGrid))
plot(effort[Which_out,],add=TRUE)  #plot which photographs will be removed. (8935 photos)

Effort_clean = effort[-Which_out,]

## remove tracks or bear sightings that occur on land or outside grid
Intersect.land = st_intersects(tracks,Land)
Intersect.grid = st_intersects(tracks,Grid_union)
IL_vec = unlist(lapply(Intersect.land,"length"))
Which_IL = which(IL_vec>0)
IG_vec = unlist(lapply(Intersect.grid,"length"))
Which_OutGrid = which(IG_vec==0)
Which_out = unique(c(Which_IL,Which_OutGrid))
tracks = tracks[-Which_out,]

Intersect.land = st_intersects(bears,Land)
Intersect.grid = st_intersects(bears,Grid_union)  #all bears okay

Intersect.land = st_intersects(bears_off,Land)
Intersect.grid = st_intersects(bears_off,Grid_union)  #all bears okay


#get a CHESS grid ID ("GridID") for each on effort photograph by overlaying points on grid
#run in blocks of 1000 (slightly faster than 10000, MUCH faster than 1)
Which_cell = rep(0,nrow(Effort_clean))
for(i in 1:826){
  Tmp_range = c(((i-1)*1000+1):(i*1000))
  Tmp_mat = unlist(st_intersects(Effort_clean[Tmp_range,],new_grid_US))
  Which_cell[Tmp_range]=Tmp_mat
}
Tmp_range=c(826001:826490)
Tmp_mat = unlist(st_intersects(Effort_clean[Tmp_range,],new_grid_US))
Which_cell[Tmp_range]=Tmp_mat
Effort_clean[["GridID"]]=Which_cell

save(Effort_clean,file='Effort_clean.Rda')

# attach gridID to tracks, bears
Which_cell = unlist(st_intersects(tracks,new_grid_US))
tracks[["GridID"]]=Which_cell
Which_cell = unlist(st_intersects(bears,new_grid_US))
bears[["GridID"]]=Which_cell
Which_cell = unlist(st_intersects(bears_off,new_grid_US))
bears_off[["GridID"]]=Which_cell

#reorder points within flightID by date / time
n_photos=nrow(Effort_clean)
#Order_all = rep(0,n_photos)
FlightSide_ID = paste(Effort_clean$flight_num,Effort_clean$camera_loc)
Effort_clean$FlightSide_ID = FlightSide_ID
FlightSide_unique_IDs = unique(FlightSide_ID)
#Coords=matrix(0,n_photos,2)
#Data_ordered = st_drop_geometry(Effort_clean)
#Data_ordered = Effort_clean


Data_grouped = Effort_clean %>% group_by(FlightSide_ID)
Data_ordered = Data_grouped %>% arrange(correct_dt,.by_group=TRUE)
Data_ordered = Data_ordered %>% ungroup(Data_ordered)
Data_ordered = st_as_sf(Data_ordered)
#cur_pl=1
#for(ifl in 1:length(FlightSide_unique_IDs)){
#  Which = which(FlightSide_ID == FlightSide_unique_IDs[ifl])
#  Tmp=Effort_clean[Which,]
#  Order = order(Effort_clean[Which,]$correct_dt)
  #Coords[cur_pl:(cur_pl+length(Which)-1),]=st_coordinates(Tmp[Order,])
  #Data_ordered[cur_pl:(cur_pl+length(Which)-1),]=st_drop_geometry(Tmp[Order,])
  #Order_all[Which]=cur_pl+Order-1
#  st_geometry(Data_ordered[cur_pl:(cur_pl+length(Which)-1),])=st_geometry(Tmp[Order,])
#  cur_pl=cur_pl+length(Which)
#}


#produce some unique identifiers for points obtained for each flightID, camera, and pass through a grid cell
Data_df = st_drop_geometry(Data_ordered)
FlightCellSide_ID = paste(Data_df[,"flight_num"],Data_df[,"GridID"],Data_df[,"camera_loc"])
FlightCell_ID = paste(Data_df[,"flight_num"],Data_df[,"GridID"])
FlightCellSide_unique_IDs = unique(FlightCellSide_ID)

#for each identifier, break flights into "legs" based on number of segments >5 seconds apart [so that different aircraft interpolations can be done for each]
for(iid in 1:length(FlightCellSide_unique_IDs)){
  Which_photos = which(FlightCellSide_ID==FlightCellSide_unique_IDs[iid])
  n_photos=length(Which_photos)
  DT<-Data_df[Which_photos,"correct_dt"]
  DT_diff=difftime(DT[2:length(DT)],DT[1:(length(DT)-1)],units="secs")
  Breaks = which(abs(DT_diff)>5)
  if(length(Breaks)==1){
    FlightCellSide_ID[Which_photos[1:Breaks]]=paste0(FlightCellSide_ID[Which_photos[1:Breaks]],'leg1')
    FlightCellSide_ID[Which_photos[(Breaks+1):length(Which_photos)]]=paste0(FlightCellSide_ID[Which_photos[(Breaks+1):length(Which_photos)]],'leg2')
  }
  if(length(Breaks)>1){
    FlightCellSide_ID[Which_photos[1:Breaks[1]]]=paste0(FlightCellSide_ID[Which_photos[1:Breaks[1]]],'leg1')
    for(ibreak in 2:(length(Breaks))){
      FlightCellSide_ID[Which_photos[(Breaks[ibreak-1]+1):Breaks[ibreak]]]=paste0(FlightCellSide_ID[Which_photos[(Breaks[ibreak-1]+1):Breaks[ibreak]]],paste0('leg',ibreak))
    } 
    FlightCellSide_ID[Which_photos[(Breaks[length(Breaks)]+1):length(Which_photos)]]=paste0(FlightCellSide_ID[Which_photos[(Breaks[length(Breaks)]+1):length(Which_photos)]],paste0('leg',length(Breaks)+1))
  }
}

#for each "leg", base plane x, y position for duplicated x,y coords based on flight track and frequency of camera firings
FlightCellSide_unique_IDs = unique(FlightCellSide_ID) 
Coords=st_coordinates(Data_ordered)
Interp_loc = Coords
for(iid in 1:length(FlightCellSide_unique_IDs)){
  Which_photos = which(FlightCellSide_ID==FlightCellSide_unique_IDs[iid])
  n_photos=length(Which_photos)
  if(n_photos>1 & n_photos<10){
    Tmp_coords = Coords[Which_photos,]
    I_dup = (Tmp_coords[2:n_photos,1]==Tmp_coords[1:(n_photos-1),1] & Tmp_coords[2:n_photos,2]==Tmp_coords[1:(n_photos-1),2])
    Duplicated = as.numeric(which(I_dup==1))+1
    if(length(Duplicated)>0){
      #first.pl = 1
      #last.pl = max(which(I_dup==0))  #first and last non-duplicated
      #number of duplicated sequences
      I.str = rep(0,length(Duplicated))
      if(length(Duplicated)>1){
        I.str[2:length(Duplicated)]=(Duplicated[2:length(Duplicated)]==(Duplicated[1:(length(Duplicated)-1)]+1))
      }
      n.strings = sum(I.str==0)
      Begin=End = which(I.str==0)  #beginning and end of each missing (well, erroneous) gps stream
      for(istr in 1:(n.strings-1))End[istr]=Begin[istr+1]-1
      End[n.strings]=length(Duplicated)
      
      #first, treat case where there is only one record with usable GPS info.
      #in this case, can't interpolate or extrapolate; simply add 150m to each lat and long so that will be adding up the total photo area
      if((n_photos-length(Duplicated))==1){
        Tmp_coords[2:n_photos,1]=Tmp_coords[1,1]+150*(1:length(Duplicated))
        Tmp_coords[2:n_photos,2]=Tmp_coords[1,2]+150*(1:length(Duplicated))
      }
      else{  #>1 record with usable GPS info: interpolation and/or extrapolation possible
        #fill in last observation if it is a duplicate
        if(Duplicated[End[n.strings]]==n_photos){
          Which.not = c(1,which(I_dup==0)+1)
          Which.not = Which.not[c(length(Which.not)-1,length(Which.not))]
          Tmp_coords[n_photos,1] = Tmp_coords[Which.not[2],1]+(n_photos-Which.not[2])*(Tmp_coords[Which.not[2],1]-Tmp_coords[Which.not[1],1])/(Which.not[2]-Which.not[1])
          Tmp_coords[n_photos,2] = Tmp_coords[Which.not[2],2]+(n_photos-Which.not[2])*(Tmp_coords[Which.not[2],2]-Tmp_coords[Which.not[1],2])/(Which.not[2]-Which.not[1])
          
          if(End[n.strings]==Begin[n.strings]){
            End=End[1:(n.strings-1)]
            Begin=Begin[1:(n.strings-1)]
            n.strings=n.strings-1
          }
          else End[n.strings]=End[n.strings]-1
        }
        
        #fill in internal missing records via interpolation if there are any
        if(n.strings>0){
          for(istr in 1:n.strings){
            begin.pl = Duplicated[Begin[istr]]-1
            end.pl = Duplicated[End[istr]]+1
            Last.coords = Tmp_coords[begin.pl,]
            Next.coords = Tmp_coords[end.pl,]
            slope.incr = end.pl - begin.pl
            Tmp_coords[Duplicated[Begin[istr]]:Duplicated[End[istr]],1]=Last.coords[1]+(Next.coords[1]-Last.coords[1])/slope.incr * c(1:(end.pl-begin.pl-1))
            Tmp_coords[Duplicated[Begin[istr]]:Duplicated[End[istr]],2]=Last.coords[2]+(Next.coords[2]-Last.coords[2])/slope.incr * c(1:(end.pl-begin.pl-1))
          }
        }
      }
    }
    Interp_loc[Which_photos,]=Tmp_coords
  }
  if(n_photos>=10){
    Tmp_coords = Coords[Which_photos,]
    I_dup = (Tmp_coords[2:n_photos,1]==Tmp_coords[1:(n_photos-1),1] & Tmp_coords[2:n_photos,2]==Tmp_coords[1:(n_photos-1),2])
    Duplicated = as.numeric(which(I_dup==1))+1
    Cur_dat = Data_ordered[Which_photos,]
    #Cur_df=data.frame(img_dt=as.numeric(difftime(Data_ordered[Which_photos,"correct_dt"],Data_ordered[Which_photos[1],"correct_dt"],units="secs")))
    #Cur_df$x = Tmp_coords[,1]
    #Cur_df$y = Tmp_coords[,2]
    #if(length(Duplicated)>0)Cur_df=Cur_df[-Duplicated,]
    #Cur_dat = Cur_dat[-Duplicated,]
    Duplicated = which(duplicated(Cur_dat$correct_dt)==TRUE)
    if(length(Duplicated)>0)Cur_dat=Cur_dat[-Duplicated,]
    Cur_df = st_drop_geometry(Cur_dat)
    Inits = list(a = c(Tmp_coords[1,1],0,Tmp_coords[1,2],0),P=diag(c(10000,10000,10000,10000)))
    #now run through crawl
    Fix=c(log(100),NA,2)
    #Fix=c(200,200,NA,NA)
    #if(iid %in% c(2133,2176,2196,2364,2573,2595,2876,5160,5898,6180,6241,6329,6634,6651,6678,6878,6903,7054,7179))Fix=c(log(150),log(150),NA,NA)
    my_fit <- crwMLE(mov.model=~1,err.model=list(x=~1),data=Cur_dat,Time.name="correct_dt",need.hess=FALSE,initial.state=Inits,fixPar=Fix,theta=c(0))
    
    if(is.character(my_fit)==TRUE){
      Fix=c(log(150),NA,2)  
      my_fit <- crwMLE(mov.model=~1,err.model=list(x=~1),data=Cur_dat,Time.name="correct_dt",need.hess=FALSE,initial.state=Inits,fixPar=Fix,theta=c(0))
    }
    elapsed = as.numeric(difftime( Cur_df[nrow(Cur_df),"correct_dt"],Cur_df[1,"correct_dt"],units="secs"))
    t_diff=elapsed/(n_photos-1)
    predTime=as.numeric(Cur_dat[1,"correct_dt"])[1]+c(0:(n_photos-1))*t_diff
    Preds=crwPredict(my_fit,predTime=predTime)
    Which_pred = which(Preds$locType=="p")
    Interp_loc[Which_photos,]=cbind(Preds$mu.x[Which_pred],Preds$mu.y[Which_pred])
  }
}
save(Interp_loc,file="Interp_loc.RDa")

#formulate new sp points object with interpolated locations
chess_geo=SpatialPointsDataFrame(coords=Interp_loc,data=Data_df,proj4str=CRS(st_crs(Data_ordered)$proj4string))
#fill in missing altitudes
Which.missing = which(is.na(chess_geo$gga_alt))  
for(i in 1:length(Which.missing))chess_geo$gga_alt[Which.missing[i]]=chess_geo$gga_alt[Which.missing[i]-1]


#formulate lookup table that gives camera angles, offset by aircraft, lens, and camera position
# easy this time around since we used a 100mm Zeiss lens for the machine vision cameras in every flight
Angles = expand.grid(camera=c("S","P","C"),lens=c(100))
Angles$vert = Angles$horiz = Angles$offset = rep(0,nrow(Angles)) 
Angles["vert"]=14
Angles["horiz"]=21
Angles["offset"]=25.5
Angles[which(Angles$camera=="C"),"offset"]=0  
Angles[,c("vert","horiz","offset")]=Angles[,c("vert","horiz","offset")]/360*2*pi  #covert to radians

get_footprint_corners_port <- function(Fl_angles,roll,bearing,Tmp_coords,alt){  
  #note: only works when abs(roll)< (horiz/2-offset)
  vert_width_far = alt/cos(Fl_angles[,"offset"]+0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  vert_width_near = alt/cos(Fl_angles[,"offset"]-0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  #horiz_width = alt*tan(0.5*Fl_angles[,"horiz"])
  horiz_dist_far = alt*tan(0.5*Fl_angles[,"horiz"]+Fl_angles[,"offset"]+roll)
  horiz_dist_near = alt*tan(Fl_angles[,"offset"]+roll-0.5*Fl_angles[,"offset"])
  #locate a point orthogonal to aircraft bearing at the outer edge of photograph that is closest (tangent) to the plane
  if(bearing<=(0.5*pi)){
    Tmp = c(-cos(0.5*pi-bearing),sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<=pi){
    Tmp = c(-cos(bearing-0.5*pi),-sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<=(1.5*pi)){
    Tmp = c(cos(1.5*pi-bearing),-sin(1.5*pi-bearing))
  }  
  if(bearing>(1.5*pi)){
    Tmp = c(cos(bearing-1.5*pi),sin(bearing-1.5*pi)) 
  }
  C_far = horiz_dist_far * Tmp + Tmp_coords
  C_near = horiz_dist_near * Tmp + Tmp_coords
  Poly_points = matrix(0,5,2)
  m=tan(bearing)
  tmp = sqrt(1+m^2)
  Tmp = rep(vert_width_far/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[1,] = C_far - Tmp
  Poly_points[2,] = C_far + Tmp
  Tmp = rep(vert_width_near/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[3,] = C_near + Tmp
  Poly_points[4,] = C_near - Tmp
  Poly_points[5,]=Poly_points[1,]
  Poly_points
}

get_footprint_corners_starboard <- function(Fl_angles,roll,bearing,Tmp_coords,alt){  
  #note: only works when abs(roll)< (horiz/2-offset)
  vert_width_far = alt/cos(Fl_angles[,"offset"]+0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  vert_width_near = alt/cos(Fl_angles[,"offset"]-0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  horiz_dist_far = alt*tan(0.5*Fl_angles[,"horiz"]+Fl_angles[,"offset"]-roll)
  horiz_dist_near = alt*tan(Fl_angles[,"offset"]-roll-0.5*Fl_angles[,"offset"])
  #locate a point orthogonal to aircraft bearing at the outer edge of photograph that is closest (tangent) to the plane
  if(bearing<=(0.5*pi)){
    Tmp = c(cos(0.5*pi-bearing),-sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<=pi){
    Tmp = c(cos(bearing-0.5*pi),sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<=(1.5*pi)){
    Tmp = c(-cos(1.5*pi - bearing),sin(1.5*pi - bearing))
  }  
  if(bearing>(1.5*pi)){
    Tmp = c(-cos(bearing-1.5*pi),-sin(bearing-1.5*pi)) 
  }
  C_far = horiz_dist_far * Tmp + Tmp_coords
  C_near = horiz_dist_near * Tmp + Tmp_coords
  Poly_points = matrix(0,5,2)
  m=tan(bearing)
  tmp = sqrt(1+m^2)
  Tmp = rep(vert_width_far/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[1,] = C_far - Tmp
  Poly_points[2,] = C_far + Tmp
  Tmp = rep(vert_width_near/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[3,] = C_near + Tmp
  Poly_points[4,] = C_near - Tmp
  Poly_points[5,]=Poly_points[1,]
  Poly_points
}

get_footprint_corners_center <- function(Fl_angles,roll,bearing,Tmp_coords,alt){  
  #note: only works when abs(roll)< (horiz/2-offset)
  vert_width_port = alt/cos(0.5*Fl_angles[,"horiz"]+roll)*sin(0.5*Fl_angles[,"vert"])
  vert_width_starboard = alt/cos(0.5*Fl_angles[,"horiz"]-roll)*sin(0.5*Fl_angles[,"vert"])
  horiz_dist_port = alt*tan(0.5*Fl_angles[,"horiz"]+roll)
  horiz_dist_starboard = alt*tan(0.5*Fl_angles[,"horiz"]-roll)
  #locate a point orthogonal to aircraft bearing at the outer edge of photograph that is closest (tangent) to the plane
  if(bearing<=(0.5*pi)){
    Tmp = c(-cos(0.5*pi-bearing),sin(0.5*pi-bearing)) 
  }
  if(bearing>(0.5*pi) & bearing<=pi){
    Tmp = c(-cos(bearing-0.5*pi),-sin(bearing-0.5*pi)) 
  }
  if(bearing>pi & bearing<=(1.5*pi)){
    Tmp = c(cos(1.5*pi-bearing),-sin(1.5*pi-bearing))
  }  
  if(bearing>(1.5*pi)){
    Tmp = c(cos(bearing-1.5*pi),sin(bearing-1.5*pi)) 
  }
  C_port = horiz_dist_port * Tmp + Tmp_coords
  C_starboard = -horiz_dist_starboard * Tmp + Tmp_coords
  Poly_points = matrix(0,5,2)
  m=tan(bearing)
  tmp = sqrt(1+m^2)
  Tmp = rep(vert_width_port/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[1,] = C_port - Tmp
  Poly_points[2,] = C_port + Tmp
  Tmp = rep(vert_width_starboard/tmp,2)
  Tmp[2]=Tmp[2]*m
  Poly_points[3,] = C_starboard + Tmp
  Poly_points[4,] = C_starboard - Tmp
  Poly_points[5,]=Poly_points[1,]
  Poly_points
}

save.image('effort_workspace.Rda')
#plot(Tmp_coords[1],Tmp_coords[2],xlim=c(Tmp_coords[1]-500,Tmp_coords[1]+500),ylim=c(Tmp_coords[2]-500,Tmp_coords[2]+500))
#points(C_far[1],C_far[2],col='orange')
#points(C_near[1],C_near[2],col='orange')
#points(C_port[1],C_port[2],col='cyan')
#points(C_starboard[1],C_starboard[2],col='cyan')
#points(Poly_points[1,1],Poly_points[1,2],col="red")
#points(Poly_points[2,1],Poly_points[2,2],col="purple")
#points(Poly_points[3,1],Poly_points[3,2],col="blue")
#points(Poly_points[4,1],Poly_points[4,2],col="green")

#  Formulate a list of times within flights for which banking turns are happening (for filter)
# base banking inference on camera for unique flight & cell combo that has the most photographs
# step 1: calculate change in bearing
FlightCell_unique_IDs = unique(FlightCell_ID) 
Photos = matrix(0,length(FlightCell_unique_IDs),3)
Delta_bearing=0
Delta_flight="dummy"
Delta_dt=chess_geo[1,]$correct_dt
n_cameras=3

#par(mfrow=c(3,3))
#for(iid in 78:86){
#  Which_photos = which(FlightCell_ID==FlightCell_unique_IDs[iid])
#  plot(chess_geo[Which_photos,])
#}

for(iid in 1:length(FlightCell_unique_IDs)){
  Which_photos = which(FlightCell_ID==FlightCell_unique_IDs[iid])
  n_photos=length(Which_photos)
  N_photos=rep(0,3)
  N_photos[1]=length(which(chess_geo[Which_photos,]$camera_loc=="P"))
  N_photos[2]=length(which(chess_geo[Which_photos,]$camera_loc=="S"))
  N_photos[3]=length(which(chess_geo[Which_photos,]$camera_loc=="C"))
  Photos[iid,]=N_photos
  max_n = max(N_photos)
  if(Photos[iid,3]==max_n)cam_name="C"  
  if(Photos[iid,1]==max_n)cam_name="P"
  if(Photos[iid,2]==max_n)cam_name="S"
  if(max_n>3){
    Which_cam = which(chess_geo[Which_photos,]$camera_loc==cam_name)
    Coords = Interp_loc[Which_photos[Which_cam],]
    Bearing=rep(0,max_n)
    Bearing = as.numeric(atan2(Coords[2:max_n,2]-Coords[1:(max_n-1),2],Coords[2:max_n,1]-Coords[1:(max_n-1),1]))
    Bearing[max_n]=Bearing[max_n-1]
    Bearing[2:(max_n-1)] = 0.5*(Bearing[1:(max_n-2)]+Bearing[2:(max_n-1)])
    #if(sum(Bearing<0)>0){  #make sure bearing is in [0,2*pi]
    #  Bearing[which(Bearing<0)]=Bearing[which(Bearing<0)]+2*pi
    #}
    Delta_dt = c(Delta_dt,chess_geo[Which_photos[Which_cam],]$correct_dt[1:(max_n-1)])
    Delta_flight=c(Delta_flight,rep(chess_geo[Which_photos[1],]$flight_num,max_n-1))
    Delta_bearing=c(Delta_bearing,Bearing[2:max_n]-Bearing[1:(max_n-1)])
  }
}
Delta_bearing=Delta_bearing[-1]
Delta_flight=Delta_flight[-1]
Delta_dt = Delta_dt[-1]

#filter based on 10 values with abs(Delta_bearing)
#first make a pass through to document sequence
I_filter=Seq=rep(0,length(Delta_bearing))
Abs_bearing=abs(Delta_bearing)
for(i in 2:length(Delta_bearing)){
  if(Abs_bearing[i]>0.01){
    Seq[i]=1
    if(Seq[i-1]>0 & abs(as.numeric(difftime(Delta_dt[i],Delta_dt[i-1],units="secs")))<5)Seq[i]=Seq[i-1]+1
  }
}
Which.eq.10 = which(Seq==10)  #375 turns 
n_turns=length(Which.eq.10)
Filter_df=data.frame(flight=rep("a",n_turns),start=rep(Delta_dt[1],n_turns),end=rep(Delta_dt[1],n_turns),stringsAsFactors=FALSE)
for(ipeak in 1:n_turns){
  cur.pl=Which.eq.10[ipeak]
  I_filter[(cur.pl-9):cur.pl]=1
  Filter_df[ipeak,"flight"]=Delta_flight[cur.pl]
  Filter_df[ipeak,"start"]=Delta_dt[cur.pl-9]
  Filter_df[ipeak,"end"]=Delta_dt[cur.pl]
  if(Seq[cur.pl+1]==11){
    flag=0
    cur.pl2=1
    while(flag==0){
      if(Seq[cur.pl+cur.pl2]>Seq[cur.pl+cur.pl2-1])I_filter[cur.pl+cur.pl2]=1
      else flag=1
      cur.pl2=cur.pl2+1
    }
    Filter_df[ipeak,"end"]=Delta_dt[cur.pl+cur.pl2-2]
  }
}
save(Filter_df,file="CHESS_turn_filter.Rda")


#filter out all photos that occur within turns
I_filter=rep(0,length(chess_geo))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(chess_geo$flight_num==Filter_df[irow,"flight"] & chess_geo$correct_dt>=Filter_df[irow,"start"] & chess_geo$correct_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}
chess_geo=chess_geo[which(I_filter==0),]  #remove 24378/826493 (3% of records)

#filter tracks that occur within turns
I_filter=rep(0,nrow(tracks))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(tracks$correct_dt>=Filter_df[irow,"start"] & tracks$correct_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}
tracks=tracks[which(I_filter==0),]  #remove 422/16690 (2.5% of records)

#filter bears that occur within turns
I_filter=rep(0,nrow(bears))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(bears$correct_dt>=Filter_df[irow,"start"] & bears$correct_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}
#no bears on turns!

save(chess_geo,file='chess_geo_filtered.Rda')
save(bears,file='bears_filtered.Rda')
save(bears_off,file='bears_off_filtered.Rda')
save(tracks,file="tracks_filtered.Rda")

#####  Now, produce footprint given flight track
Alt = .3048*chess_geo$gga_alt  #convert to meters
attr(chess_geo$correct_dt, "tzone") <- "UTC"  #convert time to UTC
Sides=c("P","S","C")
FlightCell_ID = paste(chess_geo$flight_num,chess_geo$GridID)
FlightCell_unique_IDs=unique(FlightCell_ID)
Photos = matrix(0,length(FlightCell_unique_IDs),3)
Photo_area = Flight_ID=Cell_ID=Day=rep(0,length(FlightCell_unique_IDs))
jday = floor(julian(chess_geo[["correct_dt"]], origin = as.POSIXct("2016-01-01", tz = "US/Alaska")))
chess_geo$day = as.numeric(jday - min(jday))+1

for(iid in 1:length(FlightCell_unique_IDs)){
  Which_photos = which(FlightCell_ID==FlightCell_unique_IDs[iid])
  n_photos=length(Which_photos)
  N_photos=rep(0,3)
  N_photos[1]=length(which(chess_geo[Which_photos,]$camera_loc=="P"))
  N_photos[2]=length(which(chess_geo[Which_photos,]$camera_loc=="S"))
  N_photos[3]=length(which(chess_geo[Which_photos,]$camera_loc=="C"))
  Photos[iid,]=N_photos
  flightid=chess_geo[Which_photos[1],]$flight_num
  Flight_ID[iid]=flightid
  Day[iid] = chess_geo[Which_photos[1],]$day
  Cell_ID[iid]=chess_geo[Which_photos[1],]$GridID
  Polys=vector("list",n_photos)
  cur_pl = 0
  for(iside in 1:3){
    side=Sides[iside]
    if(iside>1)cur_pl=cur_pl+N_photos[iside-1]
    Which_cam = which(chess_geo[Which_photos,]$camera_loc==Sides[iside])
    if(N_photos[iside]>0){
      Coords = chess_geo[Which_photos[Which_cam],]@coords
      Bearing=rep(0,N_photos[iside])
      if(N_photos[iside]>1){
        Bearing = as.numeric(atan2(Coords[2:N_photos[iside],2]-Coords[1:(N_photos[iside]-1),2],Coords[2:N_photos[iside],1]-Coords[1:(N_photos[iside]-1),1]))
        Bearing[N_photos[iside]]=Bearing[N_photos[iside]-1]
        Bearing[2:(N_photos[iside]-1)] = 0.5*(Bearing[1:(N_photos[iside]-2)]+Bearing[2:(N_photos[iside]-1)])
      }
      if(sum(Bearing<0)>0){  #make sure bearing is in [0,2*pi]
        Bearing[which(Bearing<0)]=Bearing[which(Bearing<0)]+2*pi
      }
      if(side=="P")footprint_fn = get_footprint_corners_port
      if(side=="C")footprint_fn = get_footprint_corners_center
      if(side=="S")footprint_fn = get_footprint_corners_starboard
      Tmp_angles = as.matrix(Angles[which(Angles[,"camera"]==side),c("vert","horiz","offset")])
      #now determine footprints
      for(iphoto in 1:N_photos[iside]){
        Polys[[cur_pl+iphoto]]=Polygons(list(Polygon(footprint_fn(Fl_angles=Tmp_angles,roll=0,bearing=Bearing[iphoto],Tmp_coords=Coords[iphoto,],alt=Alt[Which_photos[Which_cam[iphoto]]]))),paste(cur_pl+iphoto))
      }
    }
  }
  SPDF = SpatialPolygons(Polys,proj4string=CRS(st_crs(Data_ordered)$proj4string))  
  if(iid==10)save(SPDF,file="crap_spdf.Rda")
  Union = gUnionCascaded(SPDF)
  Photo_area[iid]=gArea(Union)
}  

Area_table = data.frame(flightid=Flight_ID,Grid_ID=Cell_ID,Area_m2 = Photo_area,Day=Day)
save(Area_table,file="CHESS_area_photographed.Rda")

