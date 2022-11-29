#######  
#
#    R script to calculate survey effort for CHESS (taking out turns, etc.) - reworked for 'large' 
#    (U.S. + Russia) Chukchi Grid
#
#    Last updated work with updated haul-out analysis by J. London
#    Paul Conn, 10 Mar 2021
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
effort <- sf::st_read(con, 
                         query = "SELECT DISTINCT flight_num, camera_loc, image_name, correct_dt, effort_type, gga_alt, geom 
                         FROM surv_chess.tbl_effort_raw WHERE effort_type = \'ON\'", 
                         geometry_column = "geom")
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
                            WHERE ((effort_type = \'BEAR\' OR effort_type = \ 'BEAR_KILL_SITE\')
                            OR (detection_skeyes = \'N\' AND effort_type = \'ON\'))
                            AND hotspot_type <> \'Duplicate\'",
                            geometry_column = "geom")
IDs = c("Bear01","Bear02","Bear20","Bear21","Bear22")  #bears detected visually or post hoc that didn't show up in 'on effort' infrared detection
#note this takes out cubs (their mother was detected on effort) and bears at point Barrow that weren't part of normal survey effort (polar_bear_id = 7-18 were at bone site)
bears_off = bears_off[bears_off$polar_bear_id %in% IDs,]  #take out cubs (their mother was detected on effort) and bears at the point (polar_bear_id = 7-18 were at bone site)
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

load("Chess_grid_list_all_Feb2021.RDa")

effort_df <- st_drop_geometry(effort)
Ecoords <- st_coordinates(effort)
effort_df = cbind(effort_df,Ecoords)
tracks_df <- st_drop_geometry(tracks)

grid_large <- Grid_list[[1]]

###### 2) intersect effort big chukchi grid

#remove images on land or outside grid  - this step takes a long time!!!!
AK = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ak_dcw.shp')
AK = st_transform(AK,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
Russia = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/russia_dcw.shp')
Russia = st_transform(Russia,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
Land=st_union(st_union(AK),st_union(Russia))
Grid_union = st_union(grid_large)
Intersect.land = st_intersects(effort,Land)
Intersect.grid = st_intersects(effort,Grid_union)

IL_vec = unlist(lapply(Intersect.land,"length"))
Which_IL = which(IL_vec>0)

IG_vec = unlist(lapply(Intersect.grid,"length"))
Which_OutGrid = which(IG_vec==0)

Which_out = unique(c(Which_IL,Which_OutGrid))
plot(effort[Which_out,],add=TRUE)  #plot which photographs will be removed. #7975 photos 

Effort_clean = effort[-Which_out,]
save(Effort_clean,file="effort_clean_US.RDa")

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
for(i in 1:827){
  Tmp_range = c(((i-1)*1000+1):(i*1000))
  Tmp_mat = unlist(st_intersects(Effort_clean[Tmp_range,],grid_large))
  Which_cell[Tmp_range]=Tmp_mat
}
Tmp_range=c(827001:827418)
Tmp_mat = unlist(st_intersects(Effort_clean[Tmp_range,],grid_large))
Which_cell[Tmp_range]=Tmp_mat
Effort_clean[["GridID"]]=Which_cell

save(Effort_clean,file='Effort_clean_large.Rda')

# attach gridID to tracks, bears
Which_cell = unlist(st_intersects(tracks,grid_large))
tracks[["GridID"]]=Which_cell
Which_cell = unlist(st_intersects(bears,grid_large))
bears[["GridID"]]=Which_cell
Which_cell = unlist(st_intersects(bears_off,grid_large))
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
    Cur_dat$tn = as.numeric(Cur_dat$correct_dt)
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
    my_fit <- crwMLE(mov.model=~1,err.model=list(x=~1),data=Cur_dat,Time.name="tn",need.hess=FALSE,initial.state=Inits,fixPar=Fix,theta=c(0))
    
    if(is.character(my_fit)==TRUE){
      Fix=c(log(150),NA,2)  
      my_fit <- crwMLE(mov.model=~1,err.model=list(x=~1),data=Cur_dat,Time.name="tn",need.hess=FALSE,initial.state=Inits,fixPar=Fix,theta=c(0))
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

save.image('effort_workspace_large.Rda')
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
save(Filter_df,file="CHESS_turn_filter_large.Rda")


#filter out all photos that occur within turns
I_filter=rep(0,length(chess_geo))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(chess_geo$flight_num==Filter_df[irow,"flight"] & chess_geo$correct_dt>=Filter_df[irow,"start"] & chess_geo$correct_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}
chess_geo=chess_geo[which(I_filter==0),]  #remove 24347/827418 (3% of records)

#filter tracks that occur within turns
I_filter=rep(0,nrow(tracks))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(tracks$correct_dt>=Filter_df[irow,"start"] & tracks$correct_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}
tracks=tracks[which(I_filter==0),]  #remove 854/23972 (3.5% of records)

#filter bears that occur within turns
I_filter=rep(0,nrow(bears))
for(irow in 1:nrow(Filter_df)){
  Which_photos = which(bears$correct_dt>=Filter_df[irow,"start"] & bears$correct_dt<=Filter_df[irow,"end"])
  I_filter[Which_photos]=1
}
#no bears on turns!

save(chess_geo,file='chess_geo_filtered_large.Rda')
save(bears,file='bears_filtered_large.Rda')
save(bears_off,file='bears_off_filtered_large.Rda')
save(tracks,file="tracks_filtered_large.Rda")

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
save(Area_table,file="CHESS_area_photographed_large.Rda")
save.image('CHESS_effort_workspace_pre_seals.RData')

#get seal count data
hotspots <- sf::st_read(con,
                        query = "SELECT *
                        FROM surv_chess.tbl_process",
                        geometry_column = "geom")
Animals <- hotspots[hotspots$hotspot_type=="Animal",]  
tabulate(as.factor(Animals$species_id))  #1137 bearded, 4732 ringed, 16 unknown seals - 0.3% - just ignore?
Seals = Animals[Animals$species_id %in% c("Ringed Seal","Bearded Seal","UNK Seal"),]
options(digits.secs=6)
DT = as.POSIXct(substr(Seals$timestamp,1,18),format="%Y%m%d%H%M%OS",tz="UTC")

#remove seals that aren't encountered while "on effort"
#reformat photograph names of "on effort" photos
change_photo_name = function(x){  
  Strings = unlist(strsplit(x,"_"))
  paste0("CHESS_FL",substr(Strings[3],3,nchar(Strings[3])),"_",Strings[4],"_",substr(Strings[6],3,8),"_",substr(Strings[6],9,18),"_COLOR-8-BIT.JPG")
}
Effort_photos_names = sapply(chess_geo$image_name,"change_photo_name")
Seals_on_effort = Seals[which(Seals$color_image_name %in% Effort_photos_names),]  #5892 seal recs in db; 5793 in 'on effort' photos (see below); 5371 after land, grid, and turn filters 

#calculating seals before turn filter
#kaka = sapply(effort_df$image_name,"change_photo_name")
#kaka2 = Seals[which(Seals$color_image_name %in% kaka),]  


#find grid ID, day for each seal group
n_surveyed=nrow(Area_table)
Seals_on_effort=st_transform(Seals_on_effort,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
n_seals = nrow(Seals_on_effort)
Intersects = st_intersects(Seals_on_effort,grid_large,sparse=FALSE)
Grid_seal = Day_seal = rep(0,n_seals)
DT = as.POSIXct(substr(Seals_on_effort$timestamp,1,18),format="%Y%m%d%H%M%OS",tz="UTC")
Jday = floor(julian(DT, origin = as.POSIXct("2016-01-01", tz = "US/Alaska")))
Seals_on_effort$day = as.numeric(Jday - min(Jday))+1
C_i = matrix(0,n_surveyed,3)
colnames(C_i)=c("Bearded","Ringed","Unk")

Bd_ad = Seals_on_effort[which(Seals_on_effort$species_id=="Bearded Seal" & Seals_on_effort$age_class=="Non-Pup"),]
Bd_pup = Seals_on_effort[which(Seals_on_effort$species_id=="Bearded Seal" & Seals_on_effort$age_class=="Pup"),]
Plot_df = data.frame(Age=c(rep("Adult",nrow(Bd_ad)),rep("Pup",nrow(Bd_pup))),Day=c(Bd_ad$day,Bd_pup$day))
library(ggplot2)
ggplot(Plot_df)+geom_histogram(aes(x=Day,group=Age,fill=Age))
summary(Bd_ad$day)
summary(Bd_pup$day)

Rd_ad = Seals_on_effort[which(Seals_on_effort$species_id=="Ringed Seal" & Seals_on_effort$age_class!="Pup"),]
Rd_pup = Seals_on_effort[which(Seals_on_effort$species_id=="Ringed Seal" & Seals_on_effort$age_class!="Non-Pup"),]
Plot_df = data.frame(Age=c(rep("Adult",nrow(Rd_ad)),rep("Pup",nrow(Rd_pup))),Day=c(Rd_ad$day,Rd_pup$day))
library(ggplot2)
ggplot(Plot_df)+geom_histogram(aes(x=Day,group=Age,fill=Age))
summary(Rd_ad$day)
summary(Rd_pup$day)


for(iseal in 1:n_seals){
  Grid_seal[iseal]=which(Intersects[iseal,]>0)
  #find surveyed cell 
  cur_which = which(Area_table$Grid_ID==Grid_seal[iseal] & Area_table$Day==Seals_on_effort$day[iseal])
  cur_sp_col = (Seals_on_effort$species_id[iseal]=="Bearded Seal") + 
      2*(Seals_on_effort$species_id[iseal]=="Ringed Seal") +
      3*(Seals_on_effort$species_id[iseal]=="UNK Seal")
  C_i[cur_which,cur_sp_col]=C_i[cur_which,cur_sp_col]+ Seals_on_effort$number_of_seals[iseal]
}


#Produce Day, Hour in solar hour for haulout predictions
Longitude.grid = coordinates(spTransform(rgeos::gCentroid(as(grid_large,'Spatial'),byid=TRUE), CRS("+proj=longlat +datum=WGS84")))[,1]
Longitude.sampled = Longitude.grid[Area_table$Grid_ID]
Hour.sampled = rep(0,n_surveyed)
DT_sampled = SolarT_sampled = rep(chess_geo$correct_dt[1],n_surveyed)
for(is in 1:n_surveyed){
  Cur_photos = chess_geo[which(chess_geo$GridID==Area_table[is,"Grid_ID"] & chess_geo$day==Area_table[is,"Day"]),]
  DT_sampled[is] = mean(Cur_photos$correct_dt)
  SolarT = solaR::local2Solar(Cur_photos$correct_dt,lon=Longitude.sampled[is])
  SolarT_sampled[is]=mean(SolarT)
  Hour.sampled[is] = mean(as.numeric(strftime(SolarT, format="%H")))
}
DayHour=matrix(0,n_surveyed,2)
DayHour[,1] = Area_table$Day
DayHour[,2] = Hour.sampled

# Haul-out predictions
grid.wx <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2016-04-07 00:00:00' AND '2016-06-01 00:00:00'")

Covs = data.frame(matrix(0,n_surveyed,6))
Hour = DayHour
colnames(Covs)=c("day","hour","precip","temp","pressure","wind")
#determine closest grid cell # to each location surveyed
#Distances = st_distance(st_as_sf(Effort_sp),st_centroid(st_as_sf(grid_large)))
#Which_closest = rep(0,n)
#for(i in 1:n)Which_closest[i]=which(Distances[i,]==min(Distances[i,]))
#Cell_num = grid_large[Which_closest,]$cell
Cell_num = grid_large$cell[Area_table$Grid_ID]
for(i in 1:n_surveyed){
  Cur_dat = grid.wx[which(grid.wx$cell==Cell_num[i]),]
  DT_diff = abs(Cur_dat$fdatetime_range_start-DT_sampled[i])
  cur_row = which(DT_diff == min(DT_diff))[1]
  Covs[i,3:5]=Cur_dat[cur_row,c("rast_acpcp","rast_air2m","rast_prmsl")]
  Covs[i,"wind"] = sqrt(Cur_dat[cur_row,"rast_uwnd"]^2+Cur_dat[cur_row,"rast_vwnd"]^2)
}

#transform covariates to scale used in GLMPMs
Covs[,"temp"]=(Covs[,"temp"]-270)/27
Covs[,"pressure"]=(Covs[,"pressure"]-100000)/10000
Covs[,"wind"]=Covs[,"wind"]/10
#now convert day, hour into format used in GLMPMs
Covs$day= (lubridate::yday(SolarT_sampled)-120)/10
Covs$day2=Covs$day^2
Covs$day3=Covs$day^3
Covs$hour = factor(lubridate::hour(SolarT_sampled),levels=c(0:23))
Hour = lubridate::hour(SolarT_sampled)
Sin1 = sin(pi*Hour/12)
Cos1 = cos(pi*Hour/12)
Sin2 = sin(pi*Hour/6)
Cos2 = cos(pi*Hour/6)
Sin3 = sin(pi*Hour/4)
Cos3 = cos(pi*Hour/4)
AS.vec = c("ADULT.F","ADULT.M","SUB","YOY")
L_list = vector("list",4)  
Covs$temp2 = Covs$temp
Covs$Dry = rep(0,n_surveyed)  #needed for model.matrix
Covs$sin1 = Sin1
Covs$sin2 = Sin2
Covs$sin3 = Sin3
Covs$cos1 = Cos1
Covs$cos2 = Cos2
Covs$cos3 = Cos3
Covs$Northing = 0.75  #value 'typical' of tagged bearded seals in the Chukchi   

# ------------------------------------------------------------------------------
#                   Bearded SEALS
# ------------------------------------------------------------------------------

# OLD approach using old data, age_sex
# load("c:/users/paul.conn/git/haulout/test_bearded.Rdata")
# 
# npar=nrow(test.bearded$covb)
# FE = test.bearded$fixed.effects
# for(iage in 1:4){
#   #design matrix
#   Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
#   L_list[[iage]]=model.matrix(test.bearded$fixed.formula,data=Covs)
# }
# L = matrix(0,n_surveyed*4,npar)
# counter=1
# for(i in 1:n_surveyed){
#   for(iage in 1:4){
#     L[counter,]= L_list[[iage]][i,]
#     counter=counter+1
#   }
# }
# 
# Ell = L%*%FE$estimate[-which(FE$estimate==0)]
# Cell = L%*%test.bearded$covb%*%t(L)
# Pred = plogis(Ell)
# Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n_surveyed*4,n_surveyed*4)
# Pred_var = Deriv%*%crossprod(Cell,Deriv)
# 
# #stable stage distribution combo prediction
# Pi = matrix(c(0.27,0.23,0.38,0.12),1,4)
# M = kronecker(diag(n_surveyed),Pi)
# Pred_ho_bearded = M %*% Pred
# Var_ho_bearded = M %*% tcrossprod(Pred_var,M)
# 
# #output means, var-cov matrices for HO predictions
# HO_out = list(Mu_bd= Pred_ho_bearded,Var_bd = Var_ho_bearded)             
# save(HO_out,file="Haulout_bearded_chess_US.RData")

#load("c:/users/paul.conn/git/haulout/london_analysis/berchukhaulout/data/fit_bearded.Rdata")

remotes::install_github('jmlondon/berchukFits') #London et al. 2021 availability predictions
library(berchukFits)
data(bearded_fit)
fit_bearded=bearded_fit

npar=nrow(fit_bearded$covb)
FE = fit_bearded$fixed.effects

Covs$dry=Covs$Dry
Covs$northing=Covs$Northing
L = model.matrix(fit_bearded$fixed.formula,data=Covs)

Ell = L%*%FE$estimate
Cell = L%*%fit_bearded$covb%*%t(L)
Pred_ho_bearded = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n_surveyed,n_surveyed)
Var_ho_bearded = Deriv%*%crossprod(Cell,Deriv)

#output means, var-cov matrices for HO predictions
HO_out = list(Mu_bd= Pred_ho_bearded,Var_bd = Var_ho_bearded)             
save(HO_out,file="Haulout_bearded_chess_US_Mar2021.RData")


ep = 247/258
ep2 = ep^2
varp = ep*(1-ep)/258  #detection probability
EX2 = HO_out$Mu_bd^2
Tmp = diag(as.vector(HO_out$Mu_bd^2))
n.obs = length(HO_out$Mu_bd)
for(i in 1:(n.obs-1)){
  for(j in (i+1):n.obs){
    Tmp[i,j]=Tmp[j,i]=HO_out$Mu_bd[i]*HO_out$Mu_bd[j]
  }
}
Var_bd = Tmp*varp + (varp+ep2)*HO_out$Var_bd 

library(Matrix)
Thin_bd = c(HO_out$Mu_bd*ep)

save.image("effort_workspace_ringed_start.RDa")

#### ringed seal corrections
#ringed seals predictions using GAM model
load("c:/users/paul.conn/git/haulout/ringed/big_gam_ringed_ho_model.RData")  
Covs_ringed = data.frame(sea = factor(rep("Chukchi",n_surveyed),levels=c("Bering","Chukchi")),
                         jday= Area_table$Day+97,
                         dmelt_Markus=rep(0,n_surveyed),
                         spring_temp = rep(0,n_surveyed))
for(i in 1:n_surveyed){
  Covs_ringed$dmelt_Markus[i]=Grid_list[[Area_table[i,"Day"]]]$dmelt_Markus[Area_table[i,"Grid_ID"]]
  Covs_ringed$spring_temp[i]=Grid_list[[Area_table[i,"Day"]]]$spring_temp[Area_table[i,"Grid_ID"]] 
}
Pred_ho_ringed = mgcv::predict.gam(big_model,newdata=Covs_ringed,type="response",se.fit=TRUE)

#thin_rd = 0.65*23/29*ep  #includes disturbance, p
Thin_rd = 23/29*ep*Pred_ho_ringed[[1]]  #includes disturbance, p, haul-out
Sigma_thin = vector("list",2)
Sigma_thin[[1]]=Var_bd  
Thin = c(Thin_bd,Thin_rd)

#now ringed seal variance
eb=23/29
epb2=(eb*ep)^2
varb = eb*(1-eb)/29
#var_rd = ((0.2*0.65)^2+0.65^2)*(varb+eb^2)*(varp+ep^2)-0.65^2*eb^2*ep^2  #exact formula  - used when 0.65 was asymptote
varpb = varb*varp + varb*ep^2 + varp * eb^2
EX2 = Pred_ho_ringed$fit^2
Tmp = diag(as.vector(Pred_ho_ringed$fit^2))
for(i in 1:(n.obs-1)){
  for(j in (i+1):n.obs){
    Tmp[i,j]=Tmp[j,i]=Pred_ho_ringed$fit[i]*Pred_ho_ringed$fit[j]
  }
}
Var_rd = Tmp*varpb + (varpb+epb2)*diag(Pred_ho_ringed$se.fit^2)   #assuming prediction variance independent


Sigma_thin[[2]]=Var_rd
Sigma_thin = as(bdiag(Sigma_thin),"dgTMatrix")
#compute logit scale distribution using delta method
diff_logit <- function(x) -1/(x*(x-1))
Thin_logit = log(Thin/(1-Thin))
Diff_logit = diag(diff_logit(Thin))
Sigma_logit_thin = Diff_logit %*% Sigma_thin %*% t(Diff_logit)
diag(Sigma_logit_thin)=diag(Sigma_logit_thin)+0.000001  #adding a  small nugget to  prevent numerical problems







#for extra day effect on availability...
X_day = matrix(0,n_surveyed,2)
X_day[,1] = DayHour[,1]
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2] = (X_day[,1])^2 
X_day = kronecker(diag(2),X_day)

n_cells = nrow(grid_large)
S_i = (Area_table[,"Day"]-1)*n_cells+Area_table[,"Grid_ID"] #locations and times surveyed
Area_trans = as.numeric(Area_table$Area_m2 / st_area(grid_large[1,]))
Data_chess_us = list("C_i"=C_i, "P_i"=Area_trans,"S_i"=S_i-1,"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=as(Sigma_logit_thin,'dgTMatrix'),"X_day"=X_day)
Data_chess_us$h_mean = c(mean(Thin[1:n_surveyed]),mean(Thin[(n_surveyed+1):(2*n_surveyed)]))

save(Data_chess_us,file="CHESS_us_data_TMB_Mar2021.RData")
