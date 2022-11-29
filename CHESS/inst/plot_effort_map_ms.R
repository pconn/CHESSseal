#plot_survey_area.R


library(sp)
library(rgdal)
library(nPacMaps)  #from Josh London
library(rgeos)
library(doBy)
library(ggplot2)
library(plyr)
library(RPostgreSQL)
library(sf)

load('Chess_grid_list_all.RDa')
Cur.grid = Grid_list[[1]]
load('effort_clean_US.RDa')

Russian_tracks = st_read("c:/users/paul.conn/git/chess/RussianFlightTracks/transects_by_date.shp")
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Russian_tracks = st_transform(Russian_tracks,laea_180_proj)

#remove segments on land or outside grid  
AK = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/ak_dcw.shp')
AK = st_transform(AK,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
Russia = read_sf('c:/users/paul.conn/git/BeaufortPower/shapefiles/russia_dcw.shp')
Russia = st_transform(Russia,"+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")


eez = st_read("c:/users/paul.conn/git/chess/RussianFlightTracks/AK_EEZ.shp")

#ggplot()+geom_point(data=Points.df,aes(x=GPSLONG_INT,y=GPSLAT_INT))

shelf <- SpatialLinesDataFrame(gBoundary(Shelf_break[1,]),data.frame(id=1))
eez <- SpatialLinesDataFrame(gBoundary(EEZ_Alaska),data.frame(id=1))

# fortify shelf and eez
shelf <- fortify(shelf)
eez <- fortify(eez)

# use some colorbrewer colors
red <- "#e41a1c"
blue <- "#377eb8"
orange <- "#ff7f00"
brown <- "#a65628"

Effort_subsample = Effort_clean[c(1:16548)*50,]
# make our plot
p <- ggplot() + geom_sf(data=boss_grid,fill="white",color="gray") + 
   geom_sf(data=AK[,1],fill="grey60",color="grey60",size=1.2) +
   geom_sf(data=Russia[,1],fill="grey60",color="grey60",size=1.2) +
   geom_sf(data=Russian_tracks,color="red",size=1)+
   geom_sf(data=Effort_subsample,color="blue",size=0.3)+
   coord_sf() +
   scale_y_continuous(breaks = c(63,65,67,69,71,73,77),limits=c(-500000,900000))+
   scale_x_continuous(breaks = c(176,178,180,2,4,6,8),limits=c(-2900000,-1500000))
  
  
   #geom_sf(data=eez,color="brown",size=0.3)

p <- p + labs(x="Easting",y="Northing")
p <- p + theme(text=element_text(size=16),title=element_text(size=16),axis.ticks = element_blank(),axis.text=element_blank())


#p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) #remove slash in legend
p

# save a pdf
#ggsave(file="BOSS_survey_ribbon.tiff",dpi=300)

pdf(file="BOSS_effort_map.pdf")
p
dev.off()



