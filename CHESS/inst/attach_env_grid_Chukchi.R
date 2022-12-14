library(sp)
library(grid)
library(rgdal)
library(rgeos)
library(raster)
library(nPacMaps)  #from Josh London
library(maptools)
library(gpclib)
library(automap)
library(ggplot2)
library(plyr)
#library(bass)

MEAN_ADJUST=TRUE  
#if(DEBUG)source("c:/users/paul.conn/git/bass/bass/R/bass.R")

#fun <- function(x) {
#	ifelse(x < 251, x/2.5, NA)
#}
#sic_raster <- calc(r, fun)

###some possible projections

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
		"+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

Polar_stereo<-paste("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

AK_albers_proj <-paste("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0"," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#just using an old sea ice raster to get the adjancency matrix for spatial modeling
r<-raster("//nmfs/akc-nmml/Polar//Data/Environ/SeaIce/Nimbus/Sea_Ice_Concentrations_from_Nimbus-7_SMMR_and_DMSP_SSM_I-SSMIS_Passive_Microwave/nt_20120401_f17_v01_n.bin.reproj.tif")

sic_raster <- projectRaster(r, res = c(25067.53, 25067.53),
		crs = laea_180_proj)

data(alaska_dcw)
data(russia_dcw)
#alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
#russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))

x_min <- -400000
x_max <- 1500000
y_min <- -2630000
y_max <- -1600000
#y_max <- -1400000

pep_ext <- extent(x_min, x_max, y_min, y_max)
sic_raster <- crop(sic_raster, pep_ext, snap = "near")
dim_raster=dim(sic_raster)
pep_ext <- extent(sic_raster)

#this will load current sic raster
#load('ChukchiBering_big_square_raster.Rda')


#x_max<-870000
#pep_ext <- extent(x_min, x_max, y_min, y_max)
#sic_raster <- crop(sic_raster, pep_ext, snap = "near")
#dim_raster=dim(sic_raster)
#pep_ext <- extent(sic_raster)


Grid<-readShapeSpatial("c:/users/paul.conn/git/chukchipower/shapefiles/ChukchiBeaufort_SSMI-EEZ_Grid",proj4string=CRS(Polar_stereo))

alaska_dcw <- spTransform(alaska_dcw, CRS(Polar_stereo))
russia_dcw <- spTransform(russia_dcw, CRS(Polar_stereo))


#plot(sic_raster, col = colorRampPalette(c("dark blue", "deepskyblue","skyblue", "lightskyblue", #"white"))(20), main = "SSM/I Sea-Ice Concentration and PEP BOSS Extent\n01 April 2012")
plot(Grid)
plot(alaska_dcw, col = "black", add = TRUE)
plot(russia_dcw, col = "black", add = TRUE)
#plot(pep_ext, col = "red", add = TRUE)



Land=list(alaska=alaska_dcw,russia=russia_dcw)

#the following takes awhile.  Instead, consider loading cur.Grid.Rdat
Grid_poly=add.prop.land(Grid=Grid,Land=Land)
#remove initial layer
save.image("cur.Chukchi.Grid.Rdat")
#load("cur.Bering.Grid.Rdat") 


#define separate SpPolyDFs for mainland
Area_alaska=gArea(alaska_dcw,byid=TRUE)
Area_russia=gArea(russia_dcw,byid=TRUE)
Alaska_mainland=alaska_dcw[which(Area_alaska==max(Area_alaska)),]
Russia_mainland=russia_dcw[which(Area_russia==max(Area_russia)),]

#Attach distance to mainland for each cell (distance from grid cell centroid to landscape polygon)
Grid_points=gCentroid(Grid_poly,byid=TRUE)
Dist_AK=gDistance(Grid_points,Alaska_mainland,byid=TRUE)
Dist_Rus=gDistance(Grid_points,Russia_mainland,byid=TRUE)
Dist_mainland=apply(cbind(as.vector(Dist_AK),as.vector(Dist_Rus)),1,'min')
if(MEAN_ADJUST)Dist_mainland=Dist_mainland/mean(Dist_mainland)
Grid_poly[["dist_mainland"]]=Dist_mainland

#Attach easting and northing coordinates
Centroids=gCentroid(Grid_poly,byid=TRUE)
Grid_poly[["Easting"]]=Centroids$y
Grid_poly[["Northing"]]=Centroids$x

if(MEAN_ADJUST){
  Grid_poly[["Easting"]]=Grid_poly[["Easting"]]/mean(Grid_poly[["Easting"]])
  Grid_poly[["Northing"]]=Grid_poly[["Northing"]]/mean(Grid_poly[["Northing"]])
}
Grid_poly@data=Grid_poly@data[,5:9]

#if(MEAN_ADJUST)Grid_poly$depth=Grid_poly$depth/mean(Grid_poly#$depth) #standradize by dividing by its mean
#

save(Grid_poly,file="USChukchi_Beaufort_Grid_static.Rdat")

#limit to chukchi on east side
Grid.Beaufort=Grid_poly[Grid_poly@data[,"Sea"]=="Beaufort",]
Grid.Chukchi=Grid_poly[Grid_poly@data[,"Sea"]=="Chukchi",]

Adj.Chukchi=gIntersects(Grid.Chukchi,byid=TRUE)*1.0
Adj.Beaufort=gIntersects(Grid.Beaufort,byid=TRUE)*1.0


Adj=rect_adj(x=dim_raster[2],y=dim_raster[1],byrow=TRUE)

Grid.chukchi=list(Grid=Grid.Chukchi,Adj=Adj.Chukchi)
save(Grid.chukchi,file="Grid_chukchi.Rda")
Grid.beaufort=list(Grid=Grid.Beaufort,Adj=Adj.Beaufort)
save(Grid.beaufort,file="Grid_beaufort.Rda")

