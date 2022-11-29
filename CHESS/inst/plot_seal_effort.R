# plot seal effort for CHESS

#note this is currently broken since ptolemy, etc. had to be updated from nPacMaps
#and various deprecated other things....

load("CHESS_area_photographed.Rda")
load("tracks_filtered.Rda")
load('CHESS_data.RData')




laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Bengtson<-sf::st_read("c:/users/paul.conn/git/chukchipower/shapefiles/chukchi_strata.shp")
Bengtson<-sf::st_transform(Bengtson,laea_180_proj)
Bengtson_union<-sf::st_combine(Bengtson)

library(sf)
st_write(Bengtson_union,"./shapefiles/plot1/Bengtson.shp")


### plot observations maps
library(ptolemy) 

### plot Russian + U.S. data
Russian_tracks = st_read("c:/users/paul.conn/git/chess/RussianFlightTracks/transects_by_date.shp")
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
Rus_tracks = st_transform(Russian_tracks,laea_180_proj)
track_union = st_union(st_union(tracks),st_union(Rus_tracks))
#Rus_effort = read.csv('Russian_data_upd.csv')

st_write(tracks,"./shapefiles/plot1/US_tracks.shp")


load('CHESS_grid_list_all.RDa')
Grid_sf = st_as_sf(Grid_list[[1]])
b <- sf::st_bbox(track_union)
b[c(1,2)]=b[c(1,2)]-100000
b[c(3,4)]=b[c(3,4)]+100000

st_write(Grid_sf,"./shapefiles/plot1/analysis.grid.shp")



#Rus_all = st_as_sf(Rus_sightings, coords = c("Lon", "Lat"), 
#                crs = 4326, agr = "field")

library(ggplot2)
arctic_map = ptolemy::extract_gshhg(xlims=c(180-50,180+50), ylims=c(40,90), resolution = "i", epsg=3572, simplify=TRUE)
arctic_map = st_transform(arctic_map,laea_180_proj)

st_write(arctic_map, "./shapefiles/plot1/land.shp")

Bay_mask = c(955:957,995:997,1030:1032,1060:1062)
crap = ggplot() + geom_sf(data=Bengtson_union,fill="yellow",col="yellow")
crap = crap + geom_sf(data=Grid_list[[1]][-Bay_mask,],col="khaki3",alpha=0.2)
#crap = crap + geom_sf(data=Grid_list[[1]][Bay_mask,],fill="orange")
crap = crap + geom_sf(data=tracks,size=.2)
crap = crap + geom_sf(data=Rus_tracks,size=1)
crap = crap + geom_sf(data=arctic_map,fill="grey60",size=0.2)
crap = crap + coord_sf(xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))
crap
pdf("CHESS_seal_effort.pdf")
#pdf("Effort.pdf")
crap
dev.off()

jpeg("CHESS_seal_effort.jpg")
#pdf("Effort.pdf")
crap
dev.off()

