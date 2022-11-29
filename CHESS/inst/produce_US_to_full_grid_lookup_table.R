# for US to full CHESS grid conversion
load("Chess_grid_list_all.RDa")   ###RUS + US grid
Big_grid = Grid_list[[1]]

load("US_chess_grid_list.RDa")
Little_grid = Grid_list[[1]]

Little_centroids = sf::st_centroid(Little_grid)

Inter = st_intersects(Little_centroids,Big_grid,sparse=FALSE)
Cell_lookup_table = rep(0,nrow(Inter))
for(i in 1:nrow(Inter))Cell_lookup_table[i] = which(Inter[i,]==1)

save("Cell_lookup_table",file="Cell_lookup_table.RData")


# output identifier for mostly in U.S. waters
Big_grid$I_us = rep(0,nrow(Big_grid))
Big_grid$I_us[Cell_lookup_table]=1
plot(Big_grid[,"I_us"])  #looks like the 'old' US grid didn't fully go out to current northern boundary

#Big_grid$I_us[c(1,2,11,19,41,71,105,143,221,263,264,308)]=1
Big_grid$I_us[c(1,2,11,19,41,71,105,143,221,264)]=1

plot(Big_grid[,"I_us"])  #looks like the 'old' US grid didn't fully go out to current northern boundary


EEZ = sf::st_read("C:/Users/paul.conn/git/BeaufortPower/shapefiles/alaska_EEZ_line_edit.shp")
library(ggplot2)
Big_grid$I_us = as.character(Big_grid$I_us)
ggplot()+geom_sf(data=Big_grid[,"I_us"],aes(fill=I_us))+geom_sf(data=EEZ)

I_US = which(Big_grid$I_us==1)
save(I_US,file="I_US.RData")