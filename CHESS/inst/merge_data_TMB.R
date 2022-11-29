#merge US, Russian data together into one TMB input.  
library(INLA)

load("CHESS_russia_Data_TMB_Apr2021.RData")
load("CHESS_us_data_TMB_Feb2021.RData")
load("Chess_grid_list_all_Feb2021.RDa")

#assemble melt data for dates and times surveyed (not done correctly in Data_rus)
Melt_stacked = rep(0,Data_rus$n_s*Data_rus$n_t)
for(it in 1:Data_rus$n_t)Melt_stacked[(it-1)*Data_rus$n_s+c(1:Data_rus$n_s)]=Grid_list[[it]]$dmelt_MDSDA
Melt_rus = Melt_stacked[Data_rus$S_i]
Melt_us = Melt_stacked[Data_chess_us$S_i]
Data_rus$X_s$ice2 = Data_rus$X_s$sea_ice^2
Data_rus$X_s$dist2 = Data_rus$X_s$dist_mainland^2

library(sf)
Dist_strait = log(as.numeric(st_distance(Grid_list[[1]],Grid_list[[1]][1351,]))/1000 + 25)
Data_rus$X_s$dist_strait = rep(Dist_strait,Data_rus$n_t)

n_cells = nrow(Grid_list[[1]])

CHESS_data = list(C_rus = as.matrix(Data_rus$C_i), C_us=as.matrix(Data_chess_us$C_i), Nophoto_i = Data_rus$Nophoto_i,
                  Prop_photo_i = Data_rus$Prop_photo_i, Melt_rus_i = Melt_rus, Melt_us_i = Melt_us, P_rus_i = Data_rus$P_i,
                  P_us_i = Data_chess_us$P_i, A_s = Data_rus$A_s, S_rus_i = Data_rus$S_i,
                  S_us_i = Data_chess_us$S_i, X_s = as.matrix(Data_rus$X_s), Thin_mu_rus_i = Data_rus$thin_mu_logit,
                  Thin_mu_us_i = Data_chess_us$thin_mu_logit, Thin_sigma_rus = Data_rus$Sigma_logit_thin,
                  Thin_sigma_us = Data_chess_us$Sigma_logit_thin, X_day_rus = Data_rus$X_day, X_day_us = Data_chess_us$X_day,
                  n_s = Data_rus$n_s, n_sp = as.integer(Data_rus$n_sp), n_t = Data_rus$n_t, h_mean_rus = Data_rus$h_mean[1],
                  h_mean_us = Data_chess_us$h_mean[1])
CHESS_data$X_day_rus = CHESS_data$X_day_rus[1:(nrow(CHESS_data$X_day_rus)/2),1:2]
CHESS_data$X_day_us = CHESS_data$X_day_us[1:759,1:2]  #don't need day adjustment for ringed seals since that's modeled separately 


#Chaunskaya Bay mask for bearded seals

CHESS_data$Which_cells_bearded = c(1:CHESS_data$n_s)
CHESS_data$Which_cells_bearded=CHESS_data$Which_cells_bearded[-c(955:957,995:997,1030:1032,1060:1062)]-1
CHESS_data$I_cell=rep(1,CHESS_data$n_s) 
CHESS_data$I_cell[c(955:957,995:997,1030:1032,1060:1062)]=0

save(CHESS_data,file="CHESS_seal_data_no_area_filter.RData")
           
#Remove data where < 0.1% of a grid cell was surveyed
Which = which(CHESS_data$P_us_i<0.001)
CHESS_data$C_us = CHESS_data$C_us[-Which,]
CHESS_data$Melt_us_i = CHESS_data$Melt_us[-Which]
CHESS_data$P_us_i = CHESS_data$P_us_i[-Which]
CHESS_data$S_us_i = CHESS_data$S_us_i[-Which]
Which2 = c(Which,Which+length(Melt_us))
CHESS_data$Thin_mu_us_i = CHESS_data$Thin_mu_us_i[-Which2]  #changed 2/25/21 - ringed seals now have >1 thinning par!
CHESS_data$Thin_sigma_us = CHESS_data$Thin_sigma_us[-Which2,-Which2]
CHESS_data$X_day_us = CHESS_data$X_day_us[-Which,]
Which = which(CHESS_data$P_rus_i<0.001)
CHESS_data$C_rus = CHESS_data$C_rus[-Which,]
CHESS_data$Melt_rus_i = CHESS_data$Melt_rus[-Which]
CHESS_data$P_rus_i = CHESS_data$P_rus_i[-Which]
CHESS_data$S_rus_i = CHESS_data$S_rus_i[-Which]
Which2 = c(Which,Which+length(Melt_rus))
CHESS_data$Thin_mu_rus_i = CHESS_data$Thin_mu_rus_i[-Which2]
CHESS_data$Thin_sigma_rus = CHESS_data$Thin_sigma_rus[-Which2,-Which2]
CHESS_data$X_day_rus = CHESS_data$X_day_rus[-Which,]
CHESS_data$Nophoto_i = CHESS_data$Nophoto_i[-Which]
CHESS_data$Prop_photo_i = CHESS_data$Prop_photo_i[-Which]

CHESS_data$Day_us_i = ceiling(CHESS_data$S_us_i / n_cells)
CHESS_data$Day_rus_i = ceiling(CHESS_data$S_rus_i / n_cells)


Loc_s = cbind(CHESS_data$X_s[,"easting"],CHESS_data$X_s[,"northing"]) 
Loc_s[,1] = Loc_s[,1]/mean(Loc_s[,1])
Loc_s[,2] = Loc_s[,2]/mean(Loc_s[,2])
mesh = inla.mesh.create( Loc_s )
CHESS_data$spde <- (inla.spde2.matern(mesh, alpha=1)$param.inla)[c("M0","M1","M2")]


save(CHESS_data,file="CHESS_data_Nov2021.RData")