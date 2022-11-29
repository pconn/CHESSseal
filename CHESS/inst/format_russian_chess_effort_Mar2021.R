#### format russian CHESS effort, including weather for haulout predictions for bearded seals



library('RPostgreSQL')
library('sf')
library('rgeos')
library('sp')
library('Matrix')
set.seed(12345)  #just so pseudo-zeroes will always be in the same place

load("c:/users/paul.conn/git/haulout/london_analysis/berchukhaulout/data/fit_bearded.Rdata")

load("Chess_grid_list_all_Feb2021.RDa")
grid_large <- as(Grid_list[[1]],"Spatial")
n_cells = length(grid_large)

Effort = read.csv('Russian_data_upd.csv',header=TRUE)
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
spdf <- SpatialPoints(coords = Effort[,5:6],
                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
Effort_sp = spTransform(spdf,CRS(laea_180_proj))

#plot(grid_large[,"fp"])
#plot(Effort_sp,col='blue',add=T)
#remove effort that intersects w/ U.S. BOSS Grid (won't want to double count)
I.intersect = as.vector(gIntersects(Effort_sp,gUnionCascaded(grid_large),byid=TRUE))

Effort_sp = SpatialPointsDataFrame(Effort_sp,Effort[,1:11])
Effort_sp = Effort_sp[I.intersect,]

n_surveyed=nrow(Effort_sp)

#which cell surveyed
rownames(grid_large@data)=c(1:length(grid_large))  #rename IDs
Grid_sampled = rep(0,n_surveyed)
Intersects = gIntersects(Effort_sp,grid_large,byid=TRUE)
for(i in 1:n_surveyed)Grid_sampled[i]=which(Intersects[,i]==1)

#Produce Day, Hour in solar hour for haulout predictions
Time = as.character(Effort_sp$Time)

# for(it in 1:length(Time)){  #convert time to standard format
#   Tmp = strsplit(Time[it],":")
#   Last = strsplit(Tmp[[1]][3]," ")
#   Last_tmp = Last[[1]][1]
#   if(length(Last[[1]])>1){  #there's an 'AM' or 'PM' attached
#     if(Last[[1]][2]=='PM' & as.numeric(Tmp[[1]][1])<12){
#       Tmp[[1]][1]=as.character(as.numeric(Tmp[[1]][1])+12)
#     }
#   }
#   Time[it]=paste(Tmp[[1]][1],Tmp[[1]][2],Last_tmp,sep=":")
# }
Cell_IDs=as.numeric(rownames(grid_large@data))
Longitude.sampled = coordinates(spTransform(Effort_sp,CRS("+proj=longlat +datum=WGS84")))[,1]
ymdt = paste(Effort_sp$Date,Time)
Effort_sp$dt = as.POSIXct(ymdt,format="%m/%d/%Y %H:%M",tz="Asia/Magadan")
Effort_sp$SolarT = solaR::local2Solar(Effort_sp$dt,lon=Longitude.sampled)

#Day and Hour
date.start=as.Date("2016-04-07") #Alaska time
DayHour_rus=matrix(0,nrow(Effort_sp),2)
DayHour_rus[,1] = as.numeric(as.Date(Effort_sp$SolarT)-date.start+1)  
DayHour_rus[,2] = as.numeric(strftime(Effort_sp$SolarT, format="%H"))

con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))



grid.wx <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2016-04-17 00:00:00' AND '2016-05-19 00:00:00'")

n = nrow(Effort_sp)
Covs = data.frame(matrix(0,n,6))
Hour = DayHour_rus
colnames(Covs)=c("day","hour","precip","temp","pressure","wind")
#determine closest grid cell # to each location surveyed
Distances = st_distance(st_as_sf(Effort_sp),st_centroid(st_as_sf(grid_large)))
Which_closest = rep(0,n)
for(i in 1:n)Which_closest[i]=which(Distances[i,]==min(Distances[i,]))
Cell_num = grid_large[Which_closest,]$cell
for(i in 1:length(Effort_sp)){
  Cur_dat = grid.wx[which(grid.wx$cell==Cell_num[i]),]
  DT_diff = abs(Cur_dat$fdatetime_range_start-Effort_sp[i,]$dt)
  cur_row = which(DT_diff == min(DT_diff))[1]
  Covs[i,3:5]=Cur_dat[cur_row,c("rast_acpcp","rast_air2m","rast_prmsl")]
  Covs[i,"wind"] = sqrt(Cur_dat[cur_row,"rast_uwnd"]^2+Cur_dat[cur_row,"rast_vwnd"]^2)
}

#transform covariates to scale used in GLMPMs
Covs[,"temp"]=(Covs[,"temp"]-270)/27
Covs[,"pressure"]=(Covs[,"pressure"]-100000)/10000
Covs[,"wind"]=Covs[,"wind"]/10
#now convert day, hour into format used in GLMPMs
Covs$day= (lubridate::yday(Effort_sp$dt)-120)/10
Covs$day2=Covs$day^2
Covs$day3=Covs$day^3
Covs$hour = factor(lubridate::hour(Effort_sp$SolarT),levels=c(0:23))
Hour = lubridate::hour(Effort_sp$SolarT)
Sin1 = sin(pi*Hour/12)
Cos1 = cos(pi*Hour/12)
Sin2 = sin(pi*Hour/6)
Cos2 = cos(pi*Hour/6)
Sin3 = sin(pi*Hour/4)
Cos3 = cos(pi*Hour/4)
AS.vec = c("ADULT.F","ADULT.M","SUB","YOY")
L_list = vector("list",4)  
Covs$temp2 = Covs$temp
Covs$Dry = rep(0,n)  #needed for model.matrix
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
save(HO_out,file="Haulout_dists_russia_Mar2021.RData")


Sigma_thin = vector("list",2)
ep = 68/70
ep2 = ep^2
eb = 48/51
varb = eb*(1-eb)/51
varp = ep*(1-ep)/70  #detection probability
Thin_bd = HO_out$Mu_bd*ep*eb
epb = ep*eb
epb2 = epb^2
varpb =  ep^2*varb + eb^2*varp + varp*varb  #ep * eb

EX2 = HO_out$Mu_bd^2
Tmp = diag(as.vector(HO_out$Mu_bd^2))
n.obs = length(HO_out$Mu_bd)
for(i in 1:(n.obs-1)){
  for(j in (i+1):n.obs){
    Tmp[i,j]=Tmp[j,i]=HO_out$Mu_bd[i]*HO_out$Mu_bd[j]
  }
}
Var_bd = Tmp*varpb + (varpb+epb2)*HO_out$Var_bd 

Sigma_thin[[1]]=Var_bd


load("c:/users/paul.conn/git/haulout/ringed/big_gam_ringed_ho_model.RData")  
Covs_ringed = data.frame(sea = factor(rep("Chukchi",n_surveyed),levels=c("Bering","Chukchi")),
                         jday= lubridate::yday(Effort_sp$dt),
                         dmelt_Markus=rep(0,n_surveyed),
                         spring_temp = rep(0,n_surveyed))
for(i in 1:n_surveyed){
  Covs_ringed$dmelt_Markus[i]=Grid_list[[Covs_ringed[i,"jday"]-97]]$dmelt_Markus[Grid_sampled[i]]
  Covs_ringed$spring_temp[i]=Grid_list[[Covs_ringed[i,"jday"]-97]]$spring_temp[Grid_sampled[i]] 
}
Pred_ho_ringed = mgcv::predict.gam(big_model,newdata=Covs_ringed,type="response",se.fit=TRUE)




#thin_rd = 0.65*23/29*ep  #includes disturbance, p
Thin_rd = 132/189*ep*Pred_ho_ringed[[1]]  #includes disturbance, p, haul-out
Sigma_thin = vector("list",2)
Sigma_thin[[1]]=Var_bd  
Thin = c(Thin_bd,Thin_rd)

eb=132/189
varb = eb*(1-eb)/189
epb = ep*eb
varpb =  ep^2*varb + eb^2*varp + varp*varb  #ep * eb
epb2=(eb*ep)^2
EX2 = Pred_ho_ringed$fit^2
Tmp = diag(as.vector(Pred_ho_ringed$fit^2))
for(i in 1:(n.obs-1)){
  for(j in (i+1):n.obs){
    Tmp[i,j]=Tmp[j,i]=Pred_ho_ringed$fit[i]*Pred_ho_ringed$fit[j]
  }
}
Var_rd = Tmp*varpb + (varpb+epb2)*diag(Pred_ho_ringed$se.fit^2)   #temporarily assuming prediction variance independent

Sigma_thin[[2]]=Var_rd
Sigma_thin = as.matrix(bdiag(Sigma_thin))
Sigma_thin = as(Sigma_thin,"dgTMatrix")
#compute logit scale distribution using delta method
diff_logit <- function(x) -1/(x*(x-1))
Thin_logit = log(Thin/(1-Thin))
Diff_logit = diag(diff_logit(Thin))
Sigma_logit_thin = Diff_logit %*% Sigma_thin %*% t(Diff_logit)
diag(Sigma_logit_thin)=diag(Sigma_logit_thin)+0.000001  #adding a  small nugget to  prevent numerical problems


n_species=2
n_transects = n
Mapping = Grid_sampled

#for extra day effect on availability...
X_day = matrix(0,n_transects,2)
X_day[,1] = DayHour_rus[,1]
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2] = (X_day[,1])^2 
X_day = kronecker(diag(n_species),X_day)

S_i = (DayHour_rus[,1]-1)*n_cells+Mapping

# compile habitat covariate data for each cell
n_hab_col=ncol(grid_large@data)

t_steps = length(Grid_list)
Hab_cov = data.frame(matrix(0,n_cells*t_steps,n_hab_col+2))
colnames(Hab_cov)=c(colnames(grid_large@data),"ice2","depth2")
counter=1
for(it in 1:t_steps){
  Tmp_data = as(Grid_list[[it]],"Spatial")@data
  Tmp_data$ice2 = Tmp_data$sea_ice^2
  Hab_cov[counter:(counter+n_cells-1),]=Tmp_data
  counter=counter+n_cells
}
Hab_cov$depth=Hab_cov$depth/mean(Hab_cov$depth)  
Hab_cov$depth2=Hab_cov$depth^2
Hab_cov$I_fast = (Hab_cov$fast_ice>0)
Hab_cov$I_no_ice = (Hab_cov$sea_ice<0.001)
Coords = coordinates(gCentroid(grid_large,byid=TRUE))
Coords = t(t(Coords)/colMeans(Coords))
Hab_cov$easting = rep(Coords[,1],t_steps)
Hab_cov$northing = rep(Coords[,2],t_steps)
X_s = Hab_cov[,c("dist_mainland","depth","sea_ice","fast_ice","I_fast","I_no_ice","easting","northing","snow_depth")]
X_s[,"depth"]=sign(X_s[,"depth"])*abs(X_s[,"depth"])^(1/3)
X_s[,"depth2"]=(X_s[,"depth"])^2
X_s$dist_mainland=X_s$dist_mainland/mean(X_s$dist_mainland)


## format observed data
C_i = Effort_sp@data[,c(10,9,11)]  #just bearded and ringed

#compute proportion of area surveyed using effective strip widths from Column 3, Table 1 of Chernook et al.
Width = 2*tan(44*pi/180)*Effort_sp$Altitude
Area_trans = Width*Effort_sp@data$Length/gArea(grid_large[1,])

#add in pseudo_zeros  - maybe handled with binary ice / no ice covariate instead
# n_0 = 50
# Which_ice0 = which(X_s_2012[,4]==0)
# S_i_2012 = c(S_i_2012,sample(Which_ice0,n_0))
# Which_ice0 = which(X_s_2013[,4]==0)
# S_i_2013 = c(S_i_2013,sample(Which_ice0,n_0))

# Pup_props = c(0.12,0.14,0.12,0.12)   # Pups not modeled separately for CHESS
# Pup_prop_sigma = (diag(0.2*Pup_props))^2
# Diff_logit = diag(diff_logit(Pup_props))
# Sigma_pup_logit = Diff_logit %*% Pup_prop_sigma %*% t(Diff_logit)
# 
# Pups_2012 = BE_2012@data[,13] 
# Pups_2012[Pups_2012=="-"]=0
# Pups_2012=as.numeric(as.character(Pups_2012))
Photo_i = rowSums(Effort_sp@data[,9:10])
Tot_i = rowSums(Effort_sp@data[,9:11])
Nophoto_i = Tot_i - Photo_i
Prop_photo_i = Photo_i/Tot_i
Which_undefined = which(is.na(Prop_photo_i) | Prop_photo_i==Inf)
Prop_photo_i[Which_undefined]=mean(Prop_photo_i[-Which_undefined])


Data_rus = list("C_i"=C_i, "Nophoto_i" = Nophoto_i, "melt" = grid_large$melt_jday_diff,"Prop_photo_i"=Prop_photo_i,"P_i"=Area_trans,"A_s"=1-Hab_cov[,"land_cover"],"S_i"=S_i-1,"X_s"=X_s,"thin_mu_logit"=Thin_logit,"Sigma_logit_thin"=as(Sigma_logit_thin,'dgTMatrix'),"X_day"=X_day,"n_s" = n_cells, "n_sp" = n_species,"n_t"=t_steps)
Data_rus$h_mean = c(mean(HO_out$Mu_bd),mean(Pred_ho_ringed[[1]]))

save(Data_rus,file = "CHESS_russia_Data_TMB_Mar2021.RData")

