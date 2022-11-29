# Abundance estimation for multiple species from count data 
# using spatial regression with prior distributions on detection probability at
# each location sampled.  This version uses splines for env covariate effects

#hopefully, final model where penalty is for ringed seal avail penalty if different from 0.8 at end of study

library( TMB )
library(mgcv)
library(Matrix)
#library( TMBhelper)  #install w/ devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

load('CHESS_data_Nov2021.RData')
# Compile
TmbFile = "./CHESS/src/CHESS_ho_pspline_Feb2021"  #note: edit depending on your working directory
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 

#set up GAMs

Data_df = data.frame(CHESS_data$X_s)
GAM_data = data.frame(Dummy=rep(1,nrow(Data_df)),dist_mainland=Data_df$dist_mainland,
                      depth=Data_df$depth, sea_ice = Data_df$sea_ice,
                      fast_ice=Data_df$fast_ice,I_fast=Data_df$I_fast,
                      I_no_ice=Data_df$I_no_ice,easting=Data_df$easting,
                      northing=Data_df$northing,snow_depth=Data_df$snow_depth)
#gam_setup = gam(Dummy ~ s(dist_mainland, bs = "cs",k=4) + s(sea_ice, bs = "cs",k=6) +
#                  s(fast_ice, bs = "cs",k=6) + s(easting, bs = "cs",k=4) + s(northing,bs="cs",k=3),
#                data = GAM_data,fit=FALSE)
gam_setup = gam(Dummy ~ s(sea_ice, bs = "cs",k=5) +
                  s(fast_ice, bs = "cs",k=5) + s(easting, bs = "cs",k=5) + s(northing,bs="cs",k=5) + s(snow_depth,bs="cs",k=5),
                data = GAM_data,fit=FALSE)
#gam_setup = gam(Dummy ~ s(sea_ice, bs = "cs",k=4) +
#                  s(fast_ice, bs = "cs",k=4) + s(easting, bs = "cs",k=4) + s(northing,bs="cs",k=4),
#                data = GAM_data,fit=FALSE)


#S_dist_mainland = gam_setup$smooth[[1]]$S[[1]]
S_sea_ice = gam_setup$smooth[[1]]$S[[1]]
S_fast_ice = gam_setup$smooth[[2]]$S[[1]]
S_easting = gam_setup$smooth[[3]]$S[[1]]
S_northing = gam_setup$smooth[[4]]$S[[1]]
S_snow_depth = gam_setup$smooth[[5]]$S[[1]]

#S_list = list(S_dist_mainland,S_sea_ice,S_fast_ice,S_easting,S_northing,
#              S_dist_mainland,S_sea_ice,S_fast_ice,S_easting,S_northing) #include twice since 2 species
S_list = list(S_sea_ice,S_fast_ice,S_easting,S_northing,S_snow_depth,
              S_sea_ice,S_fast_ice,S_easting,S_northing,S_snow_depth) #include twice since 2 species
#S_list = list(S_sea_ice,S_fast_ice,S_easting,S_northing,
#              S_sea_ice,S_fast_ice,S_easting,S_northing)

S_combined = .bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S

#For report, used for constructing plots----
#dist_mainland=seq(min(GAM_data$dist_mainland),max(GAM_data$dist_mainland),by = 0.1)
sea_ice = seq(min(GAM_data$sea_ice),max(GAM_data$sea_ice),by = 0.05)
fast_ice = seq(min(GAM_data$fast_ice),max(GAM_data$fast_ice),by = 0.05)
easting = seq(min(GAM_data$easting),max(GAM_data$easting),by = 0.2)
northing = seq(min(GAM_data$northing),max(GAM_data$northing),by = 0.02)
snow_depth = seq(min(GAM_data$snow_depth),max(GAM_data$snow_depth),by = 1)

#landReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(dist_mainland))
seaiceReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(sea_ice))
fasticeReport = PredictMat(gam_setup$smooth[[2]],data = data.frame(fast_ice))
eastingReport = PredictMat(gam_setup$smooth[[3]],data = data.frame(easting))
northingReport = PredictMat(gam_setup$smooth[[4]],data = data.frame(northing))
snowReport =  PredictMat(gam_setup$smooth[[5]],data = data.frame(snow_depth))

#designMatrixForReport = list(landReport,seaiceReport,fasticeReport,eastingReport,northingReport)
designMatrixForReport = list(seaiceReport,fasticeReport,eastingReport,northingReport,snowReport)
#designMatrixForReport = list(seaiceReport,fasticeReport,eastingReport,northingReport)

Melt_df = data.frame(melt=c(CHESS_data$Melt_us_i,CHESS_data$Melt_rus_i))
Melt_df$Dummy = rep(1,nrow(Melt_df))
gam_melt_setup = gam(Dummy ~ s(melt, bs = "cs",k=8),data = Melt_df,fit=FALSE)
S_melt = bdiag(gam_melt_setup$smooth[[1]]$S[[1]])
Sdim_melt = nrow(S_melt) # Find dimension of each S
melt = as.numeric(c(-68:45))
meltReport = PredictMat(gam_melt_setup$smooth[[1]],data = data.frame(melt))
designMatrixForReport_melt = meltReport

load("c:/users/paul.conn/git/BOSSrussia/Covs_ringed_wBS.RData")  #to get western bering availabilty predictions and VC matrix
designMatrix_melt_wBS = PredictMat(gam_melt_setup$smooth[[1]],data = data.frame("melt"=Covs_ringed$dmelt_MDSDA))
load("c:/users/paul.conn/git/BOSSst/Covs_ringed_eBS_2012.RData")  #to get western bering availabilty predictions and VC matrix
designMatrix_melt_eBS_2012 = PredictMat(gam_melt_setup$smooth[[1]],data = data.frame("melt"=Covs_ringed$dmelt_MDSDA))
load("c:/users/paul.conn/git/BOSSst/Covs_ringed_eBS_2013.RData")  #to get western bering availabilty predictions and VC matrix
designMatrix_melt_eBS_2013 = PredictMat(gam_melt_setup$smooth[[1]],data = data.frame("melt"=Covs_ringed$dmelt_MDSDA))
CHESS_data$designMatrix_melt_wBS = designMatrix_melt_wBS
CHESS_data$designMatrix_melt_eBS_2012 = designMatrix_melt_eBS_2012
CHESS_data$designMatrix_melt_eBS_2013 = designMatrix_melt_eBS_2013


CHESS_data = CHESS_data[-which(names(CHESS_data)=="X_s")]  #drop orig design matrix
CHESS_data$S = S_combined
CHESS_data$Sdims = Sdims
CHESS_data$X_s = as.matrix(bdiag(gam_setup$X[,-1],gam_setup$X[,-1])) #one for each species
CHESS_data$designMatrixForReport = .bdiag(designMatrixForReport)
CHESS_data$designMatrixForReport = bdiag(CHESS_data$designMatrixForReport,CHESS_data$designMatrixForReport)
CHESS_data$I_no_ice = rep(Data_df$I_no_ice,2)

n_us_i = length(CHESS_data$Melt_us_i)
n_rus_i = length(CHESS_data$Melt_rus_i)
CHESS_data$S_melt = S_melt
CHESS_data$Sdim_melt = Sdim_melt
CHESS_data$designMatrixForReport_melt = designMatrixForReport_melt
CHESS_data$designMatrix_melt_wBS = designMatrix_melt_wBS
CHESS_data$designMatrix_melt_eBS_2012 = designMatrix_melt_eBS_2012
CHESS_data$designMatrix_melt_eBS_2013 = designMatrix_melt_eBS_2013
CHESS_data$X_ringed_us_i = gam_melt_setup$X[1:n_us_i,-1]
CHESS_data$X_ringed_rus_i= gam_melt_setup$X[(n_us_i+1):(n_rus_i+n_us_i),-1]
CHESS_data$X_ringed_match = CHESS_data$designMatrixForReport_melt

CHESS_data$i_chaunskaya = 0  #indicator for whether or not to model chaunskaya bay  (unstable for bearded)
#CHESS_data$wt_ho = 1000000000.0
CHESS_data$wt_ho = 0
#CHESS_data$wt_melt=10000.0
CHESS_data$wt_ringed = 500.0

#Ho_pred = read.csv("c:/users/paul.conn/git/haulout/ringed/ringed_ho_pred_telem_sea.csv")
#Ho_pred = Ho_pred[Ho_pred$Sea=="Chukchi",]

#load("c:/users/paul.conn/git/haulout/ringed/ringed_ho_pred_telem_high_temp.RData")
#Ho_pred = HO_pred

#CHESS_data$ho_days = c(25:55)
#CHESS_data$ho_match = Ho_pred$Avg_ho[Ho_pred$Jday %in% c(121:151)]
CHESS_data$ho_melt = c(69:114)  #which melt value *entries* to penalize if different than ho_match; melt goes -68:45
#CHESS_data$ho_days = c(1:55)
#CHESS_data$ho_match = Ho_pred$Avg_ho[Ho_pred$Jday %in% c(97:151)]
#CHESS_data$ho_match = Ho_pred$Avg_ho[Ho_pred$Jday %in% c(131:151)]
CHESS_data$ho_match = rep(1.0,length(CHESS_data$ho_melt))  #these need to be the same length!
CHESS_data$ho_fix = 0  #fixes ringed seal availability to haul-out (no lair contribution) 
CHESS_data$day_start = 1  #first day to model counts for ringed seals.  


# Options
#CHESS_data$Options_vec = c("SE"=0)  #bias correction for beta, cell-level intensity?
set.seed(1111)
#CHESS_data$X_s = CHESS_data$X_s[,which(colnames(CHESS_data$X_s) %in% c("depth","depth2"))]
Beta_init = matrix(runif(CHESS_data$n_sp*ncol(CHESS_data$X_s),-2,2),CHESS_data$n_sp,ncol(CHESS_data$X_s))
Params = list("log_N"=log(c(200000,500000)),Beta = rep(0,sum(Sdims)),"thin_beta_day"=c(0,0),
              "phi_log_rus"=rnorm(3),"p_logit_rus"=rnorm(3),  #includes 'NA' as separate count type
              "phi_log_us"=rnorm(3),"p_logit_us"=rnorm(3),"thin_logit_us_i"=CHESS_data$Thin_mu_us_i,
              "thin_logit_rus_i"=CHESS_data$Thin_mu_rus_i,"p_sp_logit" = 3,"log_lambda" = rep(rep(0,length(Sdims))),
              "log_lambda_avail" = 0, "Beta_melt" = rep(0,Sdim_melt),"intercept_avail"=0
              )

# Random
Random= c("Beta","Beta_melt","log_lambda","log_lambda_avail")
#Random=NULL

# Fix parameters
Map = list()
Map[["thin_logit_us_i"]]=factor(rep(NA,length(Params$thin_logit_us_i)))
Map[["thin_logit_rus_i"]]=factor(rep(NA,length(Params$thin_logit_rus_i)))
Map[["thin_beta_day"]]=factor(rep(NA,2))
#Map[["Beta_melt"]]=factor(rep(NA,length(Params$Beta_melt)))
#Map[["intercept_avail"]]=factor(NA)
#Map[["Beta"]]=factor(matrix(c(1:2,NA,NA,3:10,NA,11,NA,12,NA,NA,13:17,NA),nrow=2)) #fix pars to zero for diff covariates
#Map[["Beta"]]=factor(matrix(c(1:2,NA,NA,3:14,NA,NA,15:18,NA,NA),nrow=2)) #fix pars to zero for diff covariates


#Map[["thin_beta_day"]]=factor(rep(NA,length(Params$thin_beta_day)))  #uncomment to fix day effect = 0

#dyn.load( dynlib(TmbFile) )
#Start_time = Sys.time()
#Obj = MakeADFun( data=CHESS_data, parameters=Params, silent=FALSE)
#Obj$fn( Obj$par )





# Make object
dyn.load( dynlib(TmbFile) )
Start_time = Sys.time()
Obj = MakeADFun( data=CHESS_data, parameters=Params, random=Random, map=Map,silent=FALSE)
Obj$fn( Obj$par )
#image(Obj$env$spHess(random=TRUE))  #look at covariance structure

# Perform estimation
start.time = Sys.time()
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))         #
Report = Obj$report()

SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )


#save some stuff helpful for predicting melt-based availability for ringed seals
#note that the traditional approach of computing variance via X * Beta * X' won't work here because of the constraint on melt; variance predictions
#using this approach vary greatly from values output by TMB in sdreport...  instead, use a lookup table approach to formulate VC matrix
#melt_avail_prereqs <- list(gam_melt_setup=gam_melt_setup,intercept_avail=SD$value[which(names(SD$value)=="intercept_avail")],Beta_melt=Report$Beta_melt,
#                           VarCov=SD$cov[397:404,397:404],melt = melt) #note VarCov includes intercept_avail + Beta_melt
Index = which(names(SD$value)=="avail_pred_wBS")
VC_melt=SD$cov[Index,Index]
diag(VC_melt)=diag(VC_melt)+0.000001  #makes it pos-def
melt_avail <- list(Pred_melt=SD$value[Index],VC_melt=VC_melt) #note VarCov includes intercept_avail + Beta_melt
save(melt_avail,file="melt_avail_wBS.RData")

Index = which(names(SD$value)=="avail_pred_eBS_2012")
VC_melt=SD$cov[Index,Index]
diag(VC_melt)=diag(VC_melt)+0.000001  #makes it pos-def
melt_avail <- list(Pred_melt=SD$value[Index],VC_melt=VC_melt) #note VarCov includes intercept_avail + Beta_melt
save(melt_avail,file="melt_avail_eBS_2012.RData")

Index = which(names(SD$value)=="avail_pred_eBS_2013")
VC_melt=SD$cov[Index,Index]
diag(VC_melt)=diag(VC_melt)+0.000001  #makes it pos-def
melt_avail <- list(Pred_melt=SD$value[Index],VC_melt=VC_melt) #note VarCov includes intercept_avail + Beta_melt
save(melt_avail,file="melt_avail_eBS_2013.RData")


#### plotting, etc.
Plot_df1 = data.frame(Melt = melt,Type = c(rep("Not fitted",68),rep("Fitted",46)),
                     Proportion = c(SD$value[which(names(SD$value)=="splineForReport_melt")]))
#Plot_df2 = data.frame(Melt = rep(c(97:151),1),
#                      Proportion = SD$value[185:239])
library(ggplot2)
Avail_plot = ggplot()+geom_line(data=Plot_df1,aes(x=Melt,y=Proportion),size=1)+ theme(text=element_text(size=14))+xlab("Day from snowmelt onset")
Avail_plot
pdf("Avail_plot.pdf")            
  Avail_plot
dev.off()

png("Avail_plot.png")
  Avail_plot
dev.off()
                                                                                            


#plot estimates for Russia and U.S.
Survey_day = c(1:55)
load("I_US.Rdata")
I_us = rep(0,1354)
I_us[I_US]=1
Bearded = Ringed = matrix(0,2,55)
for(it in 1:55){
  Bearded[1,it]=sum(I_us*Report$Z_s[1,(1354*(it-1)+1):(1354*(it-1)+1354)])
  Bearded[2,it]=sum((1-I_us)*Report$Z_s[1,(1354*(it-1)+1):(1354*(it-1)+1354)])
  Ringed[1,it]=sum(I_us*Report$Z_s[2,(1354*(it-1)+1):(1354*(it-1)+1354)])
  Ringed[2,it]=sum((1-I_us)*Report$Z_s[2,(1354*(it-1)+1):(1354*(it-1)+1354)])
}
par(mfrow=c(2,1))
plot(Survey_day,Bearded[1,],ylim=c(0,140000))
lines(Survey_day,Bearded[2,])
plot(Survey_day,Ringed[1,],ylim=c(0,700000))
lines(Survey_day,Ringed[2,])

mean(Bearded[1,])
mean(Bearded[1,])/(mean(colSums(Bearded)))  #65% in US waters mean(Bearded[1,])
mean(Ringed[1,])
mean(Ringed[1,])/(mean(colSums(Ringed)))  #80% in US waters mean(Ringed[1,])

#estimates for Irina
load("Chess_grid_list_all.RDa")
for(it in 1:CHESS_data$n_t){
  Grid_list[[it]]$bearded_est = Report$Z_s[1,((it-1)*CHESS_data$n_s+1):((it-1)*CHESS_data$n_s+CHESS_data$n_s)]
  Grid_list[[it]]$ringed_est = Report$Z_s[2,((it-1)*CHESS_data$n_s+1):((it-1)*CHESS_data$n_s+CHESS_data$n_s)]
}
save(Grid_list,file="CHESS_seal_ests_for_Irina_2022.RData")




Var_boot = c(224434357,2405072570)  #obtained by running run_CHESS_TMB_ho_pspline_boot.R, etc. then combine_boot_var.R

SE.bd = sqrt(SD$sd[length(SD$sd)-1]^2+Var_boot[1])
SE.rd = sqrt(SD$sd[length(SD$sd)]^2+Var_boot[2])

C = exp(1.96 * sqrt(log(1+(SE.bd/Report$N[1])^2)))
CI_lo<- Report$N[1]/C
CI_hi<- Report$N[1]*C
cat(paste("Bearded: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))

C = exp(1.96 * sqrt(log(1+(SE.rd/Report$N[2])^2)))
CI_lo<- Report$N[2]/C
CI_hi<- Report$N[2]*C
cat(paste("Ringed: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))


library(sf)
library(ggplot2)
plot_N_map<-function(cur.t,N,Grid,No.plot=NA){
  #require(rgeos)
  Tmp=Grid[[1]]
  Abundance=N[,cur.t]
  Cur.df=cbind(data.frame((st_coordinates(st_centroid(Tmp,byid=TRUE)))),Abundance)
  Cur.df[,1]=as.integer(Cur.df[,1])
  Cur.df[,2]=as.integer(Cur.df[,2])
  if(! is.na(No.plot)){
    Cur.df=Cur.df[-No.plot,]
  }
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(x=Easting,y=Northing,fill=Abundance)+geom_raster()+tmp.theme
  p1 = p1 + scale_fill_viridis_c()
  p1
}
load("Chess_grid_list_all.RDa")
Bay_mask = c(955:957,995:997,1030:1032,1060:1062)

#fill Z_s 
isp =2
N_s = matrix(Report$Z_s[isp,],CHESS_data$n_s,CHESS_data$n_t)
#N_s[,1]=CHESS_data$A_s[1:CHESS_data$n_s]

plot_N_map(1,N_s,Grid_list,No.plot = Bay_mask)
plot_N_map(41,N_s,Grid_list)

#plot expected vs observed counts (US)
plot(Report$E_count_obs_us_i[,2],CHESS_data$C_us[,2])

plot_N_map(35,log(N_s),Grid_list)
#time averaged estimate
for(it in 1:55){
  N_s[,1] = rowMeans(N_s)
}
Bearded_plot = plot_N_map(1,N_s,Grid_list,No.plot = Bay_mask)+ggtitle("A. Bearded seals")

isp =2
N_s = matrix(Report$Z_s[isp,],CHESS_data$n_s,CHESS_data$n_t)
plot_N_map(1,N_s,Grid_list,No.plot = Bay_mask)
plot_N_map(35,log(N_s),Grid_list)
#time averaged estimates
for(it in 1:55){
  N_s[,1] = rowMeans(N_s)
}
Ringed_plot = plot_N_map(1,N_s,Grid_list,No.plot = Bay_mask)+ggtitle("B. Ringed seals")

library(gridExtra)
p1 = ggplotGrob(Bearded_plot)
p2 = ggplotGrob(Ringed_plot)
pdf("Seal_N_plot.pdf",width=4,height=6)
  grid.arrange(p1,p2,ncol=1)
  #cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()

png("Seal_N_plot.png",width=480,height=720)
grid.arrange(p1,p2,ncol=1)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()

#Plot time-specific estimates; rows: ice, bearded, ringed, cols: Apr 7, May 3, May 31
It = c(1,27,54) 
library(RColorBrewer)
tmp_theme = themecolorRampPalette(brewer.pal(6, "Blues"))
Bearded_plots=Ringed_plots=Ice_plots=vector("list",3)
myPalette = colorRampPalette(brewer.pal(6, "Blues"))
Titles = c("7 Apr","3 May","31 May")
for(it in 1:3){
  N_s= matrix(Report$Z_s[1,],CHESS_data$n_s,CHESS_data$n_t)
  Bearded_plots[[it]]=plot_N_map(It[it],N_s,Grid_list,No.plot = Bay_mask)
  Bearded_plots[[it]]=Bearded_plots[[it]]+labs(fill='Bd N')+
    theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title = element_blank())+
    scale_fill_viridis_c(limits = c(0,6100), breaks = c(0, 2000, 4000, 6000),values=c(0,.1,.2,1))
  
  N_s= matrix(Report$Z_s[2,],CHESS_data$n_s,CHESS_data$n_t)
  Ringed_plots[[it]]=plot_N_map(It[it],N_s,Grid_list,No.plot = Bay_mask)
  Ringed_plots[[it]]=Ringed_plots[[it]]+labs(fill='Rd N')
  Ringed_plots[[it]]=Ringed_plots[[it]]+theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
  Ringed_plots[[it]]=Ringed_plots[[it]]+scale_fill_viridis_c(limits = c(0,23000), breaks = c(0, 6000, 12000, 18000),values=c(0,.01,.2,1))
  
  N_s= matrix(Grid_list[[It[it]]]$sea_ice,CHESS_data$n_s,CHESS_data$n_t)
  Ice_plots[[it]]=plot_N_map(It[it],N_s,Grid_list,No.plot = Bay_mask)
  Ice_plots[[it]]=Ice_plots[[it]]+scale_fill_gradientn(colours=rev(myPalette(100)),name="Ice",limits=c(0,1))+
    theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle(Titles[it])
  
}
pdf("Seal_N_st_plot.pdf",height=6,width=8)
grid.arrange(Ice_plots[[1]],Ice_plots[[2]],Ice_plots[[3]],Bearded_plots[[1]],Bearded_plots[[2]],
             Bearded_plots[[3]], Ringed_plots[[1]], Ringed_plots[[2]],Ringed_plots[[3]],ncol=3)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()

png("Seal_N_st_plot.png",height=640,width=800)
grid.arrange(Ice_plots[[1]],Ice_plots[[2]],Ice_plots[[3]],Bearded_plots[[1]],Bearded_plots[[2]],
             Bearded_plots[[3]], Ringed_plots[[1]], Ringed_plots[[2]],Ringed_plots[[3]],ncol=3)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()

#Data for Marilyn Myers
# N_s= matrix(Report$Z_s[1,],CHESS_data$n_s,CHESS_data$n_t)
# N_bd = rowMeans(N_s[,45:55])
# N_s= matrix(Report$Z_s[2,],CHESS_data$n_s,CHESS_data$n_t)
# N_rd = rowMeans(N_s[,45:55])
# Abund = Grid_list[[55]][,"land_cover"]
# Abund$N_bd = N_bd
# Abund$N_rd = N_rd
# Abund$D_bd = N_bd / (625*(1-Abund$land_cover))
# Abund$D_rd = N_rd / (625*(1-Abund$land_cover))
# save(Abund,file="Seal_abundance_Chukchi_May2016.RData")
# st_write(Abund,"Seal_abundance_Chukchi_May2016.shp")

#plot counts
plot_N_map_xy<-function(N,XY,leg.title="Abundance"){
  require(ggplot2)
  library(RColorBrewer)
  myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  Abundance=N
  Cur.df=cbind(data.frame(Easting=XY[,1],Northing=XY[,2],Abundance))
  colnames(Cur.df)=c("Easting","Northing","Abundance")
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  p1
}
Coords = st_coordinates(st_centroid(Grid_list[[1]]))
Coords2=Coords
Coords2[,1] = round(Coords2[,1])
Coords2[,2] = round(Coords2[,2])
Which_after_us = which(CHESS_data$Day_us_i<34)
Which_after_rus = which(CHESS_data$Day_rus_i<34)
#Count = log(c(CHESS_data$C_us[Which_after_us,2],CHESS_data$C_rus[Which_after_rus,2])+1)
Count = c(CHESS_data$C_us[Which_after_us,2],CHESS_data$C_rus[Which_after_rus,2])

Mapping = c(CHESS_data$S_us_i[Which_after_us] %% CHESS_data$n_s,CHESS_data$S_rus_i[Which_after_rus] %% CHESS_data$n_s)
Coords3 = Coords2[Mapping,]
plot_N_map_xy(N=Count,XY=Coords3,leg.title="Count")  #note that when a cell is surveyed >1 time, only one value is plotted here





##GOF
set.seed(11111)
library(tweedie)
Resids = 0*Report$E_count_obs_us_i
for(icol in 1:3){
  Resids[,icol] = ptweedie(CHESS_data$C_us[,icol]-1,mu=Report$E_count_obs_us_i[,icol],phi=Report$phi_us[icol],power=Report$power_us[icol])+runif(nrow(CHESS_data$C_us))*dtweedie(CHESS_data$C_us[,icol],mu=Report$E_count_obs_us_i[,icol],phi=Report$phi_us[icol],power=Report$power_us[icol])
}

Resids[Resids>1]=0.999
Resids_US = Resids

Resid_binned = matrix(0,10,3)
for(irow in 1:nrow(Resids)){
  for(icol in 1:3){
    Resid_binned[ceiling(Resids[irow,icol]/0.1),icol]=Resid_binned[ceiling(Resids[irow,icol]/0.1),icol]+1
  }
}
Xsq = rep(0,3)
for(i in 1:3){
  Xsq[i]=10/nrow(Resids)*sum((Resid_binned[,i]-nrow(Resids)/10)^2)
}
Pval = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Pval

Resids.df = data.frame(Residual = as.vector(Resids))
Labels1 = rep('',3)
for(i in 1:3){
  Labels1[i] = paste0('Obs = ',i) #,', p=',format(Pval[i],digits=2))
}
Labels1 = factor(Labels1,levels=Labels1)
Resids.df$Labels=rep(Labels1,each=nrow(Resids))

# myplot=ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:10)/10)+facet_wrap(~Labels)+ylab('Frequency')+xlab('Randomized quantile residual')
# pdf('GOF_CHESS_us.pdf')
# myplot
# dev.off()

##russia GOF
Resids = 0*Report$E_count_obs_rus_i
Mu = Report$E_count_obs_rus_i
Mu[which(Mu==0)]=0.00000001
for(icol in 1:3){
  Resids[,icol] = ptweedie(CHESS_data$C_rus[,icol]-1,mu=Mu[,icol],phi=Report$phi_rus[icol],power=Report$power_rus[icol])+runif(nrow(CHESS_data$C_rus))*dtweedie(CHESS_data$C_rus[,icol],mu=Mu[,icol],phi=Report$phi_rus[icol],power=Report$power_rus[icol])
}

Resids[Resids>1]=0.999

Resid_binned = matrix(0,10,3)
for(irow in 1:nrow(Resids)){
  for(icol in 1:3){
    Resid_binned[ceiling(Resids[irow,icol]/0.1),icol]=Resid_binned[ceiling(Resids[irow,icol]/0.1),icol]+1
  }
}
Xsq = rep(0,3)
for(i in 1:3){
  Xsq[i]=10/nrow(Resids)*sum((Resid_binned[,i]-nrow(Resids)/10)^2)
}
Pval = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Pval

Resids.df2 = data.frame(Residual = as.vector(Resids))
Labels1 = rep('',3)
for(i in 1:3){
  Labels1[i] = paste0('Obs = ',i) #,', p=',format(Pval[i],digits=2))
}
Labels1 = factor(Labels1,levels=Labels1)
Resids.df2$Labels=rep(Labels1,each=nrow(Resids))

Resids.df = rbind(Resids.df,Resids.df2)
Resids.df$country = c(rep("U.S.",3*nrow(CHESS_data$C_us)),rep("Russia",3*nrow(CHESS_data$C_rus)))

myplot=ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:10)/10)+facet_grid(country~Labels)+ylab('Frequency')+xlab('Randomized quantile residual')
pdf('GOF_CHESS.pdf')
myplot
dev.off()

png('GOF_CHESS.png')
myplot
dev.off()



# residual mapsplot_N_map_xy<-function(N,XY,leg.title="Abundance"){
#devtools::install_github("kwstat/pals")
plot_N_map_xy<-function(N,XY,leg.title="Abundance"){
  require(ggplot2)
  library(scales)
  #library(pals)
  #myPalette <- ocean.balance
  Abundance=N
  Cur.df=cbind(data.frame(Easting=XY[,1],Northing=XY[,2],Abundance))
  colnames(Cur.df)=c("Easting","Northing","Abundance")
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradient2(low=muted("blue"),mid="white",midpoint=0,high=muted("red"),name=leg.title)
  p1
}

Which_counts_early_us = which(CHESS_data$Day_us_i<19)
Which_counts_middle_us = which(CHESS_data$Day_us_i>=19 & CHESS_data$Day_us_i<=38)
Which_counts_late_us = which(CHESS_data$Day_us_i>45)
Which_counts_early_rus = which(CHESS_data$Day_rus_i<19)
Which_counts_middle_rus = which(CHESS_data$Day_rus_i>=19 & CHESS_data$Day_rus_i<=38)
Which_counts_late_rus = which(CHESS_data$Day_rus_i>38)
Resids_us = matrix(0,nrow(Report$E_count_obs_us_i),3)
for(isp in 1:3){
  Resids_us[,isp] = -log(dtweedie(CHESS_data$C_us[,isp],mu=Report$E_count_obs_us_i[,isp],phi=Report$phi_us[isp],power=Report$power_us[isp]))
  Resids_us[,isp] = Resids_us[,isp]*-sign(Report$E_count_obs_us_i[,isp]-CHESS_data$C_us[,isp])  #make underpredictions positive
}
Resids_rus = matrix(0,nrow(Report$E_count_obs_rus_i),3)
for(isp in 1:3){
  Mu_rus = Report$E_count_obs_rus_i[,isp]
  Mu_rus[which(Mu_rus==0)]=0.0000001
  Resids_rus[,isp] = -log(dtweedie(CHESS_data$C_rus[,isp],mu=Mu_rus,phi=Report$phi_rus[isp],power=Report$power_rus[isp]))
  Resids_rus[,isp] = Resids_rus[,isp]*-sign(Report$E_count_obs_rus_i[,isp]-CHESS_data$C_rus[,isp])  
}
#list format for residual plots is 3x3 e.g. Resid_maps[[species]][[time period]]
Resid_maps <- vector("list",3)
for(i in 1:3){
  Resid_maps[[i]]<-vector("list",3)
}

#Titles=c("Early","Middle","Late")
library(gridExtra)
library(grid)
for(isp in 1:3){
  for(it in 1:3){
    if(it==1){
      Lresid <- c(Resids_us[Which_counts_early_us,isp],Resids_rus[Which_counts_early_rus,isp])
      Mapping = c(CHESS_data$S_us_i[Which_counts_early_us] %% CHESS_data$n_s,CHESS_data$S_rus_i[Which_counts_early_rus] %% CHESS_data$n_s)
    }
    if(it==2){
      Lresid <- c(Resids_us[Which_counts_middle_us,isp],Resids_rus[Which_counts_middle_rus,isp])
      Mapping = c(CHESS_data$S_us_i[Which_counts_middle_us] %% CHESS_data$n_s,CHESS_data$S_rus_i[Which_counts_middle_rus] %% CHESS_data$n_s)
    }
    if(it==3){
      Lresid <- c(Resids_us[Which_counts_late_us,isp],Resids_rus[Which_counts_late_rus,isp])
      Mapping = c(CHESS_data$S_us_i[Which_counts_late_us] %% CHESS_data$n_s,CHESS_data$S_rus_i[Which_counts_late_rus] %% CHESS_data$n_s)
    }
    Coords3 = Coords2[Mapping,]
    cur.plot = plot_N_map_xy(N=Lresid,XY=Coords3,leg.title="Score")
    #if(isp==1)cur.plot=cur.plot+ggtitle(Titles[isp])
    Resid_maps[[isp]][[it]] = ggplotGrob(cur.plot)  #note that when a cell is surveyed >1 time, only one value is plotted here
  }
}

Resid_map_list <- vector("list",9)
counter=1
for(isp in 1:3){
  for(it in 1:3){
    Resid_map_list[[counter]]=Resid_maps[[isp]][[it]]
    counter = counter+1
  }
}
tt <- ttheme_default(base_size = 12,
                     core = list(fg_params=list(rot=90)))
combine <- rbind(tableGrob(t(c("Early","Middle","Late")), theme = ttheme_default(), rows = ""), 
                 cbind(tableGrob(c("Bearded","Ringed","Unknown"), theme = tt), 
                       arrangeGrob(grobs = Resid_map_list),  size = "last"), size = "last")

grid.newpage()
grid.draw(combine)


grid.arrange(Resid_maps[[1]][[1]],Resid_maps[[1]][[2]],Resid_maps[[1]][[3]],Resid_maps[[1]][[2]],Resid_maps[[2]][[2]],Resid_maps[[2]][[3]],Resid_maps[[3]][[1]],Resid_maps[[3]][[2]],Resid_maps[[3]][[3]],ncol=3)
  
grid.arrange()

library(gridExtra)
p1 = ggplotGrob(Bearded_plot)
p2 = ggplotGrob(Ringed_plot)
pdf("Seal_N_plot.pdf",width=4,height=6)
grid.arrange(p1,p2,ncol=1)
#cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
dev.off()




### overlay Bengtson et al. survey area
load('c:/users/paul.conn/git/chukchipower/Grid_chukchi_ExpDensities.Rda')
Grid_chukchi = st_as_sf(Grid_chukchi)
Grid_chukchi = st_centroid(st_transform(Grid_chukchi,st_crs(Grid_list[[1]])))
I.intersect = rowSums(st_intersects(Grid_list[[1]],Grid_chukchi,sparse=FALSE))
sum(Report$Z_s[1,which(I.intersect>0)])   #number in U.S. at beginning of study
sum(Report$Z_s[2,which(I.intersect>0)])   #number in U.S. study area at beginning of study

start=73117
end = 74470
b_mean = 0
r_mean = 0
for(it in 0:10){
  Which = c(start:end)[I.intersect>0]
  b_mean = b_mean + sum(Report$Z_s[1,Which])   #number in Bengtson study area at end of study
  r_mean = r_mean + sum(Report$Z_s[2,Which])   #number in Bengtson study area at end of study
  start = start-CHESS_data$n_s
  end = end-CHESS_data$n_s
}
b_mean = b_mean/11
r_mean = r_mean/11

b_mean
#use bootstrap CV to get 95% CI
CV_bd = SE.bd / Report$N[1]
C = exp(1.96 * sqrt(log(1+(CV_bd)^2)))
CI_lo<- b_mean/C
CI_hi<- b_mean*C
cat(paste("Bearded: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))

r_mean
CV_rd = SE.rd / Report$N[2]
C = exp(1.96 * sqrt(log(1+(CV_rd)^2)))
CI_lo<- r_mean/C
CI_hi<- r_mean*C
cat(paste("Ringed: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))


Cells_us = unique(CHESS_data$S_us_i %% CHESS_data$n_s)
Cells_rus = unique(CHESS_data$S_rus_i %% CHESS_data$n_s)
Cells_samp = c(Cells_us,Cells_rus)

Cells = matrix(rep(0,CHESS_data$n_s),ncol=1)
Cells[Cells_samp]=1
#Cells[Cells_us]=1
#Cells[1319]=1
#Cells[c(667,712,957,1062,1215,1242,1319,1335,1345,1346,1347)]=1
Cells = matrix(CHESS_data$X_s[1:CHESS_data$n_s,"depth"],ncol=1)
plot_N_map(1,Cells,Grid_list)
#1242, 1319


#gIVH
SD_Z0 = SD$sd[25:2732]
Which_obs = c(CHESS_data$S_rus_i,CHESS_data$S_us_i)%%CHESS_data$n_s
max_obs_sd = max(SD_Z0[Which_obs])
which(SD_Z0 > max_obs_sd)

#splines
muSpline = SD$value[names(SD$value)=="splineForReport"]
sdSpline<-SD$sd[names(SD$value)=="splineForReport"]

n_entries = length(muSpline)
Plot.df = data.frame(Covariate=rep(c(rep("sea_ice",length(sea_ice)),
                                   rep("fast_ice",length(fast_ice)),rep("easting",length(easting)),
                                   rep("northing",length(northing)),rep("snow_depth",length(snow_depth))),2),
                     x = rep(c(sea_ice,fast_ice,easting,northing,snow_depth),2),
                     Smooth = muSpline, Lower = muSpline-1.96*sdSpline, Upper = muSpline+1.96*sdSpline,
                     Species = factor(c(rep("Bearded",n_entries/2),rep("ringed",n_entries/2))))
Plot.df.bearded=Plot.df[1:(n_entries/2),]

bplot = ggplot()+geom_line(data=Plot.df.bearded,aes(x=x,y=Smooth))+
         geom_ribbon(data=Plot.df.bearded,aes(x=x,ymin=Lower,ymax=Upper),alpha=0.4)+
         facet_wrap(~Covariate,scales="free")+xlab('')
pdf('bearded_smooths.pdf')
  bplot
dev.off()

png('bearded_smooths.png')
bplot
dev.off()

Plot.df.ringed=Plot.df[(n_entries/2+1):n_entries,]

bplot = ggplot()+geom_line(data=Plot.df.ringed,aes(x=x,y=Smooth))+
  geom_ribbon(data=Plot.df.ringed,aes(x=x,ymin=Lower,ymax=Upper),alpha=0.4)+
  facet_wrap(~Covariate,scales="free")+xlab('')
pdf('ringed_smooths.pdf')
bplot
dev.off()

png('ringed_smooths.png')
bplot
dev.off()


#quick GAM mod to look at time effect
GAM_df = data.frame(easting=Data_df[CHESS_data$S_us_i,"easting"],
                    northing=Data_df[CHESS_data$S_us_i,"northing"], ice = Data_df[CHESS_data$S_us_i,"sea_ice"],
                    day = CHESS_data$Day_us_i,count = CHESS_data$C_us[,2],
                    P = CHESS_data$P_us_i)
day_gam = gam(count ~ s(easting)+s(ice)+s(northing)+s(day),offset=log(P),family="poisson",data=GAM_df)
plot(day_gam)

New_data = data.frame(easting = rep(GAM_df$easting[1],55),northing = rep(GAM_df$northing[1],55),day = c(1:55),ice=rep(GAM_df$ice[1],55))
plot(predict(day_gam,newdata=New_data,type="response"))
               
Pred_df = predict(day_gam)


Est_CHESS = vector("list",3)
names(Est_CHESS)=c("meta","bearded","ringed")
Est_CHESS$meta=list(start_date="April 7",end_date="May 31")
Est_CHESS$ringed = Est_CHESS$bearded = matrix(0,CHESS_data$n_s,CHESS_data$n_t)
for(it in 1:CHESS_data$n_t){
  Est_CHESS$bearded[,it] = Report$Z_s[1,(it-1)*CHESS_data$n_s+c(1:CHESS_data$n_s)]
  Est_CHESS$ringed[,it] = Report$Z_s[2,(it-1)*CHESS_data$n_s+c(1:CHESS_data$n_s)]
}
Est_CHESS$meta$grid = Grid_list[[1]]
save(Est_CHESS,file="Ests_CHESS.RData")

#attach estimates to grid and output for Erin R. to make maps
for(i in 1:CHESS_data$n_t){
  Grid_list[[i]]$bearded = Report$Z_s[1,(it-1)*CHESS_data$n_s+c(1:CHESS_data$n_s)]
  Grid_list[[i]]$ringed = Report$Z_s[2,(it-1)*CHESS_data$n_s+c(1:CHESS_data$n_s)]
}
save(Grid_list,file="Estimate_ice_data_for_ErinR.RData")
