# Abundance estimation for multiple species from count data 
# using spatial regression with prior distributions on detection probability at
# each location sampled.  This version uses splines for env covariate effects

#hopefully, final model where penalty is for ringed seal avail penalty if different from 0.8 at end of study

library( TMB )
library(mgcv)
library(Matrix)

load('CHESSseal/CHESS_data_Nov2021.RData')

# Compile
TmbFile = "./CHESSseal/CHESS/src/CHESS_BovengEtAl2025cpp"  #note: edit depending on your working directory
compile(paste0(TmbFile,".cpp"),DLLFLAGS="") 

#set up GAMs
Data_df = data.frame(CHESS_data$X_s)
GAM_data = data.frame(Dummy=rep(1,nrow(Data_df)),dist_mainland=Data_df$dist_mainland,
                      depth=Data_df$depth, sea_ice = Data_df$sea_ice,
                      fast_ice=Data_df$fast_ice,I_fast=Data_df$I_fast,
                      I_no_ice=Data_df$I_no_ice,easting=Data_df$easting,
                      northing=Data_df$northing,snow_depth=Data_df$snow_depth)

gam_setup = gam(Dummy ~ s(sea_ice, bs = "cs",k=5) +
                  s(fast_ice, bs = "cs",k=5) + s(easting, bs = "cs",k=5) + s(northing,bs="cs",k=5) + s(snow_depth,bs="cs",k=5),
                data = GAM_data,fit=FALSE)
S_sea_ice = gam_setup$smooth[[1]]$S[[1]]
S_fast_ice = gam_setup$smooth[[2]]$S[[1]]
S_easting = gam_setup$smooth[[3]]$S[[1]]
S_northing = gam_setup$smooth[[4]]$S[[1]]
S_snow_depth = gam_setup$smooth[[5]]$S[[1]]

S_list = list(S_sea_ice,S_fast_ice,S_easting,S_northing,S_snow_depth,
              S_sea_ice,S_fast_ice,S_easting,S_northing,S_snow_depth) #include twice since 2 species

S_combined = .bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S

#For report, used for constructing plots----
sea_ice = seq(min(GAM_data$sea_ice),max(GAM_data$sea_ice),by = 0.05)
fast_ice = seq(min(GAM_data$fast_ice),max(GAM_data$fast_ice),by = 0.05)
easting = seq(min(GAM_data$easting),max(GAM_data$easting),by = 0.2)
northing = seq(min(GAM_data$northing),max(GAM_data$northing),by = 0.02)
snow_depth = seq(min(GAM_data$snow_depth),max(GAM_data$snow_depth),by = 1)

seaiceReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(sea_ice))
fasticeReport = PredictMat(gam_setup$smooth[[2]],data = data.frame(fast_ice))
eastingReport = PredictMat(gam_setup$smooth[[3]],data = data.frame(easting))
northingReport = PredictMat(gam_setup$smooth[[4]],data = data.frame(northing))
snowReport =  PredictMat(gam_setup$smooth[[5]],data = data.frame(snow_depth))

designMatrixForReport = list(seaiceReport,fasticeReport,eastingReport,northingReport,snowReport)

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
CHESS_data$X_ringed_us_i = gam_melt_setup$X[1:n_us_i,-1]
CHESS_data$X_ringed_rus_i= gam_melt_setup$X[(n_us_i+1):(n_rus_i+n_us_i),-1]
CHESS_data$X_ringed_match = CHESS_data$designMatrixForReport_melt

CHESS_data$i_chaunskaya = 0  #indicator for whether or not to model chaunskaya bay  (unstable for bearded)
CHESS_data$wt_ho = 0
CHESS_data$wt_ringed = 500.0

CHESS_data$ho_melt = c(69:114)  #which melt value *entries* to penalize if different than ho_match; melt goes -68:45
CHESS_data$ho_match = rep(1.0,length(CHESS_data$ho_melt))  #these need to be the same length!
CHESS_data$ho_fix = 0  #fixes ringed seal availability to haul-out (no lair contribution) 
CHESS_data$day_start = 1  #first day to model counts for ringed seals.  


# Options
set.seed(1111)
Beta_init = matrix(runif(CHESS_data$n_sp*ncol(CHESS_data$X_s),-2,2),CHESS_data$n_sp,ncol(CHESS_data$X_s))
Params = list("log_N"=log(c(200000,500000)),Beta = rep(0,sum(Sdims)),"thin_beta_day"=c(0,0),
              "phi_log_rus"=rnorm(3),"p_logit_rus"=rnorm(3),  #includes 'NA' as separate count type
              "phi_log_us"=rnorm(3),"p_logit_us"=rnorm(3),"thin_logit_us_i"=CHESS_data$Thin_mu_us_i,
              "thin_logit_rus_i"=CHESS_data$Thin_mu_rus_i,"p_sp_logit" = 3,"log_lambda" = rep(rep(0,length(Sdims))),
              "log_lambda_avail" = 0, "Beta_melt" = rep(0,Sdim_melt),"intercept_avail"=0
              )

# Random
Random= c("Beta","Beta_melt","log_lambda","log_lambda_avail")

# Fix parameters
Map = list()
Map[["thin_logit_us_i"]]=factor(rep(NA,length(Params$thin_logit_us_i)))
Map[["thin_logit_rus_i"]]=factor(rep(NA,length(Params$thin_logit_rus_i)))
Map[["thin_beta_day"]]=factor(rep(NA,2))


# Make object
dyn.load( dynlib(TmbFile) )
Start_time = Sys.time()
Obj = MakeADFun( data=CHESS_data, parameters=Params, random=Random, map=Map,silent=FALSE)
Obj$fn( Obj$par )

# Perform estimation
start.time = Sys.time()
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))         #
Report = Obj$report()

SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )



#### plotting, etc.
Plot_df1 = data.frame(Melt = melt,Type = c(rep("Not fitted",68),rep("Fitted",46)),
                     Proportion = c(SD$value[which(names(SD$value)=="splineForReport_melt")]))
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
