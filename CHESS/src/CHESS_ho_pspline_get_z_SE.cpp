// this version just uses first element of ringed seal thinning parameters for US
// i.e. just one parameter for asymptote of availability model
#include <TMB.hpp>

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

template<class Type>
Type dbern(Type x, Type prob, int give_log=1){
  Type logres;
  if( x==0 ) logres = log( 1-prob );
  if( x==1 ) logres = log( prob );
  if(give_log) return logres; else return exp(logres);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace Eigen;
  using namespace density;
  
  // options vec
  //DATA_FACTOR( Options_vec );
  // Slot 0: compute SE? 
  
  // Data
  DATA_MATRIX( C_rus );       	// Matrix of responses (counts) of each species at each sampled location - Russian format
  DATA_MATRIX( C_us);      // same, US data format
  DATA_VECTOR( P_rus_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( P_us_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_IVECTOR( S_rus_i ); // Site-time index for each sample
  DATA_IVECTOR( S_us_i ); // Site-time index for each sample
  DATA_IVECTOR( Day_rus_i ); // day-of-survey (russia)
  DATA_IVECTOR( Day_us_i ); // day-of-survey
  
  
  //DATA_IMATRIX( Y_s); //indicator for which sites sampled/not sampled
  DATA_MATRIX( X_s );  //design matrix for smooths
  DATA_VECTOR(I_no_ice); //use to enforce zero abundance when ice = 0
  
  DATA_VECTOR(Thin_mu_rus_i);  //logit scale mean params 
  DATA_VECTOR(Thin_mu_us_i);   // same, US
  DATA_SPARSE_MATRIX(Thin_sigma_rus); //var-cov of logit scale thining params (i.e. Thin_mu_rus)
  DATA_SPARSE_MATRIX(Thin_sigma_us);
  
  DATA_MATRIX(X_day_rus);  // design matrix for extra day and day^2 effects
  DATA_MATRIX(X_day_us);  // design matrix for extra day and day^2 effects
  
  DATA_INTEGER(n_s); //number of cells
  DATA_INTEGER(n_sp); //number of species
  DATA_INTEGER(n_t); //number of time steps
  
  DATA_SCALAR(h_mean_rus);  //mean of haulout distributions for bearded seals surveyed in Russia
  DATA_SCALAR(h_mean_us);  //mean of haulout distributions for bearded seals surveyed in US
  
  DATA_VECTOR(Nophoto_i);  //for russian data, number of observations in each sampled cell that don't have photos
  DATA_VECTOR(Prop_photo_i); //for russian data, the proportion of cells in each strata that don't have data
  DATA_VECTOR(Melt_rus_i); // Difference in day from MSMSDA values of sea ice melt
  DATA_VECTOR(Melt_us_i);
  
  DATA_SPARSE_MATRIX(S);//Penalization block diagonal matrix 
  DATA_IVECTOR(Sdims);   //Dimensions of each smooth
  DATA_SPARSE_MATRIX(designMatrixForReport);//Design matrix for report of splines
  
  DATA_SPARSE_MATRIX(S_melt);
  DATA_INTEGER(Sdim_melt);
  DATA_MATRIX(designMatrixForReport_melt); 
  DATA_MATRIX(designMatrix_melt_wBS);
  DATA_MATRIX(designMatrix_melt_eBS_2012);
  DATA_MATRIX(designMatrix_melt_eBS_2013);
  DATA_MATRIX(X_ringed_us_i); //design matrix for ringed seal availability smooths (US)
  DATA_MATRIX(X_ringed_rus_i);
  DATA_MATRIX(X_ringed_match);  // for penalty 
    
  
  DATA_IVECTOR(Which_cells_bearded);  //initial model fitting indicated issues with Chaunskaya Bay - this is to produce estimates omitting it
  DATA_VECTOR(I_cell); //vector of zero's and one's indicating which cells can have positive abundance (use e.g. to set Chaunskaya to zero) 
 
  DATA_SCALAR(wt_ho);   //weight on haulout for bearded differeng from estimated mean
  DATA_SCALAR(wt_ringed);  //weight on penalty for max avail on last day of study
  DATA_IVECTOR(ho_melt); // which days of survey to penalize if different from ho_match
  DATA_VECTOR(ho_match); 
  DATA_INTEGER(ho_fix); //if 1, fix ringed seal haul-out values to ho_match
  
  DATA_INTEGER(day_start);  // day number to start modeling (e.g. to base inference on surveys on or after day x)
  
  PARAMETER_VECTOR(log_N);       //log total abundance for each species
  PARAMETER_VECTOR(Beta);              // spline effects on density
  PARAMETER_VECTOR( thin_beta_day);     //extra day and day^2 effects
  PARAMETER_VECTOR(phi_log_rus);  //tweedie phi
  PARAMETER_VECTOR(p_logit_rus);  // tweedie power
  PARAMETER_VECTOR(phi_log_us);  //tweedie phi
  PARAMETER_VECTOR(p_logit_us);  // tweedie power
  PARAMETER_VECTOR( thin_logit_us_i );         // thinning "parameter" for each surveyed location (assumed MVN on logit scale) - bearded & ringed are stacked
  PARAMETER_VECTOR( thin_logit_rus_i );         // thinning "parameter" for each surveyed location (assumed MVN on logit scale)
  //PARAMETER(beta0_ringed);
  //PARAMETER(beta1_ringed);   //logistic params for ringed seal availability relative to maximum
  //PARAMETER(beta2_ringed);   //logistic params for ringed seal availability relative to maximum
  PARAMETER_VECTOR(Beta_melt); //ringed spline params on availability
  PARAMETER(p_sp_logit); //probability of getting a species ID for US data (logit scale)
  PARAMETER_VECTOR(log_lambda);//Penalization parameters
  PARAMETER(log_lambda_avail); 
  PARAMETER(intercept_avail);
  
  // derived sizes
  int n_rus_i = C_rus.col(0).size();
  int n_us_i = C_us.col(0).size();
  int n_obs_types_rus = C_rus.row(0).size();
  int n_obs_types_us = C_us.row(0).size();
  int n_st = X_s.col(0).size()/n_sp;
  int n_b = X_s.row(0).size()/n_sp;
  int n_lambda = log_lambda.size();

  // global stuff
  vector<Type> N(n_sp);
  for(int isp=0;isp<n_sp;isp++)N(isp)=exp(log_N(isp));
  MVNORM_t<Type>neg_log_density_thin_rus(Thin_sigma_rus);
  MVNORM_t<Type>neg_log_density_thin_us(Thin_sigma_us);
  
  vector<Type> Z_mean(n_s);

  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  Type cur_sum;
  vector<Type> phi_us(n_sp+1);
  vector<Type> power_us(n_sp+1);
  vector<Type> phi_rus(n_sp+1);
  vector<Type> power_rus(n_sp+1);
  for(int i=0;i<(n_sp+1);i++){
    phi_us(i)= exp(phi_log_us(i));
    power_us(i)= 1.0+1.0/(1.0+exp(-p_logit_us(i)));
    phi_rus(i)= exp(phi_log_rus(i));
    power_rus(i)= 1.0+1.0/(1.0+exp(-p_logit_rus(i)));
  }
  vector<Type> lambda=exp(log_lambda);
  Type lambda_avail = exp(log_lambda_avail);

  //set up thinning matrix
  matrix<Type> Thin_rus(n_sp,n_rus_i);
  matrix<Type> Thin_rus_trans(n_sp,n_rus_i);
  matrix<Type> Thin_us(n_sp,n_us_i);
  matrix<Type> Thin_us_trans(n_sp,n_us_i);

  vector<Type> Day_effect_rus(n_rus_i);
  vector<Type> Day_effect_us(n_us_i);
  vector<Type> Day_effect_rus_big(2*n_rus_i); //for prior log likelihood (bearded + ringed)
  vector<Type> Day_effect_us_big(2*n_us_i);


  vector<Type> h_mean_obs(2);  //just for bearded seals (US, Russia)

  Day_effect_rus = X_day_rus * thin_beta_day;
  Day_effect_us = X_day_us * thin_beta_day;

  h_mean_obs = h_mean_obs.setZero();
  Day_effect_rus_big = Day_effect_rus_big.setZero();
  Day_effect_us_big = Day_effect_us_big.setZero();

  for( int i=0;i<n_rus_i;i++){
    Day_effect_rus_big(i)=Day_effect_rus(i);
    Thin_rus_trans(0,i)=1.0/(1.0+exp(-thin_logit_rus_i(i)-Day_effect_rus(i)));
    Thin_rus(0,i)=P_rus_i(i)/A_s(S_rus_i(i))*Thin_rus_trans(0,i);
    h_mean_obs(0) += Thin_rus_trans(0,i);
  }
  h_mean_obs(0) = h_mean_obs(0) / n_rus_i;

  for( int i=0;i<n_us_i;i++){
    Day_effect_us_big(i)=Day_effect_us(i);
    Thin_us_trans(0,i)=1.0/(1.0+exp(-thin_logit_us_i(i)-Day_effect_us(i)));
    Thin_us(0,i)=P_us_i(i)/A_s(S_us_i(i))*Thin_us_trans(0,i);
    h_mean_obs(1) += Thin_us_trans(0,i);
  }
  h_mean_obs(1) = h_mean_obs(1) / n_us_i;

  // spline-based ringed seal availability model
  vector<Type> Lin_pred_us_i = X_ringed_us_i * Beta_melt;
  vector<Type> Lin_pred_rus_i = X_ringed_rus_i * Beta_melt;
  for( int i=0;i<n_us_i;i++){
    if(ho_fix==1)Thin_us_trans(1,i)=1.0/(1.0+exp(-thin_logit_us_i(n_us_i+i)));  //changed to only allow setting Pr(lair)=0
    else{
      Thin_us_trans(1,i)= 1.0/((1.0+exp(-(intercept_avail+Lin_pred_us_i(i))))*(1.0+exp(-thin_logit_us_i(n_us_i+i))));  //detection * avail
    }
    Thin_us(1,i)=P_us_i(i)/A_s(S_us_i(i))*Thin_us_trans(1,i);
  }
  for( int i=0;i<n_rus_i;i++){
    if(ho_fix==1)Thin_rus_trans(1,i)=1.0/(1.0+exp(-thin_logit_rus_i(n_rus_i+i)));
    else{
      Thin_rus_trans(1,i)= 1.0/((1.0+exp(-(intercept_avail+Lin_pred_rus_i(i))))*(1.0+exp(-thin_logit_rus_i(n_rus_i+i))));  //detection * avail
    }
    Thin_rus(1,i)=P_rus_i(i)/A_s(S_rus_i(i))*Thin_rus_trans(1,i);
  }

  Type p_species = 1.0/(1.0+exp(-p_sp_logit));

  //spline prior
  int k=0;  // Counter
  for(int i=0;i<n_lambda;i++){
    int m_i = Sdims(i);
    vector<Type> beta_i = Beta.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    jnll_comp(3) -= Type(0.5)*m_i*log_lambda(i) - 0.5*lambda(i)*GMRF(S_i).Quadform(beta_i); //note from Devin: m_i would need to be rank(S) if S_i not full rank (e.g. thin plate)
    k += m_i;
    jnll_comp(3) -= (Type(-0.95)*log_lambda(i)-0.005*exp(log_lambda(i)));  //gamma(0.05,0.005) prior like in Jagam
  }
  //spline availability prior
  jnll_comp(2) -= Type(0.5)*Sdim_melt*log_lambda_avail - 0.5*lambda_avail*GMRF(S_melt).Quadform(Beta_melt); //note from Devin: m_i would need to be rank(S) if S_i not full rank (e.g. thin plate)
  jnll_comp(2) -= (Type(-0.95)*log_lambda_avail-0.005*exp(log_lambda_avail));  //gamma(0.05,0.005) prior like in Jagam


  // Predicted densities
  matrix<Type> Pi_s(n_sp,n_st);
  matrix<Type> Z_s(n_sp,n_st);
  matrix<Type> E_count_sp_rus_i(n_sp,n_rus_i);
  matrix<Type> E_count_sp_us_i(n_sp,n_us_i);

  matrix<Type> E_count_obs_us_i(n_us_i,n_obs_types_us);
  matrix<Type> E_count_obs_rus_i(n_rus_i,n_obs_types_rus);

  //expected abundance
  vector<Type> linpredZ_s=X_s*Beta+I_no_ice*Type(-20.0);
  matrix<Type> Z_pred_0(n_sp,n_s);

  for(int isp=0;isp<n_sp;isp++){
     for(int ist=0; ist<n_st; ist++){
       Pi_s(isp,ist) = A_s(ist)*exp( linpredZ_s(isp*n_st+ist) );
     }
     for(int it=0;it<n_t;it++){  //potentially set some to zero
        for(int is=0;is<n_s;is++){
            Pi_s(isp,it*n_s+is)=Pi_s(isp,it*n_s+is)*I_cell(is);
        }
     }
    for(int it=0;it<n_t;it++){
      cur_sum = 0;
      for(int is=0;is<n_s;is++){
        cur_sum = cur_sum + Pi_s(isp,it*n_s+is);
      }
      for(int is=0;is<n_s;is++){
        Pi_s(isp,it*n_s+is)=Pi_s(isp,it*n_s+is)/cur_sum;
        Z_s(isp,it*n_s+is)=N(isp)*Pi_s(isp,it*n_s+is);
      }
    }
    for(int is=0;is<n_s;is++){
      Z_pred_0(isp,is) = Z_s(isp,is);  //for getting standard errors
    }
  }
  
  int isp=0;
  for(int is=0;is<n_s;is++){
    Z_mean(is)=0;
    for(int it=0;it<n_t;it++){
      Z_mean(is)+=Z_s(isp,it*n_s+is);
    }
    Z_mean(is)=Z_mean(is)/n_t;
  }

  //estimates for bearded seals omitting Chaunskaya Bay --- not currently used for paper
  vector<Type> N_bearded(n_t);
  N_bearded.setZero();
  Type n_s_b = Which_cells_bearded.size();
  cur_sum=0;
  for(int it=0;it<n_t;it++){
    for(int is=0;is<n_s_b;is++){
      N_bearded(it)+=Z_s(0,it*n_s+Which_cells_bearded(is));
    }
  }


  // Probability of counts
  E_count_obs_us_i = E_count_obs_us_i.setZero();
  E_count_obs_rus_i = E_count_obs_rus_i.setZero();
  for(int i=0; i<n_us_i; i++){
    for(int isp=0;isp<n_sp;isp++){
      E_count_sp_us_i(isp,i)=Z_s(isp,S_us_i(i))*Thin_us(isp,i);
    }
    E_count_obs_us_i(i,0)=E_count_sp_us_i(0,i)*p_species;
    E_count_obs_us_i(i,1)=E_count_sp_us_i(1,i)*p_species;
    E_count_obs_us_i(i,2)=(E_count_sp_us_i(0,i)+E_count_sp_us_i(1,i))*(1.0-p_species);
  }
  for(int i=0; i<n_rus_i; i++){
    for(int isp=0;isp<n_sp;isp++){
      E_count_sp_rus_i(isp,i)=Z_s(isp,S_rus_i(i))*Thin_rus(isp,i);
    }
    E_count_obs_rus_i(i,0)=E_count_sp_rus_i(0,i)*Prop_photo_i(i);
    E_count_obs_rus_i(i,1)=E_count_sp_rus_i(1,i)*Prop_photo_i(i);
    E_count_obs_rus_i(i,2)=(E_count_sp_rus_i(0,i)+E_count_sp_rus_i(1,i))*(1.0-Prop_photo_i(i));
  }
  for(int i=0;i<n_us_i;i++){
    for(int itype = 0;itype<n_obs_types_us;itype++){
        if(itype==0 | (Day_us_i(i)>=day_start & itype>0)){
          jnll_comp(0) -= dtweedie( C_us(i,itype), E_count_obs_us_i(i,itype),phi_us(itype),power_us(itype), true);
      }
    }
  }
  for(int i=0;i<n_rus_i;i++){
    for(int itype = 0;itype<n_obs_types_rus;itype++){
      if(itype==0 | (Day_rus_i(i)>=day_start & itype>0)){
        jnll_comp(0) -= dtweedie( C_rus(i,itype), E_count_obs_rus_i(i,itype),phi_rus(itype),power_rus(itype), true );
      }
    }
  }

  // // Probability of thinning and misID parameters (MVN prior/penalty)

  //jnll_comp(2) = 0.0;
  //jnll_comp(1) = neg_log_density_thin_rus(thin_logit_rus_i-Thin_mu_rus_i-Day_effect_rus_big);
  //jnll_comp(1) += neg_log_density_thin_us(thin_logit_us_i-Thin_mu_us_i-Day_effect_us_big);
  jnll_comp(1) = neg_log_density_thin_rus(thin_logit_rus_i-Thin_mu_rus_i);
  jnll_comp(1) += neg_log_density_thin_us(thin_logit_us_i-Thin_mu_us_i);

  for(int ispeff = 0; ispeff<2;ispeff++) jnll_comp(1) += pow(thin_beta_day(ispeff),2.0);  //ridge reg
  jnll_comp(1) += wt_ho * pow(h_mean_obs(1) - h_mean_us, 2.0);
  jnll_comp(1) += wt_ho * pow(h_mean_obs(0) - h_mean_rus, 2.0);


  //penalty on availability differing from ho_match
  int n_ho_match = ho_match.size();
  vector<Type> Lin_pred_match = intercept_avail + X_ringed_match * Beta_melt;
  for(int imelt = 0; imelt<n_ho_match; imelt++){
    jnll_comp(2) += wt_ringed * pow(ho_match(imelt) - 1.0 / (1.0+exp(-Lin_pred_match(ho_melt(imelt)-1))),2.0);
  }

  //std::cout<<"thin dens "<<jnll_comp(1)<<'\n';

  
  // Total objective
  Type jnll = jnll_comp.sum();
  //Type jnll = 0.0;
  
  //std::cout << jnll_comp << "\n";
  //std::cout<<"Range "<<Range_eta<<"\n";
  REPORT( N );
  REPORT( Z_s );
  //REPORT( total_abundance );
  REPORT( Beta );
  REPORT( Beta_melt);
  REPORT( thin_beta_day);
  REPORT( Thin_us_trans);
  REPORT( Thin_rus_trans);
  REPORT( thin_logit_rus_i);
  REPORT( thin_logit_us_i);
  REPORT( Day_effect_us );
  REPORT( Day_effect_rus );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( E_count_sp_us_i);
  REPORT( E_count_sp_rus_i);
  REPORT( E_count_obs_us_i);
  REPORT( E_count_obs_rus_i);
  REPORT( phi_rus );
  REPORT( power_rus);
  REPORT( phi_us );
  REPORT( power_us);
  REPORT(h_mean_obs);
  //REPORT(N_bearded);
  REPORT(log_lambda);
  REPORT(Z_mean);

  // Bias correction output
  vector<Type> splineForReport = designMatrixForReport*Beta;
  vector<Type> splineForReport_melt = 1.0/(1.0+exp(-(intercept_avail+designMatrixForReport_melt*Beta_melt)));
  vector<Type> avail_pred_wBS = 1.0/(1.0+exp(-(intercept_avail+designMatrix_melt_wBS*Beta_melt)));
  vector<Type> avail_pred_eBS_2012 = 1.0/(1.0+exp(-(intercept_avail+designMatrix_melt_eBS_2012*Beta_melt)));
  vector<Type> avail_pred_eBS_2013 = 1.0/(1.0+exp(-(intercept_avail+designMatrix_melt_eBS_2013*Beta_melt)));
  ADREPORT(Z_mean);

  //ADREPORT(Beta);
  //if(Options_vec(0)==1){
  //  ADREPORT( Beta );
  //ADREPORT(Z_pred_0);
  
  return jnll;
}
