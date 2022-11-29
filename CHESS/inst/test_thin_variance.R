

n_boot = 10000
p_b = rep(0,n_boot)
p_a = p_ba = p_ba2 = matrix(0,2,n_boot)
mu_b = 0.8
sigma_b = 0.05
Mu_a = c(0.5,0.6)
Mu_ba = mu_b*Mu_a
Sigma_a = matrix(c(0.01,0.005,0.005,0.01),2,2)

for(i in 1:n_boot){

  p_a[,i]=mvtnorm::rmvnorm(1,Mu_a,Sigma_a)
  p_b[i] = rnorm(1,mu_b,sigma_b)
  p_ba[,i]=p_b[i]*p_a[,i]
  
  Tmp=matrix(0,2,2)
  for(j in 1:2){
    for(k in 1:2){
       Tmp[j,k]=Mu_a[j]*Mu_a[k]
    }
  }
  Sigma_ab = sigma_b^2*Tmp+(sigma_b^2+mu_b^2)*Sigma_a
  p_ba2[,i] = mvtnorm::rmvnorm(1,Mu_ba,Sigma_ab)
    
}