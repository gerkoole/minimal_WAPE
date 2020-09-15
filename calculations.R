# R code for "Optimal call center forecasting and staffing" by Sihan Ding and Ger Koole
# Ger Koole, September 2020

##################################################################
# example introduction
# parameters
lambda=100;sd=20;n=10^7;beta=4;alpha=0.5
c=10;cu=5;cl=1
# Erlang C SL with value of alpha
ErlangC_ServiceLevel(lambda,beta,lambda*beta+alpha*sqrt(lambda*beta),1/3)
# simulations
rates=pmax(rnorm(n,lambda,sd),0)
load=rates*beta
# WAPE
mean(abs(rates-lambda))/lambda
# full info
c*mean(load+alpha*sqrt(load))
# s for expected rate
s=lambda*beta+alpha*sqrt(lambda*beta)
s
mean(c*(load+alpha*sqrt(load))+cu*pmax(load+alpha*sqrt(load)-s,0)-cl*pmin(load+alpha*sqrt(load)-s,0))
# optimal s
qnorm(cu/(cu+cl),lambda*beta,sd)
s=qnorm(cu/(cu+cl),lambda*beta,sd)+alpha*sqrt(qnorm(cu/(cu+cl),lambda*beta,sd))
s
mean(c*(load+alpha*sqrt(load))+cu*pmax(load+alpha*sqrt(load)-s,0)-cl*pmin(load+alpha*sqrt(load)-s,0))

##################################################################
# Erlang C formulas
ErlangB_BlockingProbability = function(s,a)
{return(dpois(s,a)/ppois(s,a))}

ErlangC_DelayProbability = function(s,a)
{return(s*ErlangB_BlockingProbability(s, a) / (s - a * (1 - ErlangB_BlockingProbability(s, a))))}

# Erlang C; tta = time to answer (often 20 s = 1/3 min)
ErlangC_ServiceLevel = function(lambda,beta,s,tta)
{
  if(lambda*beta<s){return(1-ErlangC_DelayProbability(s,lambda*beta)*exp(tta*(lambda-s/beta)))
  }else{if(lambda==0){return(1)}else{return(0)}}
}

