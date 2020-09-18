# R code for "Optimal call center forecasting and staffing" by Sihan Ding and Ger Koole
# Ger Koole, September 2020

##################################################################
# example introduction
# parameters
mean_rate=100;sd=20;n=10^5;beta=4;sl=0.8;tta=1/3
c=10;cu=5;cl=1
# samples of arrival rate
lambda=pmax(rnorm(n,mean,sd),0)
# required staffing
S_lambda=vector()
for(k in 1:n)S_lambda[k]=ErlangC_Agents(lambda[k],beta,sl,tta)
# WAPE
mean(abs(lambda-mean)/lambda)
# full info
c*mean(S_lambda)
# s for expected rate
s=ErlangC_Agents(mean_rate,beta,sl,tta)
s
mean(c*S_lambda+cu*pmax(S_lambda-s,0)-cl*pmin(S_lambda-s,0))
# optimal s
staffing_rate=qnorm(cu/(cu+cl),mean_rate,sd)
staffing_rate
s=ErlangC_Agents(staffing_rate,beta,sl,tta)
s
mean(c*S_lambda+cu*pmax(S_lambda-s,0)-cl*pmin(S_lambda-s,0))

##################################################################
# figure on multiplicative vs additive forecasting model
# data obtained from Giuseppe Catanese from anonymous customer of CCmath, April 2020
actuals=c(1275,927,1000,937,1201,3240,2337,2448,2405,2978)
weeks=c(rep(1,5),rep(2,5))
days=rep(c("Mo","Tu","We","Th","Fr"),2)
add_pred=predict(lm(actuals~factor(weeks)+days))
sum(abs(actuals-add_pred))/sum(actuals)
mult_pred=exp(predict(lm(log(actuals+1)~factor(weeks)+days)))
sum(abs(actuals-mult_pred))/sum(actuals)
axis(1,at=1:10,labels=days)
plot(actuals,type="l",xaxt="n")
axis(1,at=1:10,labels=days)
lines(add_pred,lty=2,col="blue")
lines(mult_pred,lty=3, col="red")
legend(1,3000,legend=c("Actuals","Additive fit","Multiplicative fit"),lty=1:3,col=c("black","blue","red"))

##################################################################
# Erlang C formulas
ErlangB_BlockingProbability = function(s,a)
{return(dpois(s,a)/ppois(s,a))}

ErlangC_DelayProbability = function(s,a)
{return(s*ErlangB_BlockingProbability(s, a) / (s - a * (1 - ErlangB_BlockingProbability(s, a))))}

# integer value for s required
ErlangCinteger_ServiceLevel = function(lambda,beta,s,tta)
{
  if(lambda*beta<s){return(1-ErlangC_DelayProbability(s,lambda*beta)*exp(tta*(lambda-s/beta)))
  }else{if(lambda==0){return(1)}else{return(0)}}
}

# s = # agents can be non-integer, return linear interpolation between floor and ceiling(s)
ErlangC_ServiceLevel = function(lambda,beta,s,tta)
{
  return((1+floor(s)-s)*ErlangCinteger_ServiceLevel(lambda,beta,floor(s),tta)+
           (s-floor(s))*ErlangCinteger_ServiceLevel(lambda,beta,floor(s)+1,tta))
}

# returns fractional value by linear interpolation 
ErlangC_Agents = function(lambda,beta,sl,tta)
{
  s=ceiling(lambda*beta)
  while(ErlangC_ServiceLevel(lambda,beta,s,tta)<sl){s=s+1}
  return(s-(ErlangC_ServiceLevel(lambda,beta,s,tta)-sl)/
           (ErlangC_ServiceLevel(lambda,beta,s,tta)-ErlangC_ServiceLevel(lambda,beta,s-1,tta)))
}

