## Working code of the multivariate bias correction method (parametric bootstrapping)
# Author: Akira Endo (https://github.com/akira-endo)


# mo_pb: Function for multiple overimputation (MO)
# Input:
# Z = test results with potential misclassification (vector)
# p_pos = predicted probability of positive test results (vector)
# sens/spec = test sensitivity/specificity (single value; or vector if vary between individuals)
# n_iter = number of MO iterations. 100 or more recomended; more iterations may be required if the sample size is small ( < 1000)
mo_pb<-function(Z,p_pos,sens,spec,n_iter=100){
  p_pos=pmin(1,pmax(0,p_pos))
  sens=pmin(1,pmax(0,sens))
  spec=pmin(1,pmax(0,spec))
  if(!is.finite(sum(p_pos+sens+spec)))stop("invalid input")
  
  ind_pos<-as.logical(Z)
  if(!is.finite(sum(ind_pos)))stop("invalid input")
  
  ## Prob of flipping the result
  p_imp<-numeric(length(Z)) # prob for bernouli imputation

  pflip_pos<-(sens*(1-p_pos)/p_pos-(1-sens)) * ((1-spec)/(sens+spec-1))
  pflip_pos=pflip_pos[ind_pos]
  pflip_pos=pmin(1,pmax(0,pflip_pos))
  p_imp[ind_pos]=1-pflip_pos
  
  pflip_neg<-(spec*p_pos/(1-p_pos)-(1-spec)) * ((1-sens)/(sens+spec-1))
  pflip_neg=pflip_neg[!ind_pos]
  pflip_neg=pmin(1,pmax(0,pflip_neg))
  p_imp[!ind_pos]=pflip_neg
  
  size<-length(Z)
  z_imp<-sapply(1:n_iter,function(x){
    z<-(runif(size)<p_imp)
    return(z)
  })
  colnames(z_imp)=paste("iter",1:n_iter)
  return(list(z_imp=z_imp,p_imp=p_imp,Z=Z,p_pos=p_pos,sens=sens,spec=spec,n_iter=n_iter))
}


## Example ##
# Prepare sample dataset
set.seed(2019)
datasize=500
sens=0.8
spec=0.9
dataset<-data.frame(test=0,disease=0,age=sample(60,datasize,T),sex=sample(2,datasize,T)-1,vaccine=sample(2,datasize,T,prob = c(0.3,0.6))-1)
dataset$disease=rbinom(datasize,1,(0.3+0.2*(dataset$age<15))*0.4^(dataset$vaccine))
dataset$test[dataset$disease==0]=rbinom(datasize-sum(dataset$disease),1,1-spec)
dataset$test[dataset$disease==1]=rbinom(sum(dataset$disease),1,sens)

head(dataset)

# Predict probability of test positivity
fit0<-glm(test~age+sex+vaccine,family = binomial(link="logit"),dataset)
p_pos<-predict(fit0,type="response")

# Generate MO dataset
n_iter<-100
mo_data<-mo_pb(dataset$test,p_pos,sens,spec,n_iter)

# Fit GLM
mo_result<-matrix(0,8,n_iter)
for(i in 1:n_iter){
  temp_data=cbind(z_imp=mo_data$z_imp[,i],dataset)
  fit<-glm(z_imp~age+sex+vaccine,family = binomial(link="logit"),temp_data)
  mo_result[,i]=summary(fit)$coefficients[,1:2]
}

# Calculate mean and sd by Rubin's rules
mean_coef<-rowMeans(mo_result[1:4,])
names(mean_coef)=names(fit0$coefficients)
# mean of MO estimates
sd_coef<-sqrt(rowMeans(mo_result[5:8,]^2)+apply(mo_result[1:4,],1,var)*(1+1/n_iter))

# Estimates from true disease status
fit_true<-glm(disease~age+sex+vaccine,family = binomial(link="logit"),dataset)

# Compare results
summary(fit_true)$coefficients[,1:2]
summary(fit0)$coefficients[,1:2]
cbind(Estimate=mean_coef,Std.Error=sd_coef)
