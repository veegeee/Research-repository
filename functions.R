#### A. GIBBS sampler with MH step
#the arguments that must be given are a dataset, indices for x and y variables. if no indices for x variables are given
#model with zero predictors is run
sampled_posterior<-function(dataset, y, x=c(), iterations, initial_values=c(0,0,0,1), burnin=0, seed=1234){
  set.seed(seed) #assures reproducibility of results
  if(length(x)==2){#if two predictors are given this block is run
    beta0_i<-initial_values[1] #initial values are set
    beta1_i<-initial_values[2]
    beta2_i<-initial_values[3]
    var_i<-initial_values[4]
    y_vector<-as.vector(unlist(dataset[,y])) # x and y vectors are set
    x_vector1<-as.vector(unlist(dataset[,x[1]]))
    x_vector2<-as.vector(unlist(dataset[,x[2]]))}
  if(length(x)==1){ ##if one predictor is given this block is run
    beta0_i<-initial_values[1]
    beta1_i<-initial_values[2]
    beta2_i<-0 #beta2_i is set to zero because it is not estimated if one predictor given
    var_i<-initial_values[4]
    y_vector<-as.vector(unlist(dataset[,y]))
    x_vector1<-as.vector(unlist(dataset[,x[1]]))
    x_vector2<-rep(0, nrow(dataset))} #vector of x2 is set to zero if one predictor is given
  if(length(x)==0){
    beta0_i<-initial_values[1]
    beta1_i<-0
    beta2_i<-0
    var_i<-initial_values[4]
    y_vector<-as.vector(unlist(dataset[,y]))
    x_vector1<-rep(0,nrow(dataset))
    x_vector2<-rep(0, nrow(dataset))}
  #a dataframe is created in which thetas will be stored
  sampled_values<-data.frame(beta0=rep(0,iterations), beta1=rep(0,iterations), beta2=rep(0,iterations), var=rep(0,iterations))
  
  for(i in 1:(iterations+burnin)){ #thetas are sampled a desired number of times, if burn in specified then this is added to n iterations
    #the mean for the conditional posterior of the intercept is defined
    mean_beta0<-(sum(y_vector-beta1_i*x_vector1-beta2_i*x_vector2)/var_i + 0/1000)/(nrow(dataset)/var_i+1/1000)
    var_beta0<-1/(nrow(dataset)/var_i+1/1000) #the variance of the conditional posterior of the intercept is defined
    sampled_beta0<-rnorm(1, mean_beta0, sqrt(var_beta0)) # a value for the intercept is sampled
    beta0_i<-sampled_beta0 #the value of the intercept is updated and will be used for the means of the other conditional posteriors
    #in this iteration
    sampled_values[i,1]<-sampled_beta0 #the sampled value is stored in the data
    ##Metropolis-Hastings step
    if(length(x)>0){ #the block is run when there is at least one predictor
      #the mean for the proposal density is defined
      mean_beta1<-(sum(x_vector1*((y_vector-beta0_i)-beta2_i*x_vector2))/var_i + 0/1000)/(sum(x_vector1*x_vector1)/var_i+1/1000)
      var_beta1<-1/(sum(x_vector1*x_vector1)/var_i+1/1000) #the variance of the proposal density is defined
      sampled_beta1<-rnorm(1, beta1_i, sqrt(0.01))
      prior<-dnorm(sampled_beta1,0, sqrt(1000), log = TRUE) #the logarithm of the prior probability of sampled beta1 is calculated
      prior1<-dnorm(beta1_i,0, sqrt(1000), log = TRUE) #the logarithm of the prior probability of previously retained beta1 is calculated
      unif<-log(runif(1,0,1)) #value U drawn from uniform distribution and its logarithm is calculated
      #the loglikelihoods of the sampled and previously retained beta1 are calculated
      loglik<-sum(dnorm(y_vector, (beta0_i+sampled_beta1*x_vector1+beta2_i*x_vector2), sqrt(var_i), log = TRUE))
      loglik1<-sum(dnorm(y_vector, (beta0_i+beta1_i*x_vector1+beta2_i*x_vector2), sqrt(var_i), log = TRUE))
      if(unif<=(loglik+prior-loglik1-prior1)){ #the logarithm of the value U compared to the logarithm of the acceptance ratio
        #if logarithm of U is smaller or equal to logarithm of acceptance ratio the sampled value is retained
        beta1_i<-sampled_beta1
        sampled_values[i,2]<-sampled_beta1
      }else{sampled_values[i,2]<-beta1_i} #if logarithm of U is larger than logarithm of acceptance ratio the previously retained
      #beta1 is retained
    }
    
    if(length(x)==2){ #if two predictors are included this block is executed
      #the code follows the same steps as the code for beta0
      mean_beta2<-(sum(x_vector2*((y_vector-beta0_i)-beta1_i*x_vector1))/var_i + 0/1000)/(sum(x_vector2*x_vector2)/var_i+1/1000)
      var_beta2<-1/(sum(x_vector2*x_vector2)/var_i+1/1000)
      sampled_beta2<-rnorm(1, mean_beta2, sqrt(var_beta2))
      beta2_i<-sampled_beta2
      sampled_values[i,3]<-sampled_beta2}
    
    #the variance is estimated for every number of predictors
    alpha_var<-nrow(dataset)/2+0.001 #the alpha of the conditional posterior inverse gamma distribution of variance is defined
    ##the beta of the conditional posterior inverse gamma distribution of variance is defined
    beta_var<- sum((y_vector-(rep(beta0_i, nrow(dataset))+beta1_i*x_vector1+beta2_i*x_vector2))^2)/2+0.001
    sampled_var<-1/rgamma(1, alpha_var, beta_var)
    var_i<-sampled_var
    sampled_values[i,4]<-sampled_var
    #print(beta1_i)
  }
  sampled_values<-sampled_values[(burnin+1):(iterations+burnin),] #the thetas used for burn in are removed
  if(length(x)==1){ #if 1 predictor the column for the second predictor which contains zeroes is dropped
    sampled_values<-sampled_values[,c(1,2,4)]
  }
  if(length(x)==0){ #if no predictors the columns of predictors are dropped which contain zeroes
    sampled_values<-sampled_values[,c(1,4)]
  }
  return(sampled_values) #sampled thetas are returned
}


#the gelman rubin statistics are calculated for different parts of the sampled posterior
#first it is calculated for first 500 thetas, then first 1000 thetas, then first 1500 thetas etcetera
gelman_rubin<-function(chain1, chain2){
  number_of_statistics<-length(chain1)%/%500 #the length of chain divided by 500 and dropping the decimals to get 
  #number of gelman-rubin statistics that must be calculated
  gelman_rubin_statistics<-rep(0, number_of_statistics) #vector where the statistics are stored
  for(i in 1 : number_of_statistics){
    #i*500 is included to subset dynamically, first the first 500 values are used then first 1000 etcetera
    chain_1i<-chain1[1:(i*500)]
    chain_2i<-chain2[1:(i*500)]
    grand_mean<-mean(c(chain_1i, chain_2i)) #grand mean calculated
    #variance for both chains is calculated
    var_chain1i<-var(chain_1i) 
    var_chain2i<-var(chain_2i)
    within_variance<-(var_chain1i+var_chain2i)/2 #within variance calculated
    mean_chain1i<-mean(chain_1i)
    mean_chain2i<-mean(chain_2i)
    #with 2 chains, the denominator of the first term of between variance equals 1, i*500 equals the number of sampled thetas 
    between_variance<-(i*500)*(((mean_chain1i-grand_mean)^2) + ((mean_chain2i-grand_mean)^2)) #between variance calculated
    total_variance<-(((i*500)-1)/(i*500))*within_variance + (1/(i*500))*between_variance #total variance calculated
    gelman_rubin_statistics[i]<-total_variance/within_variance
  }
  return(gelman_rubin_statistics) #statistics are spit out
}


##bayesian p-value is calculated
p_value<-function(theta, data, x, y, seed= 1234){
  set.seed(seed) #assures reproducibility
  #define x and y vectors
  x1<-data[,x[1]]
  x2<-data[,x[2]]
  y<-unlist(data[,y])
  y<-as.vector(y)
  abs_differences_obs<-c()
  abs_differences_sampled<-c()
  vector_comparisons<-rep(0, nrow(theta))
  for(i in 1: nrow(theta)){ #repeat this block for each sampled theta
    mean_iteration_i<-theta[i,1] + theta[i,2]*x1 + theta[i,3]*x2 #predicted value based on theta i
    sampled_ys<-rnorm(1599, theta[i,1] + theta[i,2]*x1 + theta[i,3]*x2, sd = sqrt(theta[i,4])) #sample ys
    residuals_sampled<-sampled_ys - mean_iteration_i #create residuals for sampled ys
    residuals_observed<-y - mean_iteration_i #create residuals for observed ys
    #test statistic computed for oberved and sampled data. It is equal to absolute value of difference mean and median of residuals
    abs_difference_median_mean_sampled<-abs(mean(residuals_sampled)-median(residuals_sampled))
    abs_difference_median_mean_observed<-abs(mean(residuals_observed)-median(residuals_observed))
    #save the differences
    abs_differences_obs[i]<-abs_difference_median_mean_observed
    abs_differences_sampled[i]<-abs_difference_median_mean_sampled
    # 1 is stored in vector if difference bigger for sampled values
    vector_comparisons[i]<-abs_difference_median_mean_sampled > abs_difference_median_mean_observed
    
  }
  #p-value equals mean(vector_comparisons)
  return(list(p_value=mean(vector_comparisons), sampled_differences=abs_differences_sampled, observed_differences=abs_differences_obs))
}


#standard deviations calculated
standard_deviation<-function(theta){
  mean_theta<-mean(theta)
  squared_differences<-(theta-mean_theta)^2
  #to get standard deviation of posterior, divide the sum of squares by the number of sampled thetas, not by: number of sampled thetas - 1.
  standard_deviation<-sqrt(sum(squared_differences)/length(theta))
  return(standard_deviation)
}


#three functions are created. One for each number of predictors
#works for 2 predictors
dic1<-function(thetas, data, y,x){
  #x and y vectors are defined
  y_vector<-as.vector(unlist(data[,y]))
  x_vector1<-as.vector(unlist(data[,x[1]]))
  x_vector2<-as.vector(unlist(data[,x[2]]))
  mean_theta<-colMeans(thetas) #means calculated for the loglikelihood of the mean
  #loglikelihood of mean calculated
  loglik_mean_thetas<-sum(dnorm(y_vector, mean_theta[1]+mean_theta[2]*x_vector1+mean_theta[3]*x_vector2, sqrt(mean_theta[4]), log = TRUE))
  
  loglik_vector_thetas<-c() #loglikelihood for each theta is stored here
  
  for(i in 1:nrow(thetas)){ #for each theta repeated
    #loglikelihood for the particular theta
    loglik_vector_thetas[i]<-sum(dnorm(y_vector, thetas[i,1]+thetas[i,2]*x_vector1+thetas[i,3]*x_vector2, sqrt(thetas[i,4]), log = TRUE))
  }
  mean_loglik_thetas<-mean(loglik_vector_thetas) #the mean of the loglikelihoods of all thetas calculated
  dic<- -2*loglik_mean_thetas + 2*(-2*mean_loglik_thetas + 2*loglik_mean_thetas) #formula dic applied
  return(dic) #dic returned
}

dic2<-function(thetas, data, y,x){  #same logic as first dic function, but 1 predictors
  y_vector<-as.vector(unlist(data[,y]))
  x_vector1<-as.vector(unlist(data[,x]))
  mean_theta<-colMeans(thetas)
  loglik_mean_thetas<-sum(dnorm(y_vector, mean_theta[1]+mean_theta[2]*x_vector1, sqrt(mean_theta[3]), log = TRUE))
  
  loglik_vector_thetas<-c()
  
  for(i in 1:nrow(thetas)){
    loglik_vector_thetas[i]<-sum(dnorm(y_vector, thetas[i,1]+thetas[i,2]*x_vector1, sqrt(thetas[i,3]), log = TRUE))
  }
  mean_loglik_thetas<-mean(loglik_vector_thetas)
  dic<- -2*loglik_mean_thetas + 2*(-2*mean_loglik_thetas + 2*loglik_mean_thetas)
  return(dic)
}

#works for no predictors
dic3<-function(thetas, data, y){ #same logic as first dic function, but no predictors
  y_vector<-as.vector(unlist(data[,y]))
  mean_theta<-colMeans(thetas)
  loglik_mean_thetas<-sum(dnorm(y_vector, mean_theta[1], sqrt(mean_theta[2]), log = TRUE))
  
  loglik_vector_thetas<-c()
  
  for(i in 1:nrow(thetas)){
    loglik_vector_thetas[i]<-sum(dnorm(y_vector, thetas[i,1], sqrt(thetas[i,2]), log = TRUE))
  }
  mean_loglik_thetas<-mean(loglik_vector_thetas)
  dic<- -2*loglik_mean_thetas + 2*(-2*mean_loglik_thetas + 2*loglik_mean_thetas)
  return(dic)
}