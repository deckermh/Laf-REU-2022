# generate a ton of matrices with sigmas that:
# 
# 1.) Have the same variances of variances but different ratios
# of max to min
#
# ex.) sigmas are: 1,2,3 (1/3) then 2,3,4 (1/2) then 3,4,5 (3/5)
# then 4,5,6 (2/3) ...
# 2.) Have the same ratio of max to min but different variance of variances
#
# ex.) 1,1,10 then 1,2,10 then 1,3,10 then 1,4,10 then 1,5,10 then...
# Warning: notice that 1,4,10 has same var of var as 1,6,10, etc

#num trials
N = 5

#rule of thumb
thumb = 3

#generate examples in 1.), varying ratio of min to max
for (i in 1:N) {
  # 1,2,3 then 2,3,4 then 3,4,5 then...
  sigmas = c(sqrt(i), sqrt(i+1), sqrt(i+2))
  
  #Note that the matrix is CSH (we could also do ARH1)
  Sigma = makeCSH(3, .5, sigmas)
  
  res = results_matrix(100, 3, 10, Sigma, c(0,0,0))

  differences = diff(res, 13)
  
  #get success/failure rates
  successRates = quad_correct(differences, thumb)
  
  #writing file names and saving to csv tee hee
  ratio = sigmas[1]^2 / sigmas[3]^2
  filename = paste("quad_matrix_varying_ratios/", ratio, ".csv")
  write.csv(successRates, filename, row.names = FALSE)
}

N = 10
thumb = 3

#generate examples in 2.), varying the variance of variances
for (i in 1:N) {
  # 1,2,4 then 1,2.4,4 then 1,2.8,4 then...
  sigmas = c(sqrt(1), sqrt(2+(i-1)*4/10), sqrt(10))
  
  #Note that the matrix is CSH (we could also do ARH1)
  Sigma = makeCSH(3, .5, sigmas)
  
  res = results_matrix(100, 3, 10, Sigma, c(0,0,0))
  
  differences = diff(res, 13)
  
  #get success/failure rates
  successRates = quad_correct(differences, thumb)
  
  #writing file names and saving to csv tee hee
  sigmasSquared = c(sigmas[1]^2, sigmas[2]^2, sigmas[3]^2)
  varOfVar = var(sigmasSquared)
  filename = paste("quad_matrix_varying_variances/", varOfVar, ".csv")
  write.csv(successRates, filename, row.names = FALSE)
}




n_obs = 3
thumb = 3
N = 10
# Run experiment varying the min/max ratio from 1/10 to 1 in increments of 1/10
for (i in 1:N) {
  # 1,1,10 then 2,2,10 then 3,3,10 then...
  sigmas = c(sqrt(i), sqrt(i), sqrt(10))
  
  #Note that the matrix is CSH (we could also do ARH1)
  Sigma = makeCSH(3, .5, sigmas)
  
  res = results_matrix(100, 3, 10, Sigma, c(0,0,0))
  
  differences = diff(res, 13)
  
  #get success/failure rates
  successRates = quad_correct(differences, thumb)
  
  #writing file names and saving to csv tee hee
  ratio = sigmas[1]^2 / sigmas[3]^2
  filename = paste("quad_matrix_varying_ratios/", ratio, ".csv")
  write.csv(successRates, filename, row.names = FALSE)
}




n_obs = 3
thumb = 3
N = 10
# Run experiment varying the min/max ratio from 1/10 to 1 in increments of 1/10
for (i in 1:N) {
  # 1,10,10 then 2,10,10 then 3,10,10 then... 10,10,10
  sigmas = c(sqrt(i), sqrt(10), sqrt(10))
  
  #Note that the matrix is CSH (we could also do ARH1)
  Sigma = makeCSH(3, .5, sigmas)
  
  res = results_matrix(100, 3, 10, Sigma, c(0,0,0))
  
  differences = diff(res, 13)
  
  #get success/failure rates
  successRates = quad_correct(differences, thumb)
  
  #writing file names and saving to csv tee hee
  ratio = sigmas[1]^2 / sigmas[3]^2
  filename = paste("quad_matrix_varying_ratios/", ratio, ".csv")
  write.csv(successRates, filename, row.names = FALSE)
}










#make a plot of box plots for distribution of Type 4 success percentage in
#in AIC_CS as heteroskedasticity increases
AICHeteroBoxplots <- function() {
  AIC_CS_Data = matrix(nrow = 100, ncol = 10)
  
  #we vary min/max ratio from .1 to .8
  for (i in 1:8) {
    #generate CSH matrix with varying levels of heterosked.
    sigmas = c(sqrt(i/10), sqrt(10), sqrt(10))
    Sigma = makeCSH(3, .5, sigmas)
    
    #get quad correct for data
    res = results_matrix(100, 3, 40, Sigma, c(0, 0, 0))
    d = diff(res, 13)
    quad = quad_correct(d, 3)
    
    #get type 1, 2, 3, and 4 success/error percentages for CS fit
    quad = quad[]
    
  }
  
  #create boxplots
  boxplot(AIC_UN_Data[,1:10], main = "AIC_UN - AIC_CS")
}
