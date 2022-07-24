####~~~~~~~~~~~~~~#Covariance Matrix Generators#~~~~~~~~~~~~~~~~~~####
#Sigma Generation Functions####

##
#Make a CSH matrix
##
makeCSH <- function(n, p, sigma) {
  if (n != length(sigma)) {
    return(NULL)
  }
  
  data = matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        data[i, j] = sigma[i]*sigma[j]
      }
      else {
        data[i, j] = p*sigma[i]*sigma[j]
      }
    }
  }

  return(data)
}

##
#Make a CS Matrix
##
makeCS <- function(n, p, sigma) {
  anti_identity = matrix(1, n, n) - diag(c(rep(1, n)))   ###prep for adding psigma terms
  
  CSmatrix = matrix(0, n, n)                             ###create empty matrix
  diag(CSmatrix) = c(rep(sigma ^ 2, n))                  ###input sigma^2 on diag
  
  CSmatrix = CSmatrix + p * (sigma ^ 2) * anti_identity      ###add in psig^2 terms
  return(CSmatrix)
}

##
#Make AR(1)
##
makeAR1 <- function(n, p, sigma) {
  AR1matrix = matrix(0, n, n)                             #empty matrix
  
  for (i in 1:n) {
    for (j in 1:n) {
      AR1matrix[i, j] = (p ^ (abs(i - j))) * (sigma ^ 2)
    }
  }
  
  return(AR1matrix)
}

##
#Make an ARH1 Matrix
##
makeARH1 <-function(n, p, sigma) {
    ####sigma looks like c(x,y,...,z)
    if (length(sigma) == n) {
      ARH1matrix = matrix(0, n, n)                             ###create empty matrix
      for (i in 1:n) {
        for (j in 1:n) {
          if (i == j) {
            ARH1matrix[i, j] = sigma[i] ^ 2
          }
          else
            ARH1matrix[i, j] = (p ^ (abs(i - j))) * sigma[i] * sigma[j]
        }
      }
      return(ARH1matrix)
    }
    else
      return("error, check sigma length")
  }

##
#Make an existing matrix symmetric
##
makeSymm <- function(matr) {
  newMatrix = matr
  #for each row,
  for (i in 1:nrow(matr)) {
    for (j in i:nrow(matr)) {
      newMatrix[j, i] = newMatrix[i, j]
    }
  }
  return(newMatrix)
}

##
#Unstructured Sigmas
##

##obs 3##
##MILD##
#green
matr = matrix(c(1, .5, .7,
                0, 1, .55,
                0, 0, 1), ncol = 3, nrow = 3, byrow = TRUE)
mild_green3 = makeSymm(matr)
#red
matr = matrix(c(1, .5*sqrt(5), .7*sqrt(10),
                0, 5, .55*sqrt(50),
                0, 0, 10), ncol = 3, nrow = 3, byrow = TRUE)
mild_red3 = makeSymm(matr)

##MEDIUM##
#green
matr = matrix(c(1, .5, .7,
                0, 1, .55,
                0, 0, 1), ncol = 3, nrow = 3, byrow = TRUE)
med_green3 = makeSymm(matr)
#red
matr = matrix(c(1, .5*sqrt(5), .7*sqrt(10),
                0, 5, .55*sqrt(50),
                0, 0, 10), ncol = 3, nrow = 3, byrow = TRUE)
med_red3 = makeSymm(matr)

##SPICY##
#green
matr = matrix(c(1, .9, 0.09,
                0, 1, .5,
                0, 0, 1), ncol = 3, nrow = 3, byrow = TRUE)
spicy_green3 = makeSymm(matr)


#red
matr = matrix(c(1, .9*sqrt(5), .23,
                0, 5, .5*sqrt(50),
                0, 0, 10), ncol = 3, nrow = 3, byrow = TRUE)
spicy_red3 = makeSymm(matr)

##obs 5##

##MILD##

#green
matr = matrix(c(1, .5, .3, .2, .2,
                0, 1, .4, .2, .1,
                0, 0, 1, .5, .1,
                0, 0, 0, 1, .3,
                0, 0, 0, 0, 1), ncol = 5, nrow = 5, byrow = TRUE)
mild_green5 = makeSymm(matr)
#red
matr = matrix(c(1, .5*sqrt(2), .3*sqrt(3), .2*sqrt(4), .2*sqrt(5),
                0, 2, .4*sqrt(6), .2*sqrt(8), .1*sqrt(10),
                0, 0, 3, .5*sqrt(12), .1*sqrt(15),
                0, 0, 0, 4, .3*sqrt(20),
                0, 0, 0, 0, 5), ncol = 5, nrow = 5, byrow = TRUE)
mild_red5 = makeSymm(matr)
##MEDIUM##

#green
matr = matrix(c(1, .3, .9, .5, .4,
                0, 1, .3, .4, .2,
                0, 0, 1, .4, .3,
                0, 0, 0, 1, .1,
                0, 0, 0, 0, 1), ncol = 5, nrow = 5, byrow = TRUE)
med_green5 = makeSymm(matr)
#red
matr = matrix(c(1, .3*sqrt(2), .9*sqrt(3), .5*sqrt(4), .4*sqrt(5),
                0, 2, .3*sqrt(6), .4*sqrt(8), .2*sqrt(10),
                0, 0, 3, .4*sqrt(12), .3*sqrt(15),
                0, 0, 0, 4, .1*sqrt(20),
                0, 0, 0, 0, 5), ncol = 5, nrow = 5, byrow = TRUE)
med_red5 = makeSymm(matr)
##SPICY##

#green
matr = matrix(c(1, (-.32), .60, .2, .1,
                0, 1, .4, (-.2), .3,
                0, 0, 1, .2, .68,
                0, 0, 0, 1, (-.18),
                0, 0, 0, 0, 1), ncol = 5, nrow = 5, byrow = TRUE)
spicy_green5 = makeSymm(matr)

#red
matr = matrix(c(1, (-.32)*sqrt(2), .6*sqrt(3), .2*sqrt(4), 1,
                0, 2, .4*sqrt(6), (-.2)*sqrt(8), .3*sqrt(10),
                0, 0, 3, .2*sqrt(12), .68*sqrt(15),
                0, 0, 0, 4, (-.18)*sqrt(20),
                0, 0, 0, 0, 5), ncol = 5, nrow = 5, byrow = TRUE)
spicy_red5 = makeSymm(matr)

##obs 10#

##MILD##
#green
matr = matrix(c(1, .5, .3, .2, .3, .3, .2, .1, .1, .1, 
                0, 1, .4, .2, .2, .3, .2, .2, .1, .1,
                0, 0, 1, .5, .3, .3, .3, .2, .1, .1,
                0, 0, 0, 1, .4, .2, .3, .2, .2, .1,
                0, 0, 0, 0, 1, .4, .3, .3, .2, .2,
                0, 0, 0, 0, 0, 1, .3, .3, .2, .2,
                0, 0, 0, 0, 0, 0, 1, .3, .3, .1,
                0, 0, 0, 0, 0, 0, 0, 1, .3, .3,
                0, 0, 0, 0, 0, 0, 0, 0, 1, .3,
                0, 0, 0, 0, 0, 0, 0, 0 , 0, 1),  ncol = 10, nrow =10, byrow = TRUE)
mild_green10 = makeSymm(matr)
#red
matr = matrix(c(1, .5*sqrt(2), .3*sqrt(3), .2*sqrt(4), .3*sqrt(5), .3*sqrt(6), .2*sqrt(7), .1*sqrt(8), .1*sqrt(9), .1*sqrt(10), 
                0, 2, .4*sqrt(6), .2*sqrt(8), .2*sqrt(10), .3*sqrt(12), .2*sqrt(14), .2*sqrt(16), .1*sqrt(18), .1*sqrt(20),
                0, 0, 3, .5*sqrt(12), .3*sqrt(15), .3*sqrt(18), .3*sqrt(21), .2*sqrt(24), .1*sqrt(27), .1*sqrt(30),
                0, 0, 0, 4, .4*sqrt(20), .2*sqrt(24), .3*sqrt(28), .2*sqrt(32), .2*sqrt(36), .1*sqrt(40),
                0, 0, 0, 0, 5, .4*sqrt(30), .3*sqrt(35), .3*sqrt(40), .2*sqrt(45), .2*sqrt(50),
                0, 0, 0, 0, 0, 6, .3*sqrt(42), .3*sqrt(48), .2*sqrt(54), .2*sqrt(60),
                0, 0, 0, 0, 0, 0, 7, .3*sqrt(56), .3*sqrt(63), .1*sqrt(70),
                0, 0, 0, 0, 0, 0, 0, 8, .3*sqrt(72), .3*sqrt(80),
                0, 0, 0, 0, 0, 0, 0, 0, 9, .3*sqrt(90),
                0, 0, 0, 0, 0, 0, 0, 0 , 0, 10),  ncol = 10, nrow =10, byrow = TRUE)
mild_red10 = makeSymm(matr)
##MEDIUM##

#green
matr = matrix(c(1, .35, .8, .45, .4, .6, .65, .4, .5, .44,
                0, 1, .3, .4, .2, .3, .3, .4, .3, .35,
                0, 0, 1, .4, .3, .2, .25, .4, .2, .35,
                0, 0, 0, 1, .15, .2, .3, .3, .3, .5,
                0, 0, 0, 0, 1, .3, .4, .4, .4, .6,
                0, 0, 0, 0, 0, 1, .6, .35, .3, .4,
                0, 0, 0, 0, 0, 0, 1, .25, .4, .45,
                0, 0, 0, 0, 0, 0, 0, 1, .5, .85,
                0, 0, 0, 0, 0, 0, 0, 0, 1, .65,
                0, 0, 0, 0, 0, 0, 0, 0 , 0, 1),  ncol = 10, nrow =10, byrow = TRUE)
med_green10 = makeSymm(matr)

#red
matr = matrix(c(1, .35*sqrt(2), .8*sqrt(3), .45*sqrt(4), .4*sqrt(5), .6*sqrt(6), .65*sqrt(7), .4*sqrt(8), .5*sqrt(9), .44*sqrt(10),
                0, 2, .3*sqrt(6), .4*sqrt(8), .2*sqrt(10), .3*sqrt(12), .3*sqrt(14), .4*sqrt(16), .3*sqrt(18), .35*sqrt(20),
                0, 0, 3, .4*sqrt(12), .3*sqrt(15), .2*sqrt(18), .25*sqrt(21), .4*sqrt(24), .2*sqrt(27), .35*sqrt(30),
                0, 0, 0, 4, .15*sqrt(20), .2*sqrt(24), .3*sqrt(28), .3*sqrt(32), .3*sqrt(36), .5*sqrt(40),
                0, 0, 0, 0, 5, .3*sqrt(30), .4*sqrt(35), .4*sqrt(40), .4*sqrt(45), .6*sqrt(50),
                0, 0, 0, 0, 0, 6, .6*sqrt(42), .35*sqrt(48), .3*sqrt(54), .4*sqrt(60),
                0, 0, 0, 0, 0, 0, 7, .25*sqrt(56), .4*sqrt(63), .45*sqrt(70),
                0, 0, 0, 0, 0, 0, 0, 8, .5*sqrt(72), .85*sqrt(80),
                0, 0, 0, 0, 0, 0, 0, 0, 9, .65*sqrt(90),
                0, 0, 0, 0, 0, 0, 0, 0 , 0, 10),  ncol = 10, nrow =10, byrow = TRUE)
med_red10 = makeSymm(matr)

####Generates and fits data, varying rho each trial####
rho_experiment <- function(N, n_obs, n_sub, means, variances) {
  #lengths must be consistent
  if (n_obs == length(variances) && length(variances) == length(means)) {
    rho = 0
    #initialize sigma with rho value of 0
    Sigma = makeCSH(n_obs, variances, rho)
    
    #initialize results matrix
    results = matrix(nrow = N, ncol = 15)
    
    for (i in 1:N) {
      #put IC values into results matrix
      data = generate_data(n_obs, n_sub, Sigma, means)
      ICVector = fit_data(data, n_obs, n_sub)
      results[i,] = ICVector
      
      #increment rho by 1/N, so that rho goes from 0 to 1
      #as loop ends
      rho = rho + 1 / N
      Sigma = makeCSH(n_obs, variances, rho)
    }
    
    
    return (results)
  }
  else {
    return (NULL)
  }
}

####~~~~~~~~~~~~~~#Results Matrix Generators#~~~~~~~~~~~~~~~~~~####
####Data Generation Function####
#returns data ready to input into Fit Data function
generate_data <- function(n_obs, n_sub, Sigma, means) {
    ###where Sigma is n_obs x n_sub, means is of form c(x,..,y) w/ n_obs terms
    
    ##making general data form
    datamatrix = matrix(0, n_sub, n_obs + 1)
    colnames(datamatrix) = c("subject", paste("t", 1:(n_obs), sep = ""))
    rownames(datamatrix) = c(1:n_sub)
    datamatrix[1:n_sub, 1] = c(1:n_sub)
    
    ##inputting general data
    library(MASS)
    
    data = mvrnorm(n = n_sub, means, Sigma)
    datamatrix[1:n_sub, 2:(n_obs + 1)] = data
    
    datamatrix = as.data.frame(datamatrix, header = TRUE)  ###converting to reshapeable form
    
    clean_data = reshape(datamatrix, varying = list(2:(n_obs + 1)), direction =
                           "long")
    names(clean_data) = c("subject", "time", "observation", "id")
    return(clean_data)
  }

####Fit Data Function####
#returns results matrix (summary of all IC values for each trial)
fit_data <- function(clean_data, n_obs, n_sub) {
  
  observation = clean_data$observation
  time = clean_data$time
  id = clean_data$id
  
  time = as.factor(time)
  
  
  library(nlme)
  
  ##
  #Unstructured
  ##
  
  UNfit = gls(observation ~ time, corr = corSymm(form = ~ 1 |id), weights = varIdent(form = ~ 1 | time))

  ##
  #Simple?
  ##
  SIMfit = gls(observation ~ time)
  
  ##
  #Compound Symmetry
  ##
  CSfit = gls(observation ~ time, corr = corCompSymm(form = ~ 1 |id))
  
  ##
  #AR(1)
  ##
  AR1fit = gls(observation ~ time, corr = corAR1(form = ~ 1 | id))
  
  ##
  #CSH
  ##
  CSHfit = gls(observation ~ time, corr = corCompSymm(form = ~ 1 |id), weights = varIdent(form = ~ 1 | time))
  
  ##
  #ARH(1)
  ##
  ARH1fit = gls(observation ~ time,corr = corAR1(form = ~ 1 |id), weights = varIdent(form = ~ 1 | time))
  
  CSfit_AIC = summary(CSfit)$AIC
  AR1fit_AIC = summary(AR1fit)$AIC
  UNfit_AIC = summary(UNfit)$AIC
  CSHfit_AIC = summary(CSHfit)$AIC
  ARH1fit_AIC = summary(ARH1fit)$AIC
  SIMfit_AIC = summary(SIMfit)$AIC
  
  CSfit_BIC = summary(CSfit)$BIC
  AR1fit_BIC = summary(AR1fit)$BIC
  UNfit_BIC = summary(UNfit)$BIC
  CSHfit_BIC = summary(CSHfit)$BIC
  ARH1fit_BIC = summary(ARH1fit)$BIC
  SIMfit_BIC = summary(SIMfit)$BIC
  
  CSfit_AICc = summary(CSfit)$AIC + (2 * (2)) * (2 + 1) / ((n_sub*n_obs) - 2 -
                                                             1)
  AR1fit_AICc = summary(AR1fit)$AIC + (2 * (2)) * (2 + 1) / ((n_sub*n_obs) - 2 -
                                                               1)
  UNfit_AICc = summary(UNfit)$AIC + (2 * (n_obs + choose(n_obs, 2))) * (n_obs + choose(n_obs, 2) + 1) / ((n_sub*n_obs) - (n_obs + choose(n_obs, 2)) -
                                                                                                           1)
  CSHfit_AICc = summary(CSHfit)$AIC + (2 * (n_obs + 1)) * ((n_obs + 1)  + 1) / ((n_sub*n_obs) - (n_obs + 1) -
                                                                                  1)
  ARH1fit_AICc = summary(ARH1fit)$AIC + (2 * (n_obs + 1)) * ((n_obs + 1)  + 1) / ((n_sub*n_obs) - (n_obs + 1) -
                                                                                    1)
  SIMfit_AICc = summary(SIMfit)$AIC + (2 * (1)) * (1 + 1) / ((n_sub*n_obs) - 1 -
                                                               1)
  
  return (
    c(
      UNfit_AIC,
      UNfit_AICc,
      UNfit_BIC,
      SIMfit_AIC,
      SIMfit_AICc,
      SIMfit_BIC,
      CSfit_AIC,
      CSfit_AICc,
      CSfit_BIC,
      AR1fit_AIC,
      AR1fit_AICc,
      AR1fit_BIC,
      CSHfit_AIC,
      CSHfit_AICc,
      CSHfit_BIC,
      ARH1fit_AIC,
      ARH1fit_AICc,
      ARH1fit_BIC
    )
  )
}

####Compilation of Data Gen and Fit, outputs results matrix####
results_matrix <- function(N, n_obs, n_sub, Sigma, means) {
  ##Implement for loop for new input to repeat trials##
  
  ## inputs:
  ## N = # of trials to be run
  ## Sigma = generated sigma matrix/matrices
  
  ICvectorlength = 18
  results = matrix(0, N, ICvectorlength)
  for (i in 1:N) {
    clean_data = generate_data(n_obs, n_sub, Sigma, means)
    ICvector = fit_data(clean_data, n_obs, n_sub)                          #returns vector of AIC and BIC
    results[i,] = ICvector
  }
  
  colnames(results) = c(
    "UNfit_AIC",
    "UNfit_AICc",
    "UNfit_BIC",
    "SIMfit_AIC",
    "SIMfit_AICc",
    "SIMfit_BIC",
    "CSfit_AIC",
    "CSfit_AICc",
    "CSfit_BIC",
    "AR1fit_AIC",
    "AR1fit_AICc",
    "AR1fit_BIC",
    "CSHfit_AIC",
    "CSHfit_AICc",
    "CSHfit_BIC",
    "ARH1fit_AIC",
    "ARH1fit_AICc",
    "ARH1fit_BIC"
  ) 
  rownames(results) = c(1:N)
  return(results)
}

####~~~~~~~~~~~~~~#Data Gen and Collection#~~~~~~~~~~~~~~~~~~####
####Manual (generate your own Sigma) Results Matrix Gen Function to Send to Job####
job_results_gen_manual <- function(N, n_obs, n_sub, Sigma, p, sigma, means, exp_type){
  ##exp type is a string e.g. "CS"
  ##p is a string e.g "(.4, .7, .8)"
  ##sigma is also a string e.g "(1, 1, 1)"
  res = results_matrix(N, n_obs, n_sub, Sigma, means)
  
  means_string = "means"
  for (mean in means){
    means_string = paste(means_string, mean, sep = "_")
  }
  
  file_name = paste("N", N, "obs", n_obs, "sub", n_sub, exp_type, "p", p, "sigma", sigma, means_string, sep = "_")
  file_name = paste(file_name, ".csv", sep = "")
  write.csv(res, file_name)
  return(file_name)
}

####Job Results Gen w/ Sigma Gen Built In####
job_results_gen <- function(N, n_obs, n_sub, exp_type, p, sigma_vect, means){
  ##exp type is a string e.g. "CS"
  ##sigma, and means should be vector type e.g. c(1, 2, 3) or c(1)
  ##for SIM, p = 0
  
  if (exp_type == "UN"){
    print("Do manually with job_results_gen_manual")
  }
  else if (exp_type == "SIM"){
    Sigma = matrix(0, n_obs, n_obs) + diag(sigma_vect^2, n_obs)
  }
  else if (exp_type == "CS"){
    Sigma = makeCS(n_obs, p, sigma_vect)
  }
  else if (exp_type == "AR1"){
    Sigma = makeAR1(n_obs, p, sigma_vect)
  }
  else if (exp_type == "CSH"){
    Sigma = makeCSH(n_obs, p, sigma_vect)
  }
  else if (exp_type == "ARH1"){
    Sigma = makeARH1(n_obs, p, sigma_vect)
  }
  
  res = results_matrix(N, n_obs, n_sub, Sigma, means)
  
  means_string = "means"
  for (mean in means){
    means_string = paste(means_string, mean, sep = "_")
  }
  
  sigma_string = "sigma"
  for (sigma in sigma_vect){
    sigma_string = paste(sigma_string, sigma, sep = "_")
  }
  
  file_name = paste("N", N, "obs", n_obs, "sub", n_sub, exp_type, sigma_string, "p", p, means_string, sep = "_")
  file_name = paste(file_name, ".csv", sep = "")
  write.csv(res, file_name)
  return(file_name)
}

####Line by Line Job Results Gen w/ Sigma Gen####
line_job_results_gen <- function(N, n_obs, n_sub, exp_type, p, sigma_vect, means){
  ##new update!! won't fail anymore :)
  
  means_string = "means"
  for (mean in means){
    means_string = paste(means_string, mean, sep = "_")
  }
  
  sigma_string = "sigma"
  for (sigma in sigma_vect){
    sigma_string = paste(sigma_string, sigma, sep = "_")
  }
  
  results = matrix(0, 0, 18)
  colnames(results) = c(
    "UNfit_AIC",
    "UNfit_AICc",
    "UNfit_BIC",
    "SIMfit_AIC",
    "SIMfit_AICc",
    "SIMfit_BIC",
    "CSfit_AIC",
    "CSfit_AICc",
    "CSfit_BIC",
    "AR1fit_AIC",
    "AR1fit_AICc",
    "AR1fit_BIC",
    "CSHfit_AIC",
    "CSHfit_AICc",
    "CSHfit_BIC",
    "ARH1fit_AIC",
    "ARH1fit_AICc",
    "ARH1fit_BIC"
  )
  
  file_name = paste("N", N, "obs", n_obs, "sub", n_sub, exp_type, sigma_string, "p", p, means_string, sep = "_")
  file_name = paste(file_name, ".csv", sep = "")
  
  write.csv(results, file_name)
  
  row_counter = 0
  
  while(row_counter < N){
    tryCatch(
      expr = {
        if (exp_type == "UN"){
          print("Do manually with job_results_gen_manual")
        }
        else if (exp_type == "SIM"){
          Sigma = matrix(0, n_obs, n_obs) + diag(sigma_vect^2, n_obs)
        }
        else if (exp_type == "CS"){
          Sigma = makeCS(n_obs, p, sigma_vect)
        }
        else if (exp_type == "AR1"){
          Sigma = makeAR1(n_obs, p, sigma_vect)
        }
        else if (exp_type == "CSH"){
          Sigma = makeCSH(n_obs, p, sigma_vect)
        }
        else if (exp_type == "ARH1"){
          Sigma = makeARH1(n_obs, p, sigma_vect)
        }
        
        
        clean_data = generate_data(n_obs, n_sub, Sigma, means)
        res = fit_data(clean_data, n_obs, n_sub)
        
        results = rbind(results, res)
      },
      error = function(e){
        print(e)
        row_counter = row_counter - 1
      },
      finally = {
        row_counter = row_counter + 1
        }
    )
  }
  
  dim_res = dim(results)
  rownames(results) = c(1:dim_res[1])
  print(file_name)
  print(dim_res[1])
  write.csv(results, file_name)
  
}

####UNSTRUCTURED Line by Line Job Results Gen w/ Sigma Gen####
UN_line_job_results_gen <- function(N, n_obs, n_sub, UN_sigma){
  ##new update!! won't fail anymore :)
  
  means = rep(0, n_obs)
  
  UN_type = deparse(substitute(UN_sigma))
  
  results = matrix(0, 0, 18)
  colnames(results) = c(
    "UNfit_AIC",
    "UNfit_AICc",
    "UNfit_BIC",
    "SIMfit_AIC",
    "SIMfit_AICc",
    "SIMfit_BIC",
    "CSfit_AIC",
    "CSfit_AICc",
    "CSfit_BIC",
    "AR1fit_AIC",
    "AR1fit_AICc",
    "AR1fit_BIC",
    "CSHfit_AIC",
    "CSHfit_AICc",
    "CSHfit_BIC",
    "ARH1fit_AIC",
    "ARH1fit_AICc",
    "ARH1fit_BIC"
  )
  
  file_name = paste("N", N, "obs", n_obs, "sub", n_sub, "UN", UN_type, sep = "_")
  file_name = paste(file_name, ".csv", sep = "")
  
  write.csv(results, file_name)
  
  row_counter = 0
  
  while(row_counter < N){
    tryCatch(
      expr = {
        
        clean_data = generate_data(n_obs, n_sub, UN_sigma, means)
        res = fit_data(clean_data, n_obs, n_sub)
        
        results = rbind(results, res)
      },
      error = function(e){
        print(e)
        row_counter = row_counter - 2
      },
      finally = {
        row_counter = row_counter + 1
      }
    )
  }
  
  rownames(results) = c(1:N)
  print(file_name)
  write.csv(results, file_name)
  
}

####Data Retrieval Process Streamlined####
data_retrieve <- function(file_name){
  ##note need to enter file name in quotes e.g. data_retrieve("file_name")
  data = read.csv(file_name, header = TRUE)
  data = as.matrix(data)
  data = data[,2:19]
  return(data)
}

####Mass Data Name Retrieve For CS Data####
homo_mass_data <- function(base_data_name){
  library(stringr)
  # "N_2500_obs_5_sub_50_CS_sigma_1_p_0.1_means_0_0_0_0_0.csv"
  str_pieces = str_split(base_data_name, "sigma_1", n = 3, simplify = TRUE)
  temp = str_split(str_pieces[2], "p_0.1", simplify = TRUE)
  str_pieces[2] = temp[1]
  str_pieces[3] = temp[2]
  
  #data list will be list of 36 filenames to be called
  data_list = c()
  for (sigma in c(1, 3, 5, 10)){
    for (p in c(1:5, 8)){
      filename = paste(str_pieces[1], "sigma_", sigma, str_pieces[2], "p_0.", p, str_pieces[3], sep = "")
      data_list = c(data_list, filename)
    }
  }
  return(data_list)
}
####Mass Data Name Retrieve For Heteroskedastic Data####
mass_data <- function(base_data_name){
  library(stringr)
  str_pieces = str_split(base_data_name, "sigma_1", n = 3, simplify = TRUE)
  temp = str_split(str_pieces[2], "p_0.1", simplify = TRUE)
  str_pieces[2] = temp[1]
  str_pieces[3] = temp[2]
  
  #data list will be list of 36 filenames to be called
  data_list = c()
  for (ratio in c(1:5, 8)){
    for (p in c(1:5, 8)){
      filename = paste(str_pieces[1], "sigma_", ratio, str_pieces[2], "p_0.", p, str_pieces[3], sep = "")
      data_list = c(data_list, filename)
    }
  }
  return(data_list)
}
####~~~~~~~~~~~~~~#Data Analysis Functions#~~~~~~~~~~~~~~~~~~####
####Data Analysis####
#outputs 1-6 ranked least to greatest for each IC
data_analysis <- function(results){
 
  dimension = dim(results)
  N = dimension[1]
  count_IC = dimension[2]
  
  sorted_results = matrix(0, N, 18)
  colnames(sorted_results) = c(paste("AIC", 1:6, sep = "_"),
                               paste("AICc", 1:6, sep = "_"),
                               paste("BIC", 1:6, sep = "_"))
  
  
  for (i in 1:N){
    
    ###AIC###
    
    all_AIC = c(results[i, seq(count_IC, from=1, by=3)])
    
    minAIC = min(all_AIC)
    name = names(all_AIC)
    
    sorted_AIC = sort(all_AIC, decreasing = FALSE)
    sorted_name = names(sorted_AIC)
    final_sorted = c()
    
    for (k in 1:6){
      type = substr(sorted_name[k], 1, nchar(sorted_name[k])-7)
      value = sorted_AIC[k]
      item = paste(type, value, sep = ", ")
      final_sorted = c(final_sorted, item)
    }
    
    sorted_results[i, 1:6] = final_sorted
    
    ###AICc###
    
    all_AICc = c(results[i, seq(count_IC, from=2, by=3)])
    minAICc = min(all_AICc)
    name = names(all_AICc)
    
    sorted_AICc = sort(all_AICc, decreasing = FALSE)
    sorted_name = names(sorted_AICc)
    final_sorted = c()
    
    for (k in 1:6){
      type = substr(sorted_name[k], 1, nchar(sorted_name[k])-8)
      value = sorted_AICc[k]
      item = paste(type, value, sep = ", ")
      final_sorted = c(final_sorted, item)
    }
    
    sorted_results[i, 7:12] = final_sorted
    
    ###BIC###
    
    all_BIC = c(results[i, seq(count_IC, from=3, by=3)])
    minBIC = min(all_BIC)
    name = names(all_BIC)
    
    sorted_BIC = sort(all_BIC, decreasing = FALSE)
    sorted_name = names(sorted_BIC)
    final_sorted = c()
    
    for (k in 1:6){
      type = substr(sorted_name[k], 1, nchar(sorted_name[k])-7)
      value = sorted_BIC[k]
      item = paste(type, value, sep = ", ")
      final_sorted = c(final_sorted, item)
    }
    
    sorted_results[i, 13:18] = final_sorted
    
  }
  
  return(sorted_results)
}
####Difference Matrix####
diff <- function(results, exp_col_num_AIC) {
  
  #@param exp_col_num_AIC: represents column index of AIC_[FIT]
  #competitor-truth
  #so want positive
  #can be 1 (UN), 4 (SIM), 7 (CS), 10 (AR1), 13 (CSH), 16 (ARH1)
  
  N = dim(results)[1]
  
  diff_matrix = matrix(nrow = N, ncol = dim(results)[2])
  
  #set col names
  colnames(diff_matrix) = colnames(results)
  
  #for each row...
  for (i in 1:N) {
    currentRow = results[i,]
    
    ###AIC###
    
    exp_value_AIC = results[i, exp_col_num_AIC]
    
    #for each AIC value in results... count up by 3's
    for (j in seq(18, from=1, by = 3)) {
      #each entry in difference matrix will be a difference between exp
      #and all other covariance structure AIC values
      diff_matrix[i,j] = currentRow[j] - exp_value_AIC
    }
    
    ###AICc###
    
    #now get AICc column index, and value
    exp_col_num_AICc = exp_col_num_AIC + 1
    exp_value_AICc = results[i, exp_col_num_AICc]
    
    #for each AICc value...
    for (r in seq(18, from = 2, by = 3)) {
      diff_matrix[i,r] = currentRow[r] - exp_value_AICc
    }
    
    ###BIC###
    
    #now get BIC column index, and value
    exp_col_num_BIC = exp_col_num_AIC + 2
    exp_value_BIC = results[i, exp_col_num_BIC]
    
    #for each BIC value
    for (k in seq(18, from = 3, by = 3)) {
      diff_matrix[i,k] = currentRow[k] - exp_value_BIC
    }
  }

  
  return (diff_matrix)
}
####Quad Correct Analysis####
quad_correct <- function(diff_matrix, thumb, returnPercents = TRUE){
  N = dim(diff_matrix)[1]
  
  quad_matrix = matrix(0, 4, 18)
  colnames(quad_matrix) = colnames(diff_matrix)
  rownames(quad_matrix) = c(1:4)
  
  for (i in 1:18){
    
    count_1 = 0
    count_2 = 0
    count_3 = 0
    count_4 = 0
    
    for (j in 1:N){
      
      if (diff_matrix[j, i]>0 && diff_matrix[j, i] < thumb){
        ##Right model, not significant
        count_1 = count_1 + 1
      }
      else if (diff_matrix[j, i] < 0 &&
               diff_matrix[j, i] > -(thumb)) {
        ##Wrong model, not significant
        count_2 = count_2 + 1
      }
      else if (diff_matrix[j, i] < 0 &&
               diff_matrix[j, i] < -(thumb)) {
        ##Wrong model, Significant
        count_3 = count_3 + 1
      }
      else if (diff_matrix[j, i] > 0 &&
               diff_matrix[j, i] > thumb) {
        ##Right model, significant
        count_4 = count_4 + 1
      }
    }
    
    if (returnPercents) {
      quad_matrix[1, i] = count_1/N
      quad_matrix[2, i] = count_2/N
      quad_matrix[3, i] = count_3/N
      quad_matrix[4, i] = count_4/N
    }
    else {
      quad_matrix[1, i] = count_1
      quad_matrix[2, i] = count_2
      quad_matrix[3, i] = count_3
      quad_matrix[4, i] = count_4
    }
  }
  
  return(quad_matrix)
}

####New Quad Correct Analysis(Different Method####
new_quad_correct <- function(diff_matrix, thumb, returnPercents = TRUE){
  N = dim(diff_matrix)[1]
  
  count_matrix = matrix(0, 4, 3)
  colnames(count_matrix) = c("AIC", "AICc", "BIC")
  rownames(count_matrix) = c("Type1", "Type2", "Type3", "Type4")
  
  AIC_count_1 = 0
  AIC_count_2 = 0
  AIC_count_3 = 0
  AIC_count_4 = 0
  
  AICc_count_1 = 0
  AICc_count_2 = 0
  AICc_count_3 = 0
  AICc_count_4 = 0
  
  BIC_count_1 = 0
  BIC_count_2 = 0
  BIC_count_3 = 0
  BIC_count_4 = 0
  
  for (i in 1:N){
      sortedAIC = sort(diff_matrix[i, seq(16, from=1, by=3)])
      sortedAICc = sort(diff_matrix[i, seq(17, from=2, by=3)])
      sortedBIC = sort(diff_matrix[i, seq(18, from=3, by=3)])

      ###AIC###
      if (sortedAIC[1] < (-thumb)){
        AIC_count_3 = AIC_count_3 + 1
      }
      if (sortedAIC[1] >= (-thumb) && sortedAIC[1] < 0){
        AIC_count_2 = AIC_count_2 + 1
      }
      if (sortedAIC[1] == 0){
        if (sortedAIC[2] <= thumb){
          AIC_count_1 = AIC_count_1 + 1
        }
        if (sortedAIC[2] > thumb){
          AIC_count_4 = AIC_count_4 + 1
        }
      }
      
      ###AICc###
      if (sortedAICc[1] < (-thumb)){
        AICc_count_3 = AICc_count_3 + 1
      }
      if (sortedAICc[1] >= (-thumb) && sortedAICc[1] < 0){
        AICc_count_2 = AICc_count_2 + 1
      }
      if (sortedAICc[1] == 0){
        if (sortedAICc[2] <= thumb){
          AICc_count_1 = AICc_count_1 + 1
        }
        if (sortedAICc[2] > thumb){
          AICc_count_4 = AICc_count_4 + 1
        }
      }
      
      ###BIC###
      if (sortedBIC[1] < (-thumb)){
        BIC_count_3 = BIC_count_3 + 1
      }
      if (sortedBIC[1] >= (-thumb) && sortedBIC[1] < 0){
        BIC_count_2 = BIC_count_2 + 1
      }
      if (sortedBIC[1] == 0){
        if (sortedBIC[2] <= thumb){
          BIC_count_1 = BIC_count_1 + 1
        }
        if (sortedBIC[2] > thumb){
          BIC_count_4 = BIC_count_4 + 1
        }
      }
  }
  
  count_matrix[1,] = c(AIC_count_1, AICc_count_1, BIC_count_1)
  count_matrix[2,] = c(AIC_count_2, AICc_count_2, BIC_count_2)
  count_matrix[3,] = c(AIC_count_3, AICc_count_3, BIC_count_3)
  count_matrix[4,] = c(AIC_count_4, AICc_count_4, BIC_count_4)
  
  return(count_matrix)
}

####Quad Correct Bar Plot Gen####
quad_correct_graph <- function(quad_data){
  ###outputs a pdf titled "(quad_data_name)_barplots.pdf"
  
  index = c(1,4,7,10,13,16)
  all_names = colnames(quad_data)
  
  dataset = deparse(substitute(quad_data))
  pdf(file=paste(dataset, "_", "barplots", ".pdf", sep = ""))
  
  layout(mat = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2))
  
  for (i in index){
    names = c(all_names[i], all_names[i+1], all_names[i+2])
    barplot(
      width = c(2),
      quad_data[1:4, c(i, (i + 1), (i + 2))],
      names.arg = names,
      col = c("khaki1", "royalblue1", "violetred1", "palegreen1"),
      ylim = c(0, 1.5)
    )
    legend(
      "top",
      legend = c("Type 1", "Type 2", "Type 3", "Type 4"),
      fill = c("khaki1", "royalblue1", "violetred1", "palegreen1"),
      ncol = 4
    )
  }
  
  
  layout(mat = matrix(c(1, 2, 3), nrow = 3, ncol = 1))
  
  for (i in 1:3){
    names = c(all_names[i],
              all_names[i + 3],
              all_names[i + 6],
              all_names[i + 9],
              all_names[i + 12],
              all_names[i + 15])
    barplot(
      width = c(2),
      quad_data[1:4, c(i, (i + 3), (i + 6), (i + 9), (i + 12), (i + 15))],
      names.arg = names,
      col = c("khaki1", "royalblue1", "violetred1", "palegreen1"),
      ylim = c(0, 1.5)
    )
    legend(
      "top",
      legend = c("Type 1", "Type 2", "Type 3", "Type 4"),
      fill = c("khaki1", "royalblue1", "violetred1", "palegreen1"),
      ncol = 4
    )
  }
  
  dev.off()
  
}
####Distribution Histogram Gen Functions####
AICDistribution_manual <- function(N, n_subs, Sigma, means) {
  #matrix with AIC scores for different fits as columns
  AIC_scores = matrix(nrow = N, ncol = 6)
  names = c("Unstructured", "Simple", "CS", "AR1", "CSH", "ARH1")
  
  #for each row in AIC_scores...
  for (i in 1:N) {
    n_obs = dim(Sigma)[1]
    
    #generate data
    data = generate_data(n_obs, n_subs, Sigma, means)
    #generate fits
    fits = fit_data(data, n_obs, n_subs)
    
    UN_AIC = fits[1]
    SIM_AIC = fits[4]
    CS_AIC = fits[7]
    AR1_AIC = fits[10]
    CSH_AIC = fits[13]
    ARH1_AIC = fits[16]
    
    #collect AIC scores for each fit
    AICVector = c(UN_AIC, SIM_AIC, CS_AIC, AR1_AIC, CSH_AIC, ARH1_AIC)
    
    #put the AIC scores into AIC_scores matrix
    for (j in 1:6) {
      AIC_scores[i, j] = AICVector[j]
    }
  }
  
  #histogram each column
  for (i in 1:6) {
    hist(AIC_scores[, i], xlab = paste("AIC values for", names[i], sep = " "), main = paste("Distribution of AIC", names[i], "Values", sep = " "))
  }
}

AICcDistribution_manual <- function(N, n_subs, Sigma, means) {
  #matrix with AICc scores for different fits as columns
  AICc_scores = matrix(nrow = N, ncol = 6)
  names = c("Unstructured", "Simple", "CS", "AR1", "CSH", "ARH1")
  
  #for each row in AICc_scores...
  for (i in 1:N) {
    n_obs = dim(Sigma)[1]
    
    #generate data
    data = generate_data(n_obs, n_subs, Sigma, means)
    #generate fits
    fits = fit_data(data, n_obs, n_subs)
    
    UN_AICc = fits[2]
    SIM_AICc = fits[5]
    CS_AICc = fits[8]
    AR1_AICc = fits[11]
    CSH_AICc = fits[14]
    ARH1_AICc = fits[17]
    
    #collect AICc scores for each fit
    AICcVector = c(UN_AICc, SIM_AICc, CS_AICc, AR1_AICc, CSH_AICc, ARH1_AICc)
    
    #put the AICc scores into AICc_scores matrix
    for (j in 1:6) {
      AICc_scores[i, j] = AICcVector[j]
    }
  }
  
  #histogram each column
  for (i in 1:6) {
    hist(AICc_scores[, i], xlab = paste("AICc values for", names[i], sep = " "), main = paste("Distribution of AICc", names[i], "Values", sep = " "))
  }
}

BICDistribution_manual <- function(N, n_subs, Sigma, means) {
  #matrix with BIC scores for different fits as columns
  BIC_scores = matrix(nrow = N, ncol = 6)
  names = c("Unstructured", "Simple", "CS", "AR1", "CSH", "ARH1")
  
  #for each row in BIC_scores...
  for (i in 1:N) {
    n_obs = dim(Sigma)[1]
    
    #generate data
    data = generate_data(n_obs, n_subs, Sigma, means)
    #generate fits
    fits = fit_data(data, n_obs, n_subs)
    
    UN_BIC = fits[3]
    SIM_BIC = fits[6]
    CS_BIC = fits[9]
    AR1_BIC = fits[12]
    CSH_BIC = fits[15]
    ARH1_BIC = fits[18]
    
    #collect BIC scores for each fit
    BICVector = c(UN_BIC, SIM_BIC, CS_BIC, AR1_BIC, CSH_BIC, ARH1_BIC)
    
    #put the BIC scores into AIC_scores matrix
    for (j in 1:6) {
      BIC_scores[i, j] = BICVector[j]
    }
  }
  
  #histogram each column
  for (i in 1:6) {
    hist(BIC_scores[, i], xlab = paste("BIC values for", names[i], sep = " "), main = paste("Distribution of BIC", names[i], "Values", sep = " "))
  }
}

AICDistribution <- function(results) {
  N = dim(results)[1]
  names = c("Unstructured", "Simple", "CS", "AR1", "CSH", "ARH1")
  
  AICMatrix = matrix(nrow = N, ncol = 6)
  for (i in 1:N) {
    k = 1
    for (j in seq(18, from = 1, by = 3)) {
      AICMatrix[i, k] = results[i, j]
      k = k + 1
    }
  }
  
  #histogram each column
  for (i in 1:6) {
    hist(AICMatrix[, i], xlab = paste("AIC values for", names[i], sep = " "), main = paste("Distribution of AIC", names[i], "Values", sep = " "))
  }
}

AICcDistribution <- function(results) {
  N = dim(results)[1]
  names = c("Unstructured", "Simple", "CS", "AR1", "CSH", "ARH1")
  
  AICcMatrix = matrix(nrow = N, ncol = 6)
  for (i in 1:N) {
    k = 1
    for (j in seq(18, from = 3, by = 3)) {
      AICcMatrix[i, k] = results[i, j]
      k = k + 1
    }
  }
  
  #histogram each column
  for (i in 1:6) {
    hist(AICcMatrix[, i], xlab = paste("AICc values for", names[i], sep = " "), main = paste("Distribution of AICc", names[i], "Values", sep = " "))
  }
}

BICDistribution <- function(results) {
  N = dim(results)[1]
  names = c("Unstructured", "Simple", "CS", "AR1", "CSH", "ARH1")
  
  BICMatrix = matrix(nrow = N, ncol = 6)
  for (i in 1:N) {
    k = 1
    for (j in seq(18, from = 2, by = 3)) {
      BICMatrix[i, k] = results[i, j]
      k = k + 1
    }
  }
  
  #histogram each column
  for (i in 1:6) {
    hist(BICMatrix[, i], xlab = paste("BIC values for", names[i], sep = " "), main = paste("Distribution of BIC", names[i], "Values", sep = " "))
  }
}


####Makes Box Plots for Experiment Shown in Thesis (RHO)####
#section 6
AICRhoBoxplots <- function() {
  AIC_UN_Data = matrix(nrow = 100, ncol = 10)
  AIC_SIM_Data = matrix(nrow = 100, ncol = 10)
  AIC_CS_Data = matrix(nrow = 100, ncol = 10)
  AIC_AR1_Data = matrix(nrow = 100, ncol = 10)
  AIC_CSH_Data = matrix(nrow = 100, ncol = 10)
  AIC_ARH1_Data = matrix(nrow = 100, ncol = 10)
  
  #we vary rho from 0 to .9
  for (i in 0:9) {
    #these are the parameters detailed in section 6
    Sigma = makeCS(3, i / 10, 1.2)
    res = results_matrix(100, 3, 40, Sigma, c(0, 0, 0))
    differences = diff(res, 7)
    
    AIC_UN_Data[, i+1] = differences[, 1]
    AIC_SIM_Data[, i+1] = differences[, 4]
    AIC_CS_Data[, i+1] = differences[, 7]
    AIC_AR1_Data[, i+1] = differences[, 10]
    AIC_CSH_Data[, i+1] = differences[, 13]
    AIC_ARH1_Data[, i+1] = differences[, 16]
  }
  
  #create boxplots
  boxplot(AIC_UN_Data[,1:10], main = "AIC_UN - AIC_CS")
  boxplot(AIC_SIM_Data[,1:10], main = "AIC_SIM - AIC_CS")
  boxplot(AIC_CS_Data[,1:10], main = "AIC_CS - AIC_CS")
  boxplot(AIC_AR1_Data[,1:10], main = "AIC_AR1 - AIC_CS")
  boxplot(AIC_CSH_Data[,1:10], main = "AIC_CSH - AIC_CS")
  boxplot(AIC_ARH1_Data[,1:10], main = "AIC_ARH1 - AIC_CS")
}

####Makes Box Plots for Experiment Shown in Thesis (SIGMA)####
#section 6
AICSigmaBoxplots <- function() {
  AIC_UN_Data = matrix(nrow = 100, ncol = 10)
  AIC_SIM_Data = matrix(nrow = 100, ncol = 10)
  AIC_CS_Data = matrix(nrow = 100, ncol = 10)
  AIC_AR1_Data = matrix(nrow = 100, ncol = 10)
  AIC_CSH_Data = matrix(nrow = 100, ncol = 10)
  AIC_ARH1_Data = matrix(nrow = 100, ncol = 10)
  
  #we vary sigma from 1 to 1.9. Rho = .5
  for (i in 0:9) {
    #these are the parameters detailed in section 6
    Sigma = makeCS(3, .5, 1+i/10)
    res = results_matrix(100, 3, 40, Sigma, c(0, 0, 0))
    differences = diff(res, 7)
    
    AIC_UN_Data[, i+1] = differences[, 1]
    AIC_SIM_Data[, i+1] = differences[, 4]
    AIC_CS_Data[, i+1] = differences[, 7]
    AIC_AR1_Data[, i+1] = differences[, 10]
    AIC_CSH_Data[, i+1] = differences[, 13]
    AIC_ARH1_Data[, i+1] = differences[, 16]
  }
  
  #create boxplots
  boxplot(AIC_UN_Data[,1:10], main = "AIC_UN - AIC_CS")
  boxplot(AIC_SIM_Data[,1:10], main = "AIC_SIM - AIC_CS")
  boxplot(AIC_CS_Data[,1:10], main = "AIC_CS - AIC_CS")
  boxplot(AIC_AR1_Data[,1:10], main = "AIC_AR1 - AIC_CS")
  boxplot(AIC_CSH_Data[,1:10], main = "AIC_CSH - AIC_CS")
  boxplot(AIC_ARH1_Data[,1:10], main = "AIC_ARH1 - AIC_CS")
}

####Overlapping Histograms####
overlap_histograms <- function(data, exp_col_num_AIC){
  ##generates 5 double histograms for each IC which compare expected type
  ##distribution to each other distribution
  ##outputs as a pdf called "dataname_histograms.pdf"
  
  dataset = deparse(substitute(data))
  
  c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
  c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
  
  pdf(file=paste(dataset, "_", "histograms", ".pdf", sep = ""))
  
  layout(mat = matrix(c(1, 2, 3, 4, 5, 0), nrow = 3, ncol = 2))
  
  ###AIC###
  AIC_cols = c(1, 4, 7, 10, 13, 16)
  
  for (col in AIC_cols){
    if (col != exp_col_num_AIC){
      min = min(c(data[,exp_col_num_AIC], data[,col])) - 10
      max = max(c(data[,exp_col_num_AIC], data[,col])) + 10
      
      breakpoints = pretty(min:max, n = 25)
      names = colnames(data)
      
      hist1 = hist(data[, exp_col_num_AIC], breaks = breakpoints, plot = FALSE)
      hist2 = hist(data[, col], breaks = breakpoints, plot = FALSE)
      plot(
        hist1,
        xlab = paste("Blue:", names[exp_col_num_AIC], "Pink:", names[col], sep = " "),
        main = paste(dataset, "Visualization of AIC Comparison", sep = ": "),
        col = c1
      )
      plot(hist2, col = c2, add = TRUE)
    }
    else if (col == exp_col_num_AIC){
      
    }
  }
  
  ###AICc###
  AICc_cols = c(2, 5, 8, 11, 14, 17)
  exp_col_num_AICc = exp_col_num_AIC + 1
  
  for (col in AICc_cols){
    if (col != exp_col_num_AICc){
      min = min(c(data[,exp_col_num_AICc], data[,col])) - 10
      max = max(c(data[,exp_col_num_AICc], data[,col])) + 10
      
      breakpoints = pretty(min:max, n = 25)
      names = colnames(data)
      
      hist1 = hist(data[, exp_col_num_AICc], breaks = breakpoints, plot = FALSE)
      hist2 = hist(data[, col], breaks = breakpoints, plot = FALSE)
      plot(
        hist1,
        xlab = paste("Blue:", names[exp_col_num_AICc], "Pink:", names[col], sep = " "),
        main = paste(dataset, "Visualization of AICc Comparison", sep = ": "),
        col = c1
      )
      plot(hist2, col = c2, add = TRUE)
    }
    else if (col == exp_col_num_AICc){
      
    }
  }
  
  ###BIC###
  BIC_cols = c(3, 6, 9, 12, 15, 18)
  exp_col_num_BIC = exp_col_num_AIC + 2
  
  for (col in BIC_cols){
    if (col != exp_col_num_BIC){
      min = min(c(data[,exp_col_num_BIC], data[,col])) - 10
      max = max(c(data[,exp_col_num_BIC], data[,col])) + 10
      
      breakpoints = pretty(min:max, n = 25)
      names = colnames(data)
      
      hist1 = hist(data[, exp_col_num_BIC], breaks = breakpoints, plot = FALSE)
      hist2 = hist(data[, col], breaks = breakpoints, plot = FALSE)
      plot(
        hist1,
        xlab = paste("Blue:", names[exp_col_num_BIC], "Pink:", names[col], sep = " "),
        main = paste(dataset, "Visualization of BIC Comparison", sep = ": "),
        col = c1
      )
      plot(hist2, col = c2, add = TRUE)
    }
    else if (col == exp_col_num_BIC){
      
    }
  }
  dev.off()
}
####Plot 3 and 4 as a Function of a Desired Lever####
plot34 <- function(data_list, x_vect, x_vect_var_name, exp_col_num_AIC, thumb){
  ###data_list is a list(data1, data2, ... dataN) list of all sets of data of interest
  ###IMPORTANT ~ save ur data_list and name it before use in this function
  ###b/c this name is used to create the pdf file name
  ###x_vect_var_name examples ("sigma", "rho", etc.) this becomes the x lab of each plot
  ###x_vect is a vect of the variable you want to plot the types 3 and 4 by
  ###for example for varying p would input  x_vect = c(.1, .2, .3, .4, .5, .8)
  ###which plots lines for 3 and 4 as a function of p
  ###length of data_list and x_vect should be equal
  
  data_name = deparse(substitute(data_list))
  
  pdf(file=paste(data_name, "_", "type34", ".pdf", sep = ""))
  graphics::layout(mat = matrix(c(1, 3, 5, 2, 4, 6), nrow = 3, ncol = 2))
  
  len = length(data_list)
  
  #AIC_col = c(1, 4, 7, 10, 13, 16)
  
  ###rows of these matrices will go top to bottom p = .1, .2, ... , .5, .8
  type3_data = matrix(0, len, 18)
  colnames(type3_data) = c(
    "AIC_UN",
    "AIC_SIM",
    "AIC_CS",
    "AIC_AR1",
    "AIC_CSH",
    "AIC_ARH1",
    "AICc_UN",
    "AICc_SIM",
    "AICc_CS",
    "AICc_AR1",
    "AICc_CSH",
    "AICc_ARH1",
    "BIC_UN",
    "BIC_SIM",
    "BIC_CS",
    "BIC_AR1",
    "BIC_CSH",
    "BIC_ARH1"
  )
  
  type4_data = matrix(0, len, 18)
  colnames(type4_data) = colnames(type3_data)
  
  row_count = 1
  for (data in data_list){
    diff = diff(data, exp_col_num_AIC)
    quad = quad_correct(diff, thumb)

    ##
    #type 3
    ##
    
    #AIC
    type3_data[row_count, 1] = quad[3, 1]
    type3_data[row_count, 2] = quad[3, 4]
    type3_data[row_count, 3] = quad[3, 7]
    type3_data[row_count, 4] = quad[3, 10]
    type3_data[row_count, 5] = quad[3, 13]
    type3_data[row_count, 6] = quad[3, 16]
    
    #AICc
    type3_data[row_count, 7] = quad[3, 2]
    type3_data[row_count, 8] = quad[3, 5]
    type3_data[row_count, 9] = quad[3, 8]
    type3_data[row_count, 10] = quad[3, 11]
    type3_data[row_count, 11] = quad[3, 14]
    type3_data[row_count, 12] = quad[3, 17]
    
    #BIC
    type3_data[row_count, 13] = quad[3, 3]
    type3_data[row_count, 14] = quad[3, 6]
    type3_data[row_count, 15] = quad[3, 9]
    type3_data[row_count, 16] = quad[3, 12]
    type3_data[row_count, 17] = quad[3, 15]
    type3_data[row_count, 18] = quad[3, 18]
    
    ##
    #type 4
    ##
    
    #AIC
    type4_data[row_count, 1] = quad[4, 1]
    type4_data[row_count, 2] = quad[4, 4]
    type4_data[row_count, 3] = quad[4, 7]
    type4_data[row_count, 4] = quad[4, 10]
    type4_data[row_count, 5] = quad[4, 13]
    type4_data[row_count, 6] = quad[4, 16]
    
    #AICc
    type4_data[row_count, 7] = quad[4, 2]
    type4_data[row_count, 8] = quad[4, 5]
    type4_data[row_count, 9] = quad[4, 8]
    type4_data[row_count, 10] = quad[4, 11]
    type4_data[row_count, 11] = quad[4, 14]
    type4_data[row_count, 12] = quad[4, 17]
    
    #BIC
    type4_data[row_count, 13] = quad[4, 3]
    type4_data[row_count, 14] = quad[4, 6]
    type4_data[row_count, 15] = quad[4, 9]
    type4_data[row_count, 16] = quad[4, 12]
    type4_data[row_count, 17] = quad[4, 15]
    type4_data[row_count, 18] = quad[4, 18]
    
    row_count = row_count + 1
  }
  
  
  
  ##
  #AIC PLOTS
  ##
  
  #Type3
  matplot(
    x_vect,
    type3_data[, 1:6],
    main = "AIC Type 3",
    xlab = x_vect_var_name,
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  #Type4
  matplot(
    x_vect,
    type4_data[, 1:6],
    main = "AIC Type 4",
    xlab = x_vect_var_name,
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"), 
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  ##
  #AICc PLOTS
  ##
  
  #Type3
  matplot(
    x_vect,
    type3_data[, 7:12],
    main = "AICc Type 3",
    xlab = x_vect_var_name,
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    box.lty = 0,
    ncol = 3,
    cex = .6
  )
  
  #Type4
  matplot(
    x_vect,
    type4_data[, 7:12],
    main = "AICc Type 4",
    xlab = x_vect_var_name,
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"), 
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  ##
  #BIC PLOTS
  ##
  
  #Type3
  matplot(
    x_vect,
    type3_data[, 13:18],
    main = "BIC Type 3",
    xlab = x_vect_var_name,
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  #Type4
  matplot(
    x_vect,
    type4_data[, 13:18],
    main = "BIC Type 4",
    xlab = x_vect_var_name,
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"), 
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  dev.off()
}

####Plot 3 and 4 as a Function of RULE OF THUMB (also attatches winner plots)####
thumb_plot34 <- function(data, data_name, exp_col_num_AIC){
  ###data_list is a list(data1, data2, ... dataN) list of all sets of data of interest
  ###IMPORTANT ~ save ur data_list and name it before use in this function
  ###b/c this name is used to create the pdf file name
  ###thumb_vect is a vect of the thumb vals want to investigate
  
  #data_name = deparse(substitute(data))
  
  
  pdf(file=paste(data_name, "_", "thumb", ".pdf", sep = ""))
  graphics::layout(mat = matrix(c(1, 3, 5, 2, 4, 6), nrow = 3, ncol = 2))

  
  thumb_vect = seq(0, 7, .5)
  
  ###rows of these matrices will go top to bottom p = .1, .2, ... , .5, .8
  type3_data = matrix(0, length(thumb_vect), 18)
  colnames(type3_data) = c(
    "AIC_UN",
    "AIC_SIM",
    "AIC_CS",
    "AIC_AR1",
    "AIC_CSH",
    "AIC_ARH1",
    "AICc_UN",
    "AICc_SIM",
    "AICc_CS",
    "AICc_AR1",
    "AICc_CSH",
    "AICc_ARH1",
    "BIC_UN",
    "BIC_SIM",
    "BIC_CS",
    "BIC_AR1",
    "BIC_CSH",
    "BIC_ARH1"
  )
  
  type4_data = matrix(0, length(thumb_vect), 18)
  colnames(type4_data) = colnames(type3_data)
  
  row_count = 1
  diff = diff(data, exp_col_num_AIC)
  
  for (thumb in thumb_vect){
    quad = quad_correct(diff, thumb)
    
    ##
    #type 3
    ##
    
    #AIC
    type3_data[row_count, 1] = quad[3, 1]
    type3_data[row_count, 2] = quad[3, 4]
    type3_data[row_count, 3] = quad[3, 7]
    type3_data[row_count, 4] = quad[3, 10]
    type3_data[row_count, 5] = quad[3, 13]
    type3_data[row_count, 6] = quad[3, 16]
    
    #AICc
    type3_data[row_count, 7] = quad[3, 2]
    type3_data[row_count, 8] = quad[3, 5]
    type3_data[row_count, 9] = quad[3, 8]
    type3_data[row_count, 10] = quad[3, 11]
    type3_data[row_count, 11] = quad[3, 14]
    type3_data[row_count, 12] = quad[3, 17]
    
    #BIC
    type3_data[row_count, 13] = quad[3, 3]
    type3_data[row_count, 14] = quad[3, 6]
    type3_data[row_count, 15] = quad[3, 9]
    type3_data[row_count, 16] = quad[3, 12]
    type3_data[row_count, 17] = quad[3, 15]
    type3_data[row_count, 18] = quad[3, 18]
    
    ##
    #type 4
    ##
    
    #AIC
    type4_data[row_count, 1] = quad[4, 1]
    type4_data[row_count, 2] = quad[4, 4]
    type4_data[row_count, 3] = quad[4, 7]
    type4_data[row_count, 4] = quad[4, 10]
    type4_data[row_count, 5] = quad[4, 13]
    type4_data[row_count, 6] = quad[4, 16]
    
    #AICc
    type4_data[row_count, 7] = quad[4, 2]
    type4_data[row_count, 8] = quad[4, 5]
    type4_data[row_count, 9] = quad[4, 8]
    type4_data[row_count, 10] = quad[4, 11]
    type4_data[row_count, 11] = quad[4, 14]
    type4_data[row_count, 12] = quad[4, 17]
    
    #BIC
    type4_data[row_count, 13] = quad[4, 3]
    type4_data[row_count, 14] = quad[4, 6]
    type4_data[row_count, 15] = quad[4, 9]
    type4_data[row_count, 16] = quad[4, 12]
    type4_data[row_count, 17] = quad[4, 15]
    type4_data[row_count, 18] = quad[4, 18]
    
    
    row_count = row_count + 1
  }
  
  
  
  ##
  #AIC PLOTS
  ##
  
  #Type3
  matplot(
    thumb_vect,
    type3_data[, 1:6],
    main = "AIC Type 3",
    xlab = "Rule of Thumb",
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  # abline(v=2, col="red")
  # abline(v=3, col="blue")
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  #Type4
  matplot(
    thumb_vect,
    type4_data[, 1:6],
    main = "AIC Type 4",
    xlab = "Rule of Thumb",
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  # abline(v=2, col="red")
  # abline(v=3, col="blue")
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"), 
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  ##
  #AICc PLOTS
  ##
  
  #Type3
  matplot(
    thumb_vect,
    type3_data[, 7:12],
    main = "AICc Type 3",
    xlab = "Rule of Thumb",
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  # abline(v=2, col="red")
  # abline(v=3, col="blue")
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    box.lty = 0,
    ncol = 3,
    cex = .6
  )
  
  #Type4
  matplot(
    thumb_vect,
    type4_data[, 7:12],
    main = "AICc Type 4",
    xlab = "Rule of Thumb",
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  # abline(v=2, col="red")
  # abline(v=3, col="blue")
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"), 
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  ##
  #BIC PLOTS
  ##
  
  #Type3
  matplot(
    thumb_vect,
    type3_data[, 13:18],
    main = "BIC Type 3",
    xlab = "Rule of Thumb",
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  abline(v=1, col="red")
  abline(v=3, col="green")
  abline(v=5, col="blue")
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  #Type4
  matplot(
    thumb_vect,
    type4_data[, 13:18],
    main = "BIC Type 4",
    xlab = "Rule of Thumb",
    ylab = "Proportion",
    pch = 19,
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    type = "b"
  )
  abline(v=1, col="red")
  abline(v=3, col="green")
  abline(v=5, col="blue")
  legend(
    "topleft",
    legend = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1"),
    fill = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"), 
    box.lty=0,
    ncol = 3,
    cex = .6
  )
  
  ####################3
  graphics::layout(mat = matrix(c(1, 2, 3), nrow = 3, ncol = 1))
  
  sorted_results = data_analysis(data)
  N = dim(data)[1]
  
  winner_matrix = matrix(0, 3, 6)
  colnames(winner_matrix) = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1")
  rownames(winner_matrix) = c("AIC", "AICc", "BIC")
  
  for (row in 1:N){
    for (col in c(1:3)){
      col_num = c(1, 7, 13)
      string = sorted_results[row, col_num[col]]
      string = str_split(string, ", ", simplify = TRUE)
      if (string[1] == "UN"){
        winner_matrix[col, 1] = winner_matrix[col, 1] + 1
      }
      if (string[1] == "SIM"){
        winner_matrix[col, 2] = winner_matrix[col, 2] + 1
      }
      if (string[1] == "CS"){
        winner_matrix[col, 3] = winner_matrix[col, 3] + 1
      }
      if (string[1] == "AR1"){
        winner_matrix[col, 4] = winner_matrix[col, 4] + 1
      }
      if (string[1] == "CSH"){
        winner_matrix[col, 5] = winner_matrix[col, 5] + 1
      }
      if (string[1] == "ARH1"){
        winner_matrix[col, 6] = winner_matrix[col, 6] + 1
      }
    }
  }
  print(winner_matrix)
  plot1 = barplot(
    winner_matrix[1, ],
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    main = "AIC Winner Dist.", 
    ylim = c(0,2500)
  )
  text(plot1, winner_matrix[1,] + 170, winner_matrix[1,], font=2, col= "black")
  plot2 = barplot(
    winner_matrix[2, ],
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    main = "AICc Winner Dist.", 
    ylim = c(0,2500)
  )
  text(plot2, winner_matrix[2,] + 170, winner_matrix[2,], font=2, col= "black")
  plot3 = barplot(
    winner_matrix[3, ],
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    main = "BIC Winner Dist.", 
    ylim = c(0,2500)
  )
  text(plot3, winner_matrix[3,] + 170, winner_matrix[3,], font=2, col= "black")
  
  
  
  dev.off()
}
####Winner Proportion Barplots####
winner_barplots <- function(data, data_name){
  pdf(file=paste("winners", "_", data_name, ".pdf", sep = ""))
  graphics::layout(mat = matrix(c(1, 2, 3), nrow = 3, ncol = 1))
  
  sorted_results = data_analysis(data)
  N = dim(data)[1]
  
  winner_matrix = matrix(0, 3, 6)
  colnames(winner_matrix) = c("UN", "SIM", "CS", "AR1", "CSH", "ARH1")
  rownames(winner_matrix) = c("AIC", "AICc", "BIC")
  
  for (row in 1:N){
    for (col in c(1:3)){
      col_num = c(1, 7, 13)
      string = sorted_results[row, col_num[col]]
      string = str_split(string, ", ", simplify = TRUE)
      if (string[1] == "UN"){
        winner_matrix[col, 1] = winner_matrix[col, 1] + 1
      }
      if (string[1] == "SIM"){
        winner_matrix[col, 2] = winner_matrix[col, 2] + 1
      }
      if (string[1] == "CS"){
        winner_matrix[col, 3] = winner_matrix[col, 3] + 1
      }
      if (string[1] == "AR1"){
        winner_matrix[col, 4] = winner_matrix[col, 4] + 1
      }
      if (string[1] == "CSH"){
        winner_matrix[col, 5] = winner_matrix[col, 5] + 1
      }
      if (string[1] == "ARH1"){
        winner_matrix[col, 6] = winner_matrix[col, 6] + 1
      }
    }
  }
  print(winner_matrix)
  plot1 = barplot(
    winner_matrix[1, ],
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    main = "AIC Winner Dist.", 
    ylim = c(0,2500)
  )
  text(plot1, winner_matrix[1,] + 170, winner_matrix[1,], font=2, col= "black")
  plot2 = barplot(
    winner_matrix[2, ],
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    main = "AICc Winner Dist.", 
    ylim = c(0,2500)
  )
  text(plot2, winner_matrix[2,] + 170, winner_matrix[2,], font=2, col= "black")
  plot3 = barplot(
    winner_matrix[3, ],
    col = c("violetred1", "orange", "green", "royalblue1", "purple", "pink"),
    main = "BIC Winner Dist.", 
    ylim = c(0,2500)
  )
  text(plot3, winner_matrix[3,] + 170, winner_matrix[3,], font=2, col= "black")
  
  dev.off()
}
####3D Plots for heterosked. vs. rho vs. type3/4####
graph_3D <- function(base_data_name, AIC_exp_num, thumb){
  ###!!!!!!!important if this is ur first time using
  ###make sure to un-comment the install packages and library lines
  
  ##base_data_name looks like this
  ##"N_2500_obs_5_sub_50_ARH1_sigma_1_10_10_10_10_p_0.1_means_0_0_0_0_0.csv"
  ##as long as you have all 36 necessary data files on ur desktop this func
  ##will retrieve all neccesary data just from this base name
  ##rn this only works for hetero models if obs is same and are just 
  ##varying rho and heteroskedasticity ratio
  
  #install.packages("plotly")
  library(plotly)
  #install.packages("stringr")
  library("stringr")
  
  type3_AIC = matrix(0, 6, 6)
  colnames(type3_AIC) = c("p .1" , "p .2", "p .3", "p .4", "p .5", "p .8")
  rownames(type3_AIC) = c("ratio .1" , "ratio .2", "ratio .3", "ratio .4", "ratio .5", "ratio .8")
  
  type3_AIC_UN = type3_AIC
  type3_AIC_SIM = type3_AIC
  type3_AIC_CS = type3_AIC
  type3_AIC_AR1 = type3_AIC
  type3_AIC_CSH = type3_AIC
  type3_AIC_ARH1 = type3_AIC
  
  type4_AIC_UN = type3_AIC
  type4_AIC_SIM = type3_AIC
  type4_AIC_CS = type3_AIC
  type4_AIC_AR1 = type3_AIC
  type4_AIC_CSH = type3_AIC
  type4_AIC_ARH1 = type3_AIC
  
  type3_AICc_UN = type3_AIC
  type3_AICc_SIM = type3_AIC
  type3_AICc_CS = type3_AIC
  type3_AICc_AR1 = type3_AIC
  type3_AICc_CSH = type3_AIC
  type3_AICc_ARH1 = type3_AIC
  
  type4_AICc_UN = type3_AIC
  type4_AICc_SIM = type3_AIC
  type4_AICc_CS = type3_AIC
  type4_AICc_AR1 = type3_AIC
  type4_AICc_CSH = type3_AIC
  type4_AICc_ARH1 = type3_AIC
  
  type3_BIC_UN = type3_AIC
  type3_BIC_SIM = type3_AIC
  type3_BIC_CS = type3_AIC
  type3_BIC_AR1 = type3_AIC
  type3_BIC_CSH = type3_AIC
  type3_BIC_ARH1 = type3_AIC
  
  type4_BIC_UN = type3_AIC
  type4_BIC_SIM = type3_AIC
  type4_BIC_CS = type3_AIC
  type4_BIC_AR1 = type3_AIC
  type4_BIC_CSH = type3_AIC
  type4_BIC_ARH1 = type3_AIC
  
  ##data list fills in reading style left to right top to bottom
  ## eg p1_het1 ... p8_het1, p1_het2, ... 
  
  str_pieces = str_split(base_data_name, "sigma_1", n = 3, simplify = TRUE)
  temp = str_split(str_pieces[2], "p_0.1", simplify = TRUE)
  str_pieces[2] = temp[1]
  str_pieces[3] = temp[2]
  
  #data list will be list of 36 filenames to be called
  data_list = c()
  for (ratio in c(1:5, 8)){
    for (p in c(1:5, 8)){
      filename = paste(str_pieces[1], "sigma_", ratio, str_pieces[2], "p_0.", p, str_pieces[3], sep = "")
      data_list = c(data_list, filename)
    }
  }
  
  ##this is ugly sorry y'all
  col_index = rep(1:6, 6)
  row_index = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 6), rep(6, 6))
  for (i in 1:36){
    file = data_list[i]
    data = data_retrieve(file)
    
    diff = diff(data, AIC_exp_num)
    quad = quad_correct(diff, thumb)
    # print(quad)
    
    type3_AIC_UN[col_index[i], row_index[i]] = quad[3, 1]
    type3_AIC_SIM[col_index[i], row_index[i]] = quad[3, 4]
    type3_AIC_CS[col_index[i], row_index[i]] = quad[3, 7]
    type3_AIC_AR1[col_index[i], row_index[i]] = quad[3, 10]
    type3_AIC_CSH[col_index[i], row_index[i]] = quad[3, 13]
    type3_AIC_ARH1[col_index[i], row_index[i]] = quad[3, 16]
    
    type4_AIC_UN[col_index[i], row_index[i]] = quad[4, 1]
    type4_AIC_SIM[col_index[i], row_index[i]] = quad[4, 4]
    type4_AIC_CS[col_index[i], row_index[i]] = quad[4, 7]
    type4_AIC_AR1[col_index[i], row_index[i]] = quad[4, 10]
    type4_AIC_CSH[col_index[i], row_index[i]] = quad[4, 13]
    type4_AIC_ARH1[col_index[i], row_index[i]] = quad[4, 16]
    
    type3_AICc_UN[col_index[i], row_index[i]] = quad[3, 2]
    type3_AICc_SIM[col_index[i], row_index[i]] = quad[3, 5]
    type3_AICc_CS[col_index[i], row_index[i]] = quad[3, 8]
    type3_AICc_AR1[col_index[i], row_index[i]] = quad[3, 11]
    type3_AICc_CSH[col_index[i], row_index[i]] = quad[3, 14]
    type3_AICc_ARH1[col_index[i], row_index[i]] = quad[3, 17]
    
    type4_AICc_UN[col_index[i], row_index[i]] = quad[4, 2]
    type4_AICc_SIM[col_index[i], row_index[i]] = quad[4, 5]
    type4_AICc_CS[col_index[i], row_index[i]] = quad[4, 8]
    type4_AICc_AR1[col_index[i], row_index[i]] = quad[4, 11]
    type4_AICc_CSH[col_index[i], row_index[i]] = quad[4, 14]
    type4_AICc_ARH1[col_index[i], row_index[i]] = quad[4, 17]
    
    type3_BIC_UN[col_index[i], row_index[i]] = quad[3, 3]
    type3_BIC_SIM[col_index[i], row_index[i]] = quad[3, 6]
    type3_BIC_CS[col_index[i], row_index[i]] = quad[3, 9]
    type3_BIC_AR1[col_index[i], row_index[i]] = quad[3, 12]
    type3_BIC_CSH[col_index[i], row_index[i]] = quad[3, 15]
    type3_BIC_ARH1[col_index[i], row_index[i]] = quad[3, 18]
    
    type4_BIC_UN[col_index[i], row_index[i]] = quad[4, 3]
    type4_BIC_SIM[col_index[i], row_index[i]] = quad[4, 6]
    type4_BIC_CS[col_index[i], row_index[i]] = quad[4, 9]
    type4_BIC_AR1[col_index[i], row_index[i]] = quad[4, 12]
    type4_BIC_CSH[col_index[i], row_index[i]] = quad[4, 15]
    type4_BIC_ARH1[col_index[i], row_index[i]] = quad[4, 18]
  }
  
  axx <- list(
    title = "SD Ratio Min/Max (Heterosked.)"
  )
  
  axy <- list(
    title = "Rho"
  )
  
  axz <- list(
    title = "Proportion"
  )
  
  ###TYPE 3###
  ##AIC UN##
  fig1 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_UN) 
  fig1 <- fig1 %>% add_surface()
  fig1 <- fig1 %>% layout(title = "Type 3 AIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC SIM##
  fig2 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_SIM) 
  fig2 <- fig2 %>% add_surface()
  fig2 <- fig2 %>% layout(title = "Type 3 AIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC CS##
  fig3 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_CS) 
  fig3 <- fig3 %>% add_surface()
  fig3 <- fig3 %>% layout(title = "Type 3 AIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC AR1##
  fig4 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_AR1) 
  fig4 <- fig4 %>% add_surface()
  fig4 <- fig4 %>% layout(title = "Type 3 AIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC CSH##
  fig5 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_CSH) 
  fig5 <- fig5 %>% add_surface()
  fig5 <- fig5 %>% layout(title = "Type 3 AIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC ARH1##
  fig6 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_ARH1) 
  fig6 <- fig6 %>% add_surface()
  fig6 <- fig6 %>% layout(title = "Type 3 AIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 4###
  ##AIC UN##
  fig7 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_UN) 
  fig7 <- fig7 %>% add_surface()
  fig7 <- fig7 %>% layout(title = "Type 4 AIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC SIM##
  fig8 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_SIM) 
  fig8 <- fig8 %>% add_surface()
  fig8 <- fig8 %>% layout(title = "Type 4 AIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC CS##
  fig9 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_CS) 
  fig9 <- fig9 %>% add_surface()
  fig9 <- fig9 %>% layout(title = "Type 4 AIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC AR1##
  fig10 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_AR1) 
  fig10 <- fig10 %>% add_surface()
  fig10 <- fig10 %>% layout(title = "Type 4 AIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC CSH##
  fig11 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_CSH) 
  fig11 <- fig11 %>% add_surface()
  fig11 <- fig11 %>% layout(title = "Type 4 AIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC ARH1##
  fig12 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_ARH1) 
  fig12 <- fig12 %>% add_surface()
  fig12 <- fig12 %>% layout(title = "Type 4 AIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 3###
  ##AICc UN##
  fig13 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_UN) 
  fig13 <- fig13 %>% add_surface()
  fig13 <- fig13 %>% layout(title = "Type 3 AICc UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc SIM##
  fig14 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_SIM) 
  fig14 <- fig14 %>% add_surface()
  fig14 <- fig14 %>% layout(title = "Type 3 AICc SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CS##
  fig15 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_CS) 
  fig15 <- fig15 %>% add_surface()
  fig15 <- fig15 %>% layout(title = "Type 3 AICc CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc AR1##
  fig16 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_AR1) 
  fig16 <- fig16 %>% add_surface()
  fig16 <- fig16 %>% layout(title = "Type 3 AICc AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CSH##
  fig17 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_CSH) 
  fig17 <- fig17 %>% add_surface()
  fig17 <- fig17 %>% layout(title = "Type 3 AICc CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc ARH1##
  fig18 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_ARH1) 
  fig18 <- fig18 %>% add_surface()
  fig18 <- fig18 %>% layout(title = "Type 3 AICc ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 4###
  ##AICc UN##
  fig19 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_UN) 
  fig19 <- fig19 %>% add_surface()
  fig19 <- fig19 %>% layout(title = "Type 4 AICc UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc SIM##
  fig20 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_SIM) 
  fig20 <- fig20 %>% add_surface()
  fig20 <- fig20 %>% layout(title = "Type 4 AICc SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CS##
  fig21 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_CS) 
  fig21 <- fig21 %>% add_surface()
  fig21 <- fig21 %>% layout(title = "Type 4 AICc CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc AR1##
  fig22 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_AR1) 
  fig22 <- fig22 %>% add_surface()
  fig22 <- fig22 %>% layout(title = "Type 4 AICc AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CSH##
  fig23 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_CSH) 
  fig23 <- fig23 %>% add_surface()
  fig23 <- fig23 %>% layout(title = "Type 4 AICc CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc ARH1##
  fig24 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_ARH1) 
  fig24 <- fig24 %>% add_surface()
  fig24 <- fig24 %>% layout(title = "Type 4 AICc ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 3###
  ##BIC UN##
  fig25 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_UN) 
  fig25 <- fig25 %>% add_surface()
  fig25 <- fig25 %>% layout(title = "Type 3 BIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC SIM##
  fig26 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_SIM) 
  fig26 <- fig26 %>% add_surface()
  fig26 <- fig26 %>% layout(title = "Type 3 BIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CS##
  fig27 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_CS) 
  fig27 <- fig27 %>% add_surface()
  fig27 <- fig27 %>% layout(title = "Type 3 BIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC AR1##
  fig28 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_AR1) 
  fig28 <- fig28 %>% add_surface()
  fig28 <- fig28 %>% layout(title = "Type 3 BIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CSH##
  fig29 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_CSH) 
  fig29 <- fig29 %>% add_surface()
  fig29 <- fig29 %>% layout(title = "Type 3 BIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC ARH1##
  fig30 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_ARH1) 
  fig30 <- fig30 %>% add_surface()
  fig30 <- fig30 %>% layout(title = "Type 3 BIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 4###
  ##BIC UN##
  fig31 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_UN) 
  fig31 <- fig31 %>% add_surface()
  fig31 <- fig31 %>% layout(title = "Type 4 BIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC SIM##
  fig32 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_SIM) 
  fig32 <- fig32 %>% add_surface()
  fig32 <- fig32 %>% layout(title = "Type 4 BIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CS##
  fig33 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_CS) 
  fig33 <- fig33 %>% add_surface()
  fig33 <- fig33 %>% layout(title = "Type 4 BIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC AR1##
  fig34 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_AR1) 
  fig34 <- fig34 %>% add_surface()
  fig34 <- fig34 %>% layout(title = "Type 4 BIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CSH##
  fig35 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_CSH) 
  fig35 <- fig35 %>% add_surface()
  fig35 <- fig35 %>% layout(title = "Type 4 BIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC ARH1##
  fig36 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_ARH1) 
  fig36 <- fig36 %>% add_surface()
  fig36 <- fig36 %>% layout(title = "Type 4 BIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  #uncomment if want widget links for viewing and sharing outside R
  htmlwidgets::saveWidget(as_widget(fig1), "type3_AIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig2), "type3_AIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig3), "type3_AIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig4), "type3_AIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig5), "type3_AIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig6), "type3_AIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig7), "type4_AIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig8), "type4_AIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig9), "type4_AIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig10), "type4_AIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig11), "type4_AIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig12), "type4_AIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig13), "type3_AICc_UN.html")
  htmlwidgets::saveWidget(as_widget(fig14), "type3_AICc_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig15), "type3_AICc_CS.html")
  htmlwidgets::saveWidget(as_widget(fig16), "type3_AICc_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig17), "type3_AICc_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig18), "type3_AICc_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig19), "type4_AICc_UN.html")
  htmlwidgets::saveWidget(as_widget(fig20), "type4_AICc_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig21), "type4_AICc_CS.html")
  htmlwidgets::saveWidget(as_widget(fig22), "type4_AICc_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig23), "type4_AICc_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig24), "type4_AICc_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig25), "type3_BIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig26), "type3_BIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig27), "type3_BIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig28), "type3_BIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig29), "type3_BIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig30), "type3_BIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig31), "type4_BIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig32), "type4_BIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig33), "type4_BIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig34), "type4_BIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig35), "type4_BIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig36), "type4_BIC_ARH1.html")
  
  #uncomment if want to see in R instead of widget
  # return(
  # list(
  #   fig1,
  #   fig2,
  #   fig3,
  #   fig4,
  #   fig5,
  #   fig6,
  #   fig7,
  #   fig8,
  #   fig9,
  #   fig10,
  #   fig11,
  #   fig12,
  #   fig13,
  #   fig14,
  #   fig15,
  #   fig16,
  #   fig17,
  #   fig18,
  #   fig19,
  #   fig20,
  #   fig21,
  #   fig22,
  #   fig23,
  #   fig24,
  #   fig25,
  #   fig26,
  #   fig28,
  #   fig29,
  #   fig30,
  #   fig31,
  #   fig32,
  #   fig33,
  #   fig34,
  #   fig35,
  #   fig36
  # )
  # )
}








######temp#####
for (d in data){
  dat = data_retrieve(d)
  thumb_plot34(dat, d, 13)
}

#
# Description
#
# @param threshold: could be .98, .9, .8. 
# @param dataList: List of names of CSV files
#
thumb_cutoffs <- function (threshold, dataList, exp_col_num_AIC) {

  #rules of thumb to go off of
  thumb_vect = seq(0, 7, .5)
  
  #retrieve all data
  for (data in dataList) {
    #get each data set in list of file names
    #d = data_retrieve(data)
    thumb_vect = seq(0, 7, .5)
    
    ###rows of these matrices will go top to bottom p = .1, .2, ... , .5, .8
    type3_data = matrix(0, length(thumb_vect), 18)
    colnames(type3_data) = c(
      "AIC_UN",
      "AIC_SIM",
      "AIC_CS",
      "AIC_AR1",
      "AIC_CSH",
      "AIC_ARH1",
      "AICc_UN",
      "AICc_SIM",
      "AICc_CS",
      "AICc_AR1",
      "AICc_CSH",
      "AICc_ARH1",
      "BIC_UN",
      "BIC_SIM",
      "BIC_CS",
      "BIC_AR1",
      "BIC_CSH",
      "BIC_ARH1"
    )
    
    #matrix where columns are ICs and models (18 cols)
    #and rows are thumb value (.5, 1, 1.5, ...)
    type4_data = matrix(0, length(thumb_vect), 18)
    colnames(type4_data) = colnames(type3_data)
    
    #
   
    #get diff matrix
    
    for (p in c(.1, .2, .3, .4, .5, .8)){
      diff = diff(data, exp_col_num_AIC)
      
      row_count = 1
      for (thumb in thumb_vect){
        quad = quad_correct(diff, thumb)
        
        ##
        #type 3
        ##
        
        #AIC
        type3_data[row_count, 1] = quad[3, 1]
        type3_data[row_count, 2] = quad[3, 4]
        type3_data[row_count, 3] = quad[3, 7]
        type3_data[row_count, 4] = quad[3, 10]
        type3_data[row_count, 5] = quad[3, 13]
        type3_data[row_count, 6] = quad[3, 16]
        
        #AICc
        type3_data[row_count, 7] = quad[3, 2]
        type3_data[row_count, 8] = quad[3, 5]
        type3_data[row_count, 9] = quad[3, 8]
        type3_data[row_count, 10] = quad[3, 11]
        type3_data[row_count, 11] = quad[3, 14]
        type3_data[row_count, 12] = quad[3, 17]
        
        #BIC
        type3_data[row_count, 13] = quad[3, 3]
        type3_data[row_count, 14] = quad[3, 6]
        type3_data[row_count, 15] = quad[3, 9]
        type3_data[row_count, 16] = quad[3, 12]
        type3_data[row_count, 17] = quad[3, 15]
        type3_data[row_count, 18] = quad[3, 18]
        
        ##
        #type 4
        ##
        
        #AIC
        type4_data[row_count, 1] = quad[4, 1]
        type4_data[row_count, 2] = quad[4, 4]
        type4_data[row_count, 3] = quad[4, 7]
        type4_data[row_count, 4] = quad[4, 10]
        type4_data[row_count, 5] = quad[4, 13]
        type4_data[row_count, 6] = quad[4, 16]
        
        #AICc
        type4_data[row_count, 7] = quad[4, 2]
        type4_data[row_count, 8] = quad[4, 5]
        type4_data[row_count, 9] = quad[4, 8]
        type4_data[row_count, 10] = quad[4, 11]
        type4_data[row_count, 11] = quad[4, 14]
        type4_data[row_count, 12] = quad[4, 17]
        
        #BIC
        type4_data[row_count, 13] = quad[4, 3]
        type4_data[row_count, 14] = quad[4, 6]
        type4_data[row_count, 15] = quad[4, 9]
        type4_data[row_count, 16] = quad[4, 12]
        type4_data[row_count, 17] = quad[4, 15]
        type4_data[row_count, 18] = quad[4, 18]
        
        
        row_count = row_count + 1
      }
    }
    print(type3_data)
    print(type4_data)
    }
    
}


