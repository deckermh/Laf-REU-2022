###########Sigma Generation Functions###########
#Hi
##
#Make a CSH matrix
##
makeCSH <- function(n, sigmas, rho) {
  if (n != length(sigmas)) {
    return(NULL)
  }
  
  data = array(dim = c(n, n))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        data[i, j] = sigmas[i] * sigmas[j]
      }
      else {
        data[i, j] = rho * sigmas[i] * sigmas[j]
      }
    }
  }
  result = matrix(data,
                  nrow = n,
                  ncol = n,
                  byrow = T)
  
  return(result)
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
makeARH1 <-
  function(n, p, sigma) {
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

#######################################################

##
# Data Generation Function
##

generate_data <-
  function(n_obs, n_sub, Sigma, means) {
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

#######################################################

###
#Fit Data Function
###

fit_data <- function(clean_data, n_obs, n_sub) {
  observation = clean_data$observation
  time = clean_data$time
  id = clean_data$id
  
  time = as.factor(time)
  
  
  library(nlme)
  
  ##
  #Unstructured
  ##
  UNfit = gls(observation ~ time, corr = corSymm(form = ~ 1 |id),weights = varIdent(form = ~ 1 | time))
  summary(UNfit)
  
  ##
  #Simple?
  ##
  
  ##
  #Compound Symmetry
  ##
  CSfit = gls(observation ~ time, corr = corCompSymm(form = ~ 1 |id))
  summary(CSfit)
  
  ##
  #AR(1)
  ##
  AR1fit = gls(observation ~ time, corr = corAR1(form = ~ 1 | id))
  summary(AR1fit)
  
  ##
  #CSH
  ##
  CSHfit = gls(observation ~ time, corr = corCompSymm(form = ~ 1 |id),weights = varIdent(form = ~ 1 | time))
  
  ##
  #ARH(1)
  ##
  ARH1fit = gls(observation ~ time,corr = corAR1(form = ~ 1 |id),weights = varIdent(form = ~ 1 | time))
  
  CSfit_AIC = summary(CSfit)$AIC
  AR1fit_AIC = summary(AR1fit)$AIC
  UNfit_AIC = summary(UNfit)$AIC
  CSHfit_AIC = summary(CSHfit)$AIC
  ARH1fit_AIC = summary(ARH1fit)$AIC
  
  CSfit_BIC = summary(CSfit)$BIC
  AR1fit_BIC = summary(AR1fit)$BIC
  UNfit_BIC = summary(UNfit)$BIC
  CSHfit_BIC = summary(CSHfit)$BIC
  ARH1fit_BIC = summary(ARH1fit)$BIC
  
  CSfit_AICc = summary(CSfit)$AIC + (2 * (n_obs)) * (n_obs + 1) / (n_sub - n_obs -
                                                                     1)
  AR1fit_AICc = summary(AR1fit)$AIC + (2 * (n_obs)) * (n_obs + 1) / (n_sub - n_obs -
                                                                       1)
  UNfit_AICc = summary(UNfit)$AIC + (2 * (n_obs)) * (n_obs + 1) / (n_sub - n_obs -
                                                                     1)
  CSHfit_AICc = summary(CSHfit)$AIC + (2 * (n_obs)) * (n_obs + 1) / (n_sub - n_obs -
                                                                       1)
  ARH1fit_AICc = summary(ARH1fit)$AIC + (2 * (n_obs)) * (n_obs + 1) / (n_sub - n_obs -
                                                                         1)
  
  return (
    c(
      UNfit_AIC,
      UNfit_AICc,
      UNfit_BIC,
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

#######################################################

###
#Full Monster Code###
###


results_matrix <- function(N, n_obs, n_sub, Sigma, means) {
  ## inputs:
  ## N = # of trials to be run
  ## Sigma = generated sigma matrix/matrices
  
  ICvectorlength = 15
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

#Generates and fits data, varying rho each trial
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
