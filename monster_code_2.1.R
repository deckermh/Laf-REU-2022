####Sigma Generation Functions####

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




####Data Generation Function####
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


####Fit Data Function####
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
  SIMfit_AIC = summary(SIMfit)$AIC
  
  CSfit_BIC = summary(CSfit)$BIC
  AR1fit_BIC = summary(AR1fit)$BIC
  UNfit_BIC = summary(UNfit)$BIC
  CSHfit_BIC = summary(CSHfit)$BIC
  ARH1fit_BIC = summary(ARH1fit)$BIC
  SIMfit_BIC = summary(SIMfit)$BIC
  
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
  SIMfit_AICc = summary(SIMfit)$AIC + (2 * (n_obs)) * (n_obs + 1) / (n_sub - n_obs -
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



####Data Analysis####

data_analysis <- function(results){
  dimension = dim(results)
  N = dimension[1]
  count_IC = dimension[2]
  analysis_matrix = matrix(0, N, 4)
  colnames(analysis_matrix) = c("minAIC", "minAICc", "minBIC", "AIC-BIC")
  rownames(analysis_matrix) = c(1:N)
  
  sorted_results = matrix(0, N, 18)
  colnames(sorted_results) = c(paste("AIC", 1:6, sep = "_"), paste("AICc", 1:6, sep = "_"), paste("BIC", 1:6, sep = "_"))
  
  
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
    
    model_index = 1
    for (term in all_AIC){
      if (term != minAIC){
        model_index = model_index + 1
      }
      else{
        break
      }
    }
    min_type = name[model_index]
    min_type = substr(min_type, 1, nchar(min_type)-7)
    final_input = paste(min_type, minAIC, sep = ", ")
    analysis_matrix[i,1] = final_input
    
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
    
    model_index = 1
    for (term in all_AICc){
      if (term != minAICc){
        model_index = model_index + 1
      }
      else{
        break
      }
    }
    min_type = name[model_index]
    min_type = substr(min_type, 1, nchar(min_type)-8)
    final_input = paste(min_type, minAICc, sep = ", ")
    analysis_matrix[i,2] = final_input
    
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
    
    model_index = 1
    for (term in all_BIC){
      if (term != minBIC){
        model_index = model_index + 1
      }
      else{
        break
      }
    }
    min_type = name[model_index]
    min_type = substr(min_type, 1, nchar(min_type)-7)
    final_input = paste(min_type, minBIC, sep = ", ")
    analysis_matrix[i,3] = final_input
    
    ##Difference between AIC and BIC##
    diff = abs(minAIC - minBIC)
    analysis_matrix[i,4] = diff
  }
  
  analysis_matrix = cbind(sorted_results, analysis_matrix)
  
  return(analysis_matrix)
}




####Full Monster Code####

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
  results = data_analysis(results)
  return(results)
}

