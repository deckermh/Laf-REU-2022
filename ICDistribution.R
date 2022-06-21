AICDistribution <- function(N, n_subs, Sigma, means) {
  #matrix with AIC scores for different fits as columns
  AIC_scores = matrix(nrow = N, ncol = 6)
  
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
    hist(AIC_scores[,i])
  }
}


AICcDistribution <- function(N, n_subs, Sigma, means) {
  #matrix with AICc scores for different fits as columns
  AICc_scores = matrix(nrow = N, ncol = 6)
  
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
    hist(AICc_scores[,i])
  }
}

BICDistribution <- function(N, n_subs, Sigma, means) {
  #matrix with BIC scores for different fits as columns
  BIC_scores = matrix(nrow = N, ncol = 6)
  
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
    hist(BIC_scores[,i])
  }
}