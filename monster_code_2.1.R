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
  
  CSfit_AICc = summary(CSfit)$AIC + (2 * (2)) * (2 + 1) / (n_sub - 2 -
                                                             1)
  AR1fit_AICc = summary(AR1fit)$AIC + (2 * (2)) * (2 + 1) / (n_sub - 2 -
                                                               1)
  UNfit_AICc = summary(UNfit)$AIC + (2 * (n_obs + choose(n_obs, 2))) * (n_obs + choose(n_obs, 2) + 1) / (n_sub - (n_obs + choose(n_obs, 2)) -
                                                                                                           1)
  CSHfit_AICc = summary(CSHfit)$AIC + (2 * (n_obs + 1)) * ((n_obs + 1)  + 1) / (n_sub - ( n_obs + 1) -
                                                                                  1)
  ARH1fit_AICc = summary(ARH1fit)$AIC + (2 * (n_obs + 1)) * ((n_obs + 1)  + 1) / (n_sub - (n_obs + 1) -
                                                                                    1)
  SIMfit_AICc = summary(SIMfit)$AIC + (2 * (1)) * (1 + 1) / (n_sub - 1 -
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
  ##same as original job_results_gen but fills in csv one line at a time
  ##to save completed trials that were done before a failed trial
  ##the way it is written now it just overwrites the previous csv with new
  ##matrix with a new row rather than appending, may need to change this later!!
  
  means_string = "means"
  for (mean in means){
    means_string = paste(means_string, mean, sep = "_")
  }
  
  sigma_string = "sigma"
  for (sigma in sigma_vect){
    sigma_string = paste(sigma_string, sigma, sep = "_")
  }
  
  results = matrix(0, N, 18)
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
  
  for (i in 1:N){
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
    
    results[i,] = res
    
    write.csv(results, file_name)
  }
  
  return(file_name)
}
####Data Retrieval Process Streamlined####
data_retrieve <- function(file_name){
  ##note need to enter file name in quotes e.g. data_retrieve("file_name")
  data = read.csv(file_name, header = TRUE)
  data = as.matrix(data)
  data = data[,2:19]
  return(data)
}

####~~~~~~~~~~~~~~#Data Analysis Functions#~~~~~~~~~~~~~~~~~~####
####Data Analysis####
#outputs 1-6 ranked least to greatest for each IC
data_analysis <- function(results){
  dimension = dim(results)
  N = dimension[1]
  count_IC = dimension[2]
  
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
quad_correct <- function(diff_matrix, thumb){
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
    
    quad_matrix[1, i] = count_1/N
    quad_matrix[2, i] = count_2/N
    quad_matrix[3, i] = count_3/N
    quad_matrix[4, i] = count_4/N
    
  }
  
  return(quad_matrix)
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
      
      hist1 = hist(data[,exp_col_num_AIC], breaks = breakpoints, plot = FALSE)
      hist2 = hist(data[,col], breaks = breakpoints, plot = FALSE)
      plot(hist1, xlab = paste("Blue:", names[exp_col_num_AIC], "Pink:", names[col], sep = " "), main = paste(dataset, "Visualization of AIC Comparison", sep = ": "), col = c1)
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
      
      hist1 = hist(data[,exp_col_num_AICc], breaks = breakpoints, plot = FALSE)
      hist2 = hist(data[,col], breaks = breakpoints, plot = FALSE)
      plot(hist1, xlab = paste("Blue:", names[exp_col_num_AICc], "Pink:", names[col], sep = " "), main = paste(dataset, "Visualization of AICc Comparison", sep = ": "), col = c1)
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
      
      hist1 = hist(data[,exp_col_num_BIC], breaks = breakpoints, plot = FALSE)
      hist2 = hist(data[,col], breaks = breakpoints, plot = FALSE)
      plot(hist1, xlab = paste("Blue:", names[exp_col_num_BIC], "Pink:", names[col], sep = " "), main = paste(dataset, "Visualization of BIC Comparison", sep = ": "), col = c1)
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
  layout(mat = matrix(c(1, 3, 5, 2, 4, 6), nrow = 3, ncol = 2))
  
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
    main = paste("Percentage of AIC Type 3 Results Varying", x_vect_var_name, sep = " "),
    xlab = x_vect_var_name,
    ylab = "Percentage",
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
    main = paste("Percentage of AIC Type 4 Results Varying", x_vect_var_name, sep = " "),
    xlab = x_vect_var_name,
    ylab = "Percentage",
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
    main = paste("Percentage of AICc Type 3 Results Varying", x_vect_var_name, sep = " "),
    xlab = x_vect_var_name,
    ylab = "Percentage",
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
    main = paste("Percentage of AICc Type 4 Results Varying", x_vect_var_name, sep = " "),
    xlab = x_vect_var_name,
    ylab = "Percentage",
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
    main = paste("Percentage of BIC Type 3 Results Varying", x_vect_var_name, sep = " "),
    xlab = x_vect_var_name,
    ylab = "Percentage",
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
    main = paste("Percentage of BIC Type 4 Results Varying Rho", x_vect_var_name, sep = " "),
    xlab = x_vect_var_name,
    ylab = "Percentage",
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

