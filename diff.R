#diff
#
#@param exp_col_num_AIC: represents column index of AIC_[FIT]
#can be 1 (UN), 4 (SIM), 7 (CS), 10 (AR1), 13 (CSH), 16 (ARH1)
diff <- function(results, exp_col_num_AIC) {
  N = dim(results)[1]

  diff_matrix = matrix(nrow = N, ncol = 18)
  
  #for each row...
  for (i in 1:N) {
    currentRow = results[i,]
    
    exp_value_AIC = results[i, exp_col_num_AIC]
    
    #for each AIC value in results... count up by 3's
    for (j in seq(18, from=1, by = 3)) {
      #each entry in difference matrix will be a difference between exp
      #and all other covariance structure AIC values
      diff_matrix[i,j] = currentRow[j] - exp_value_AIC
    }
    
    #now get AICc column index, and value
    exp_col_num_AICc = exp_col_num_AIC + 1
    exp_value_AICc = results[i, exp_col_num_AICc]
    
    #for each AICc value...
    for (j in seq(18, from = 2, by = 3)) {
      diff_matrix[i,j] = currentRow[j] - exp_value_AICc
    }
    
    #now get BIC column index, and value
    exp_col_num_BIC = exp_col_num_AIC + 2
    exp_value_BIC = results[i, exp_col_num_BIC]
    
    #for each BIC value
    for (j in seq(18, from = 3, by = 3)) {
      diff_matrix[i,j] = currentRow[j] - exp_value_BIC
    }
  }
  
  #set col names
  colnames(diff_matrix) = colnames(results)
  
  return (diff_matrix)
}