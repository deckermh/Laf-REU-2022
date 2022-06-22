###Function to know if the IC chose correctly

#thumb = variable to use for "rule of thumb"

# 1 = Success
# 0 = Fail

correct <- function(exp_col_num_AIC, diff_matrix, thumb) {
  N = dim(diff_matrix)[1]
  
  exp_col_num_AICc = exp_col_num_AIC + 1
  exp_col_num_BIC = exp_col_num_AIC + 2
  
  for (i in 1:N) {
    ##Failures and successes for AIC
    for (j in seq(18, from = 1, by = 3)) {
      if (j == exp_col_num_AIC) {
        
      }
      else{
        if (diff_matrix[i, j] > 0 &&
            diff_matrix[i, j] < thumb) {
          ##Right model, not significant
          print(1)
        }
        else if (diff_matrix[i, j] < 0 &&
                 diff_matrix[i, j] > -(thumb)) {
          ##Wrong model, not significant
          print(2)
        }
        else if (diff_matrix[i, j] < 0 &&
                 diff_matrix[i, j] < -(thumb)) {
          ##Wrong model, Significant
          print(3)
        }
        else if (diff_matrix[i, j] > 0 &&
                 diff_matrix[i, j] > thumb) {
          ##Right model, significant
          print(4)
        }
        else {
          ##difference of 0
          print(5)
        }
      }
    }
    
    #Failure and successes for AICc
    
    for (j in seq(18, from = 2, by = 3)) {
      if (j == exp_col_num_AICc) {
        
      }
      else{
        if (diff_matrix[i, j] > 0 &&
            diff_matrix[i, j] < thumb) {
          ##Right model, not significant
          print(1)
        }
        else if (diff_matrix[i, j] < 0 &&
                 diff_matrix[i, j] > -(thumb)) {
          ##Wrong model, not significant
          print(2)
        }
        else if (diff_matrix[i, j] < 0 &&
                 diff_matrix[i, j] < -(thumb)) {
          ##Wrong model, Significant
          print(3)
        }
        else if (diff_matrix[i, j] > 0 &&
                 diff_matrix[i, j] > thumb) {
          ##Right model, signifigant
          print(4)
        }
      }
    }
    
    #Failure and successes for BIC
    for (j in seq(18, from = 3, by = 3)) {
      if (j == exp_col_num_BIC) {
        
      }
      if (diff_matrix[i, j] > 0 &&
          diff_matrix[i, j] < thumb) {
        ##Right model, not significant
        print(1)
      }
      else if (diff_matrix[i, j] < 0 &&
               diff_matrix[i, j] > -(thumb)) {
        ##Wrong model, not significant
        print(2)
      }
      else if (diff_matrix[i, j] < 0 &&
               diff_matrix[i, j] < -(thumb)) {
        ##Wrong model, Significant
        print(3)
      }
      else if (diff_matrix[i, j] > 0 &&
               diff_matrix[i, j] > thumb) {
        ##Right model, significant
        print(4)
      }
    }
  }
}
