data_analysis <- function(results){
  dimension = dim(results)
  N = dimension[1]
  count_IC = dimension[2]
  analysis_matrix = matrix(0, N, 3)
  colnames(analysis_matrix) = c("minAIC", "minAICc", "minBIC")
  rownames(analysis_matrix) = c(1:N)
  
  for (i in 1:N){
    
    ###AIC###
    
    all_AIC = c(results[i, seq(count_IC, from=1, by=3)])
    minAIC = min(all_AIC)
    name = names(all_AIC)
    
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
    final_input = paste(min_type, minAIC, sep = ", ")
    analysis_matrix[i,1] = final_input
    
    ###AICc###
    
    all_AICc = c(results[i, seq(count_IC, from=2, by=3)])
    minAICc = min(all_AICc)
    name = names(all_AICc)
    
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
    final_input = paste(min_type, minAICc, sep = ", ")
    analysis_matrix[i,2] = final_input
    
    ###BIC###
    
    all_BIC = c(results[i, seq(count_IC, from=3, by=3)])
    minBIC = min(all_BIC)
    name = names(all_BIC)
    
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
    final_input = paste(min_type, minBIC, sep = ", ")
    analysis_matrix[i,3] = final_input
  }
      
  return(analysis_matrix)
}
  