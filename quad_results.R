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