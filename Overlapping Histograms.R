overlap_histograms <- function(data, exp_col_num_AIC){
  ##generates 5 double histograms for each IC which compare expected type
  ##distribution to each other distribution
  
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