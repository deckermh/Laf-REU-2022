#make a plot of box plots for distribution of Type 4 success percentage in
#in IC_CS as heteroskedasticity increases
#
#@param N - number of trials
#@param type - the type of error/success to boxplot (type 1, 2, 3, or 4)
#@param maxCount - the amount of total trials to run to see when each IC fails
#@param rho - the level of correlation for the cov matrix
#
ICHeteroBoxplots <- function(N, type, maxCount = 100, rho = .7) {
  if (type != 1 && type != 2 && type != 3 && type != 4) {
    return(NULL)
  }
  
  #These matrices will store a ton of IC failure counts.
  #the first row will be N type 4 percentages at .1 heteroskedasticity
  #the second row will be at .2 heteroskedasticity, etc.
  AIC_CS_Data = matrix(nrow = 8, ncol = N)
  AICc_CS_Data = matrix(nrow = 8, ncol = N)
  BIC_CS_Data = matrix(nrow = 8, ncol = N)
  
  #we vary min/max ratio from .1 to .8
  for (i in 1:8) {
    for (j in 1:N) {
      #generate CSH matrix with varying levels of heterosked.
      sigmas = c(sqrt(i), sqrt(10), sqrt(10))
      Sigma = makeCSH(3, rho, sigmas)
      
      #get quad correct for data
      res = results_matrix(maxCount, 3, 40, Sigma, c(0, 0, 0))
      d = diff(res, 13)
      #rule of thumb is 3
      quad = quad_correct(d, 3, returnPercents = F)
      
      #get success count for IC_CS fit
      AIC_CS_count = quad[type, 7]
      AICc_CS_count = quad[type, 8]
      BIC_CS_count = quad[type, 9]
      
      #add success count to all matrices 
      AIC_CS_Data[i, j] = AIC_CS_count
      AICc_CS_Data[i, j] = AICc_CS_count
      BIC_CS_Data[i, j] = BIC_CS_count
    }
  }
  
  write.matrix(AIC_CS_Data, "./SimpleHeterosked/AIC_CS_Data_SIMPLE.csv")
  write.matrix(AICc_CS_Data, "./SimpleHeterosked/AICc_CS_Data_SIMPLE.csv")
  write.matrix(BIC_CS_Data, "./SimpleHeterosked/BIC_CS_Data_SIMPLE.csv")
  
  #create boxplots
  boxplot(
    AIC_CS_Data[1, ],
    AIC_CS_Data[2, ],
    AIC_CS_Data[3, ],
    AIC_CS_Data[4, ],
    AIC_CS_Data[5, ],
    AIC_CS_Data[6, ],
    AIC_CS_Data[7, ],
    AIC_CS_Data[8, ],
    main = paste("AIC Type ", type, "Counts Vs. Heteroskedasticity"),
    xlab = "Min/Max Ratio",
    ylab = "Type 4 Count (Out of 100)",
    names = seq(.1, .8, by = .1)
    )
  
  boxplot(
    AICc_CS_Data[1, ],
    AICc_CS_Data[2, ],
    AICc_CS_Data[3, ],
    AICc_CS_Data[4, ],
    AICc_CS_Data[5, ],
    AICc_CS_Data[6, ],
    AICc_CS_Data[7, ],
    AICc_CS_Data[8, ],
    main = paste("AICc Type ", type, "Counts Vs. Heteroskedasticity"),
    xlab = "Min/Max Ratio",
    ylab = "Type 4 Count (Out of 100)",
    names = seq(.1, .8, by = .1)
  )
  
  boxplot(
    BIC_CS_Data[1, ],
    BIC_CS_Data[2, ],
    BIC_CS_Data[3, ],
    BIC_CS_Data[4, ],
    BIC_CS_Data[5, ],
    BIC_CS_Data[6, ],
    BIC_CS_Data[7, ],
    BIC_CS_Data[8, ],
    main = paste("BIC Type ", type, "Counts Vs. Heteroskedasticity"),
    xlab = "Min/Max Ratio",
    ylab = "Type 4 Count (Out of 100)",
    names = seq(.1, .8, by = .1)
  )
}

