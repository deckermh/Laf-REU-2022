#make a plot of box plots for distribution of Type 4 success percentage in
#in IC_CS as heteroskedasticity increases
#
#@param N - number of trials
#@param type - the type of error/success to boxplot (type 1, 2, 3, or 4)
#@param IC - the type of IC to plot (0 = AIC, 1 = AICc, 2 = BIC)
#@param rho - the level of correlation for the cov matrix
#
ICHeteroBoxplots <- function(N, type, IC, rho = .7) {
  if (type != 1 && type != 2 && type != 3 && type != 4) {
    return(NULL)
  }
  if (IC != 0 && IC != 1 && IC != 2) {
    return(NULL)
  }
  
  #This matrix will store a ton of IC failure counts.
  #the first row will be N type 4 percentages at .1 heteroskedasticity
  #the second row will be at .2 heteroskedasticity, etc.
  IC_CS_Data = matrix(nrow = 8, ncol = N)
  
  #we vary min/max ratio from .1 to .8
  for (i in 1:8) {
    for (j in 1:N) {
      #generate CSH matrix with varying levels of heterosked.
      sigmas = c(sqrt(i), sqrt(10), sqrt(10))
      Sigma = makeCSH(3, rho, sigmas)
      
      #get quad correct for data (NOTE N SHOULD BE 100)
      res = results_matrix(100, 3, 40, Sigma, c(0, 0, 0))
      d = diff(res, 13)
      #rule of thumb is 3
      quad = quad_correct(d, 3, returnPercents = F)
      
      #get failure count for IC_CS fit
      IC_CS_count = quad[type, 7 + IC]
      IC_CS_Data[i, j] = IC_CS_count
    }
  }
  
  #create boxplots
  boxplot(
    IC_CS_Data[1, ],
    IC_CS_Data[2, ],
    IC_CS_Data[3, ],
    IC_CS_Data[4, ],
    IC_CS_Data[5, ],
    IC_CS_Data[6, ],
    IC_CS_Data[7, ],
    IC_CS_Data[8, ],
    main = paste("Type ", type, "Percentage for Varying Heteroskedasticity")
  )
}