
#
# Plots 4 line graphs: type 1, 2, 3, and 4 rates as a lever increases.
#
# @param modelNum: a quad correct matrix has 18 columns, each corresponding to a 
# model (ex. AIC_UN). The "model" parameter is a column number 1 thru 18
# specifying which model you want to plot.
#
# @param quads: a list of quad correct matrices, each corresponding to a lever
# value (ex. 10 quad correct matrices each corresponding to rho = 0.1, 0.2, ..., 1)
#
# @param leverValues: a list of lever values corresponding to each quad correct
# matrix (see "quads" description above).
#
# @param leverName: the name of the lever on the x-axis (ex. "rho")
#
type1234_plot <- function(modelNum, quads, leverValues, leverName = "") {
  # we have a bunch of 4x3 quad correct matrices. The rows
  # are type 1 2 3 4, and the columns are each IC. I'm
  # going to take the AIC columns for each quad correct, and
  # make them into a single matrix so I can plot it.
  # then I repeat with AICc and BIC.
  AICscores = matrix(nrow = 4, ncol = length(leverValues))
  AICcscores = matrix(nrow = 4, ncol = length(leverValues))
  BICscores = matrix(nrow = 4, ncol = length(leverValues))
  
  i = 1
  # for each quad matrix in "quads", take the column of the desired model, and
  # add it to rateMatrix
  for (q in quads) {
    rateMatrix[, i] = q[, modelNum]
    
    #increment i so we can update next column of rateMatrix
    i = i + 1
  }
  
  #transpose so we can plot columns against each other
  rateMatrix = t(rateMatrix) 
  
  #plot type 1 vs type 2 vs type 3 vs type 4 as lever changes
  matplot(rateMatrix, type = c("b"), pch = 1, col = 1:4, xlab = leverName, xaxt = "n")
  axis(side = 1, at = 1:4, labels = leverValues)
  legend("topleft", legend = 1:4, col=1:4, pch=1)
}