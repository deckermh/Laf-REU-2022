graph_3D <- function(base_data_name, AIC_exp_num, thumb){
  ###!!!!!!!important if this is ur first time using
  ###make sure to un-comment the install packages and library lines
  
  ##base_data_name looks like this
  ##"N_2500_obs_5_sub_50_ARH1_sigma_1_10_10_10_10_p_0.1_means_0_0_0_0_0.csv"
  ##as long as you have all 36 necessary data files on ur desktop this func
  ##will retrieve all neccesary data just from this base name
  ##rn this only works for hetero models if obs is same and are just 
  ##varying rho and heteroskedasticity ratio
 
  #install.packages("plotly")
  library(plotly)
  #install.packages("stringr")
  library("stringr")
  
  type3_AIC = matrix(0, 6, 6)
  colnames(type3_AIC) = c("p .1" , "p .2", "p .3", "p .4", "p .5", "p .8")
  rownames(type3_AIC) = c("ratio .1" , "ratio .2", "ratio .3", "ratio .4", "ratio .5", "ratio .8")
  
  type3_AIC_UN = type3_AIC
  type3_AIC_SIM = type3_AIC
  type3_AIC_CS = type3_AIC
  type3_AIC_AR1 = type3_AIC
  type3_AIC_CSH = type3_AIC
  type3_AIC_ARH1 = type3_AIC
  
  type4_AIC_UN = type3_AIC
  type4_AIC_SIM = type3_AIC
  type4_AIC_CS = type3_AIC
  type4_AIC_AR1 = type3_AIC
  type4_AIC_CSH = type3_AIC
  type4_AIC_ARH1 = type3_AIC
  
  #change below to be like above if get to that point
  # type3_AICc = matrix(0, 6, 6)
  # colnames(type3_AICc) = colnames(type3_AIC)
  # rownames(type3_AICc) = rownames(type3_AIC)
  # 
  # type4_AICc = matrix(0, 6, 6)
  # colnames(type4_AICc) = colnames(type3_AIC)
  # rownames(type4_AICc) = rownames(type3_AIC)
  # 
  # type3_BIC = matrix(0, 6, 6)
  # colnames(type3_BIC) = colnames(type3_AIC)
  # rownames(type3_BIC) = rownames(type3_AIC)
  # 
  # type4_BIC = matrix(0, 6, 6)
  # colnames(type4_BIC) = colnames(type3_AIC)
  # rownames(type4_BIC) = rownames(type3_AIC)
  
  ##data list fills in reading style left to right top to bottom
  ## eg p1_het1 ... p8_het1, p1_het2, ... 
  
  str_pieces = str_split(base_data_name, "sigma_1", n = 3, simplify = TRUE)
  temp = str_split(str_pieces[2], "p_0.1", simplify = TRUE)
  str_pieces[2] = temp[1]
  str_pieces[3] = temp[2]
  
  #data list will be list of 36 filenames to be called
  data_list = c()
  for (ratio in c(1:5, 8)){
    for (p in c(1:5, 8)){
      filename = paste(str_pieces[1], "sigma_", ratio, str_pieces[2], "p_0.", p, str_pieces[3], sep = "")
      data_list = c(data_list, filename)
    }
  }
  
  ##this is ugly sorry y'all
  col_index = rep(1:6, 6)
  row_index = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 6), rep(6, 6))
  for (i in 1:36){
    file = data_list[i]
    data = data_retrieve(file)
    
    diff = diff(data, AIC_exp_num)
    quad = quad_correct(diff, thumb)
    # print(quad)
    
    type3_AIC_UN[col_index[i], row_index[i]] = quad[3, 1]
    type3_AIC_SIM[col_index[i], row_index[i]] = quad[3, 4]
    type3_AIC_CS[col_index[i], row_index[i]] = quad[3, 7]
    type3_AIC_AR1[col_index[i], row_index[i]] = quad[3, 10]
    type3_AIC_CSH[col_index[i], row_index[i]] = quad[3, 13]
    type3_AIC_ARH1[col_index[i], row_index[i]] = quad[3, 16]

    type4_AIC_UN[col_index[i], row_index[i]] = quad[4, 1]
    type4_AIC_SIM[col_index[i], row_index[i]] = quad[4, 4]
    type4_AIC_CS[col_index[i], row_index[i]] = quad[4, 7]
    type4_AIC_AR1[col_index[i], row_index[i]] = quad[4, 10]
    type4_AIC_CSH[col_index[i], row_index[i]] = quad[4, 13]
    type4_AIC_ARH1[col_index[i], row_index[i]] = quad[4, 16]
    
  }
  
  axx <- list(
    title = "Rho"
  )
  
  axy <- list(
    title = "Heteroskedasticity"
  )
  
  axz <- list(
    title = "Proportion"
  )
  
  ###TYPE 3###
  ##AIC UN##
  fig1 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_UN) 
  fig1 <- fig1 %>% add_surface()
  fig1 <- fig1 %>% layout(title = "Type 3 AIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

  ##AIC SIM##
  fig2 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_SIM) 
  fig2 <- fig2 %>% add_surface()
  fig2 <- fig2 %>% layout(title = "Type 3 AIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

  ##AIC CS##
  fig3 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_CS) 
  fig3 <- fig3 %>% add_surface()
  fig3 <- fig3 %>% layout(title = "Type 3 AIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

  ##AIC AR1##
  fig4 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_AR1) 
  fig4 <- fig4 %>% add_surface()
  fig4 <- fig4 %>% layout(title = "Type 3 AIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

  ##AIC CSH##
  fig5 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_CSH) 
  fig5 <- fig5 %>% add_surface()
  fig5 <- fig5 %>% layout(title = "Type 3 AIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

  ##AIC ARH1##
  fig6 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_ARH1) 
  fig6 <- fig6 %>% add_surface()
  fig6 <- fig6 %>% layout(title = "Type 3 AIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 4###
  ##AIC UN##
  fig7 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_UN) 
  fig7 <- fig7 %>% add_surface()
  fig7 <- fig7 %>% layout(title = "Type 4 AIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC SIM##
  fig8 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_SIM) 
  fig8 <- fig8 %>% add_surface()
  fig8 <- fig8 %>% layout(title = "Type 4 AIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC CS##
  fig9 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_CS) 
  fig9 <- fig9 %>% add_surface()
  fig9 <- fig9 %>% layout(title = "Type 4 AIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC AR1##
  fig10 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_AR1) 
  fig10 <- fig10 %>% add_surface()
  fig10 <- fig10 %>% layout(title = "Type 4 AIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC CSH##
  fig11 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_CSH) 
  fig11 <- fig11 %>% add_surface()
  fig11 <- fig11 %>% layout(title = "Type 4 AIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AIC ARH1##
  fig12 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AIC_ARH1) 
  fig12 <- fig12 %>% add_surface()
  fig12 <- fig12 %>% layout(title = "Type 4 AIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  
    
  htmlwidgets::saveWidget(as_widget(fig1), "type3_AIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig2), "type3_AIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig3), "type3_AIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig4), "type3_AIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig5), "type3_AIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig6), "type3_AIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig7), "type4_AIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig8), "type4_AIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig9), "type4_AIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig10), "type4_AIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig11), "type4_AIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig12), "type4_AIC_ARH1.html")

  
  #return(list(fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12))
}