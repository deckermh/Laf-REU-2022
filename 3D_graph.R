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
  
  type3_AICc_UN = type3_AIC
  type3_AICc_SIM = type3_AIC
  type3_AICc_CS = type3_AIC
  type3_AICc_AR1 = type3_AIC
  type3_AICc_CSH = type3_AIC
  type3_AICc_ARH1 = type3_AIC
  
  type4_AICc_UN = type3_AIC
  type4_AICc_SIM = type3_AIC
  type4_AICc_CS = type3_AIC
  type4_AICc_AR1 = type3_AIC
  type4_AICc_CSH = type3_AIC
  type4_AICc_ARH1 = type3_AIC
  
  type3_BIC_UN = type3_AIC
  type3_BIC_SIM = type3_AIC
  type3_BIC_CS = type3_AIC
  type3_BIC_AR1 = type3_AIC
  type3_BIC_CSH = type3_AIC
  type3_BIC_ARH1 = type3_AIC
  
  type4_BIC_UN = type3_AIC
  type4_BIC_SIM = type3_AIC
  type4_BIC_CS = type3_AIC
  type4_BIC_AR1 = type3_AIC
  type4_BIC_CSH = type3_AIC
  type4_BIC_ARH1 = type3_AIC
  
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
    
    type3_AICc_UN[col_index[i], row_index[i]] = quad[3, 2]
    type3_AICc_SIM[col_index[i], row_index[i]] = quad[3, 5]
    type3_AICc_CS[col_index[i], row_index[i]] = quad[3, 8]
    type3_AICc_AR1[col_index[i], row_index[i]] = quad[3, 11]
    type3_AICc_CSH[col_index[i], row_index[i]] = quad[3, 14]
    type3_AICc_ARH1[col_index[i], row_index[i]] = quad[3, 17]
    
    type4_AICc_UN[col_index[i], row_index[i]] = quad[4, 2]
    type4_AICc_SIM[col_index[i], row_index[i]] = quad[4, 5]
    type4_AICc_CS[col_index[i], row_index[i]] = quad[4, 8]
    type4_AICc_AR1[col_index[i], row_index[i]] = quad[4, 11]
    type4_AICc_CSH[col_index[i], row_index[i]] = quad[4, 14]
    type4_AICc_ARH1[col_index[i], row_index[i]] = quad[4, 17]
    
    type3_BIC_UN[col_index[i], row_index[i]] = quad[3, 3]
    type3_BIC_SIM[col_index[i], row_index[i]] = quad[3, 6]
    type3_BIC_CS[col_index[i], row_index[i]] = quad[3, 9]
    type3_BIC_AR1[col_index[i], row_index[i]] = quad[3, 12]
    type3_BIC_CSH[col_index[i], row_index[i]] = quad[3, 15]
    type3_BIC_ARH1[col_index[i], row_index[i]] = quad[3, 18]
    
    type4_BIC_UN[col_index[i], row_index[i]] = quad[4, 3]
    type4_BIC_SIM[col_index[i], row_index[i]] = quad[4, 6]
    type4_BIC_CS[col_index[i], row_index[i]] = quad[4, 9]
    type4_BIC_AR1[col_index[i], row_index[i]] = quad[4, 12]
    type4_BIC_CSH[col_index[i], row_index[i]] = quad[4, 15]
    type4_BIC_ARH1[col_index[i], row_index[i]] = quad[4, 18]
  }
  
  axx <- list(
    title = "Heteroskedasticity"
  )
  
  axy <- list(
    title = "Rho"
  )
  
  axz <- list(
    title = "Proportion"
  )
  
  ###TYPE 3###
  ##AIC UN##
  fig1 <- plot_ly(x = c(.8, .5, .4, .3, .2, .1), y = c(.1, .2, .3, .4, .5, .8), z = type3_AIC_UN) 
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
  
  ###TYPE 3###
  ##AICc UN##
  fig13 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_UN) 
  fig13 <- fig13 %>% add_surface()
  fig13 <- fig13 %>% layout(title = "Type 3 AICc UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc SIM##
  fig14 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_SIM) 
  fig14 <- fig14 %>% add_surface()
  fig14 <- fig14 %>% layout(title = "Type 3 AICc SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CS##
  fig15 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_CS) 
  fig15 <- fig15 %>% add_surface()
  fig15 <- fig15 %>% layout(title = "Type 3 AICc CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc AR1##
  fig16 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_AR1) 
  fig16 <- fig16 %>% add_surface()
  fig16 <- fig16 %>% layout(title = "Type 3 AICc AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CSH##
  fig17 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_CSH) 
  fig17 <- fig17 %>% add_surface()
  fig17 <- fig17 %>% layout(title = "Type 3 AICc CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc ARH1##
  fig18 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_AICc_ARH1) 
  fig18 <- fig18 %>% add_surface()
  fig18 <- fig18 %>% layout(title = "Type 3 AICc ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 4###
  ##AICc UN##
  fig19 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_UN) 
  fig19 <- fig19 %>% add_surface()
  fig19 <- fig19 %>% layout(title = "Type 4 AICc UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc SIM##
  fig20 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_SIM) 
  fig20 <- fig20 %>% add_surface()
  fig20 <- fig20 %>% layout(title = "Type 4 AICc SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CS##
  fig21 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_CS) 
  fig21 <- fig21 %>% add_surface()
  fig21 <- fig21 %>% layout(title = "Type 4 AICc CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc AR1##
  fig22 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_AR1) 
  fig22 <- fig22 %>% add_surface()
  fig22 <- fig22 %>% layout(title = "Type 4 AICc AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc CSH##
  fig23 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_CSH) 
  fig23 <- fig23 %>% add_surface()
  fig23 <- fig23 %>% layout(title = "Type 4 AICc CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##AICc ARH1##
  fig24 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_AICc_ARH1) 
  fig24 <- fig24 %>% add_surface()
  fig24 <- fig24 %>% layout(title = "Type 4 AICc ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 3###
  ##BIC UN##
  fig25 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_UN) 
  fig25 <- fig25 %>% add_surface()
  fig25 <- fig25 %>% layout(title = "Type 3 BIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC SIM##
  fig26 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_SIM) 
  fig26 <- fig26 %>% add_surface()
  fig26 <- fig26 %>% layout(title = "Type 3 BIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CS##
  fig27 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_CS) 
  fig27 <- fig27 %>% add_surface()
  fig27 <- fig27 %>% layout(title = "Type 3 BIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC AR1##
  fig28 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_AR1) 
  fig28 <- fig28 %>% add_surface()
  fig28 <- fig28 %>% layout(title = "Type 3 BIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CSH##
  fig29 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_CSH) 
  fig29 <- fig29 %>% add_surface()
  fig29 <- fig29 %>% layout(title = "Type 3 BIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC ARH1##
  fig30 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type3_BIC_ARH1) 
  fig30 <- fig30 %>% add_surface()
  fig30 <- fig30 %>% layout(title = "Type 3 BIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ###TYPE 4###
  ##BIC UN##
  fig31 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_UN) 
  fig31 <- fig31 %>% add_surface()
  fig31 <- fig31 %>% layout(title = "Type 4 BIC UN", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC SIM##
  fig32 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_SIM) 
  fig32 <- fig32 %>% add_surface()
  fig32 <- fig32 %>% layout(title = "Type 4 BIC SIM", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CS##
  fig33 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_CS) 
  fig33 <- fig33 %>% add_surface()
  fig33 <- fig33 %>% layout(title = "Type 4 BIC CS", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC AR1##
  fig34 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_AR1) 
  fig34 <- fig34 %>% add_surface()
  fig34 <- fig34 %>% layout(title = "Type 4 BIC AR1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC CSH##
  fig35 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_CSH) 
  fig35 <- fig35 %>% add_surface()
  fig35 <- fig35 %>% layout(title = "Type 4 BIC CSH", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  ##BIC ARH1##
  fig36 <- plot_ly(x = c(.1, .2, .3, .4, .5, .8), y = c(.1, .2, .3, .4, .5, .8), z = type4_BIC_ARH1) 
  fig36 <- fig36 %>% add_surface()
  fig36 <- fig36 %>% layout(title = "Type 4 BIC ARH1", scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  #uncomment if want widget links for viewing and sharing outside R
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
  
  htmlwidgets::saveWidget(as_widget(fig13), "type3_AIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig14), "type3_AIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig15), "type3_AIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig16), "type3_AIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig17), "type3_AIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig18), "type3_AIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig19), "type4_AIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig20), "type4_AIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig21), "type4_AIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig22), "type4_AIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig23), "type4_AIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig24), "type4_AIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig25), "type3_BIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig26), "type3_BIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig27), "type3_BIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig28), "type3_BIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig29), "type3_BIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig30), "type3_BIC_ARH1.html")
  
  htmlwidgets::saveWidget(as_widget(fig31), "type4_BIC_UN.html")
  htmlwidgets::saveWidget(as_widget(fig32), "type4_BIC_SIM.html")
  htmlwidgets::saveWidget(as_widget(fig33), "type4_BIC_CS.html")
  htmlwidgets::saveWidget(as_widget(fig34), "type4_BIC_AR1.html")
  htmlwidgets::saveWidget(as_widget(fig35), "type4_BIC_CSH.html")
  htmlwidgets::saveWidget(as_widget(fig36), "type4_BIC_ARH1.html")

  #uncomment if want to see in R instead of widget
  return(
  list(
    fig1,
    fig2,
    fig3,
    fig4,
    fig5,
    fig6,
    fig7,
    fig8,
    fig9,
    fig10,
    fig11,
    fig12,
    fig13,
    fig14,
    fig15,
    fig16,
    fig17,
    fig18,
    fig19,
    fig20,
    fig21,
    fig22,
    fig23,
    fig24,
    fig25,
    fig26,
    fig28,
    fig29,
    fig30,
    fig31,
    fig32,
    fig33,
    fig34,
    fig35,
    fig36
  )
  )
}