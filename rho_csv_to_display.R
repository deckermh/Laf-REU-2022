rho1 = read.csv("Varying Rho IC scores/N_2500_obs_10_sub_100_AR1_sigma_1_p_0.1_means_0_0_0_0_0_0_0_0_0_0.csv")
rho2 = read.csv("Varying Rho IC scores/N_2500_obs_10_sub_100_AR1_sigma_1_p_0.2_means_0_0_0_0_0_0_0_0_0_0.csv")
rho3 = read.csv("Varying Rho IC scores/N_2500_obs_10_sub_100_AR1_sigma_1_p_0.3_means_0_0_0_0_0_0_0_0_0_0.csv")
rho4 = read.csv("Varying Rho IC scores/N_2500_obs_10_sub_100_AR1_sigma_1_p_0.4_means_0_0_0_0_0_0_0_0_0_0.csv")
rho5 = read.csv("Varying Rho IC scores/N_2500_obs_10_sub_100_AR1_sigma_1_p_0.5_means_0_0_0_0_0_0_0_0_0_0.csv")
rho8 = read.csv("Varying Rho IC scores/N_2500_obs_10_sub_100_AR1_sigma_1_p_0.8_means_0_0_0_0_0_0_0_0_0_0.csv")

rho1 = as.matrix(rho1)
rho2 = as.matrix(rho2)
rho3 = as.matrix(rho3)
rho4 = as.matrix(rho4)
rho5 = as.matrix(rho5)
rho8 = as.matrix(rho8)

rho1 = rho1[,2:19]
rho2 = rho2[,2:19]
rho3 = rho3[,2:19]
rho4 = rho4[,2:19]
rho5 = rho5[,2:19]
rho8 = rho8[,2:19]

rho1 = diff(rho1, 10)
rho2 = diff(rho2, 10)
rho3 = diff(rho3, 10)
rho4 = diff(rho4, 10)
rho5 = diff(rho5, 10)
rho8 = diff(rho8, 10)


q1 = new_quad_correct(rho1, 0)
q2 = new_quad_correct(rho2, 0)
q3 = new_quad_correct(rho3, 0)
q4 = new_quad_correct(rho4, 0)
q5 = new_quad_correct(rho5, 0)
q8 = new_quad_correct(rho8, 0)

quads = list(q1, q2, q3, q4, q5, q8)

type1234_plot(quads, c(.1, .2, .3, .4, .5, .8), "rho")