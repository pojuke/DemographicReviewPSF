########################################################################################################################
########################################################################################################################
#### Figure Box 2 in Ke et al. -- Time will tell: the temporal and demographic contexts of plant--soil microbe interactions
#### Sensitivity analysis of the parameters in the annual-perennial plant competition model
#### Parameter values are obtained from Uricchio et al. (2019) Am Nat
########################################################################################################################
########################################################################################################################



########################################################################################################################
########################################################################################################################
#### Install the required R packages
########################################################################################################################
########################################################################################################################
#### load the package into active memory
library('tidyverse')
library("magick")
library("cowplot")
library("patchwork")



########################################################################################################################
########################################################################################################################
#### The demographic model for annual-perennial competition
########################################################################################################################
########################################################################################################################
#### Density dependence of annual plant seed production
get.fecundity.a <- function(params, n){
  with(params, unname(lambda.a*g.eff.a /
                        (1 + alpha.ap*n['A.p'] + alpha.aa*g.eff.a*n['N.a'] + a.a)))
}


#### Density dependence of perennial plant seed production
get.fecundity.p <- function(params, n){
  with(params, unname(lambda.p /
                        (1 + alpha.pp*n['A.p'] + alpha.pa*g.eff.a*n['N.a'] + a.p)))
}


#### Density dependence of perennial plant seedling survival
get.germination.p <- function(params, n){
  with(params, unname(g.eff.p*v /
                        (1 + beta.pAp*n['A.p'] + beta.pNp*g.eff.p*n['N.p'] + beta.pa*g.eff.a*n['N.a'])))
}


#### Annual-perennial competition model  
iterate.model <- function(params, n) {
  # For parameters and current model state, iterate the model one time step
  # Arguments:
  #   params -- Named list with all model parameters
  #   n -- Named vector with abundances (N.p, A.p, N.a)
  with(params, {
    N.p <- get.fecundity.p(params, n)*n['A.p'] + s.p*(1-g.eff.p)*n['N.p']
    A.p <- get.germination.p(params, n)*n['N.p'] + xi*n['A.p']
    N.a <- get.fecundity.a(params, n)*n['N.a'] + s.a*(1-g.eff.a)*n['N.a']
    return(c(N.p = unname(N.p), A.p = unname(A.p), N.a = unname(N.a)))
  })
}


#### Simulation function of the annual-perennial competition model 
DemographicModel <- function(params, initial.n, steps){
  
  #### Create result dataframe and fill in initial state
  #### This differs from the annual plant sensitivity analysis, where initial frame has to be supplied
  result <- matrix(NA, steps+1, 3)
  colnames(result) <- c('N.p', 'A.p', 'N.a')
  if(is.null(names(initial.n))){
    result[1, ] <- initial.n
  } else{
    result[1, names(initial.n)] <- initial.n
  }
  
  #### Iterate the model 
  for(i in 1:steps){
    result[i+1, ] <- iterate.model(params, result[i,])[colnames(result)]
  }
  
  #### Clean up and return results
  result <- data.frame(result) %>%
    mutate(time=row_number()-1) %>% 
    select(time, everything()) %>% 
    return()
  
}


#### Default parameters from Uricchio et al. (2019) Am Nat 
#### Plant pair: Elymus glaucus (EG; perennial) & Bromus diandrus (BD; annual)
parms <- list(
  lambda.p = 282.127, #perennial adult fecundity
  g.eff.p = 0.125, #perennial effective germination 
  s.p = 0.515, #perennial seed survival
  xi = 0.920, #perennial adult survival
  v = 0.292, #perennial seedling survival/maturation
  lambda.a = 47.594, #annual fecundity
  g.eff.a = 0.168, #annual effective germination
  s.a = 0.443, #annual seed survival
  a.a = 0.0, a.p = 0.0, #foliar pathogen effect
  alpha.aa = 0.066, alpha.ap = 0.143, alpha.pp = 0.018, alpha.pa = 0.104, #competition reduction on fecundity
  beta.pNp = 0.086, beta.pAp = 0.063, beta.pa = 0.002 #competition reduction on perennial seed survival/maturation
)


#### Simulation and parameter setting
N <- 3
Time <- 100
Init.Np <- 0.1
Init.Ap <- 0.1
Init.Na <- 0.1


#### Initial data frame setup
initial.vec <- c(N.p = Init.Np, A.p = Init.Ap, N.a = Init.Na)
Frame <- matrix(c(c(0:Time), 
                  rep(c(Init.Np, Init.Ap, Init.Na), each = (Time+1))), 
                nrow = Time+1, ncol = N+1, byrow = F,
                dimnames = list(as.character(0:Time), 
                                c("Time", paste("N", c(1:N), sep = ""))))



########################################################################################################################
########################################################################################################################
### The function for one-parameter sensitivity analysis of the annual-perennial plant model
########################################################################################################################
########################################################################################################################
Sensitivity <- function(parameters = parms, 
                        initial = initial.vec, 
                        pchange_FOCAL_par = 0.20,
                        pchange_nonfocal_par_lower = 0.05,
                        pchange_nonfocal_par_upper = 0.05,
                        n_generation = 100, 
                        n_sim = 100){

  
  #### Extract parameters and reconstruct locally within function 
  lambda.p.LOCAL <- parameters$lambda.p #perennial adult fecundity
  g.eff.p.LOCAL <- parameters$g.eff.p #perennial effective germination 
  s.p.LOCAL <- parameters$s.p #perennial seed survival
  xi.LOCAL <- parameters$xi #perennial adult survival
  v.LOCAL <- parameters$v #perennial seedling survival/maturation
  lambda.a.LOCAL <- parameters$lambda.a #annual fecundity
  g.eff.a.LOCAL <- parameters$g.eff.a #annual effective germination
  s.a.LOCAL <- parameters$s.a #annual seed survival
  a.a.LOCAL <- parameters$a.a
  a.p.LOCAL <- parameters$a.p #foliar pathogen effect
  alpha.aa.LOCAL <- parameters$alpha.aa
  alpha.ap.LOCAL <- parameters$alpha.ap 
  alpha.pp.LOCAL <- parameters$alpha.pp 
  alpha.pa.LOCAL <- parameters$alpha.pa #competition reduction on fecundity
  beta.pNp.LOCAL <- parameters$beta.pNp 
  beta.pAp.LOCAL <- parameters$beta.pAp 
  beta.pa.LOCAL <- parameters$beta.pa #competition reduction on perennial seed survival/maturation
  N <- length(initial.vec)
  Time <- n_generation
  
  
  #### Simulation with no pathogens affecting the perennial 
  out_mat <- DemographicModel(parameters, initial.vec, n_generation)
  Np_equl_default <- 1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
  # Np_equl_default <- sum(out_mat[Time, 2:N])
  
  
  #### Simulation with all perennial parameters being affected by pathogens
  parms_temp <- list(
    lambda.p = lambda.p.LOCAL * (1 - pchange_FOCAL_par), #perennial adult fecundity
    g.eff.p = g.eff.p.LOCAL * (1 - pchange_FOCAL_par), #perennial effective germination 
    s.p = s.p.LOCAL * (1 - pchange_FOCAL_par), #perennial seed survival
    xi = xi.LOCAL * (1 - pchange_FOCAL_par), #perennial adult survival
    v = v.LOCAL * (1 - pchange_FOCAL_par), #perennial seedling survival/maturation
    lambda.a = lambda.a.LOCAL, #annual fecundity
    g.eff.a  = g.eff.a.LOCAL, #annual effective germination
    s.a = s.a.LOCAL, #annual seed survival
    a.a = a.a.LOCAL, a.p = a.p.LOCAL, #foliar pathogen effect
    alpha.aa = alpha.aa.LOCAL, alpha.ap = alpha.ap.LOCAL, alpha.pp = alpha.pp.LOCAL, alpha.pa = alpha.pa.LOCAL, #competition reduction on fecundity
    beta.pNp = beta.pNp.LOCAL, beta.pAp = beta.pAp.LOCAL, beta.pa = beta.pa.LOCAL #competition reduction on perennial seed survival/maturation
  )
  out_mat <- DemographicModel(parms_temp, initial.vec, n_generation)
  Np_equl_all <- 1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
  # Np_equl_all <- sum(out_mat[Time, 2:N]) 
  
  
  #### One parameter changed at once
  #### (1) Sensitivity of g.eff.p -- perennial effective germination
  g_sens <- sapply(1:n_sim, function(x){
    parms_temp <- list(
      lambda.p = lambda.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult fecundity
      g.eff.p = g.eff.p.LOCAL * (1 - pchange_FOCAL_par), #perennial effective germination 
      s.p = s.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seed survival
      xi = xi.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult survival
      v = v.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seedling survival/maturation
      lambda.a = lambda.a.LOCAL, #annual fecundity
      g.eff.a  = g.eff.a.LOCAL, #annual effective germination
      s.a = s.a.LOCAL, #annual seed survival
      a.a = a.a.LOCAL, a.p = a.p.LOCAL, #foliar pathogen effect
      alpha.aa = alpha.aa.LOCAL, alpha.ap = alpha.ap.LOCAL, alpha.pp = alpha.pp.LOCAL, alpha.pa = alpha.pa.LOCAL, #competition reduction on fecundity
      beta.pNp = beta.pNp.LOCAL, beta.pAp = beta.pAp.LOCAL, beta.pa = beta.pa.LOCAL #competition reduction on perennial seed survival/maturation
    )
    out_mat <- DemographicModel(parms_temp, initial.vec, n_generation)
    1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
    # sum(out_mat[Time, 2:N]) 
  })
  
  
  #### One parameter changed at once
  #### (2) Sensitivity of lambda -- perennial adult fecundity
  lambda_sens <- sapply(1:n_sim, function(x){
    parms_temp <- list(
      lambda.p = lambda.p.LOCAL * (1 - pchange_FOCAL_par), #perennial adult fecundity
      g.eff.p = g.eff.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial effective germination 
      s.p = s.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seed survival
      xi = xi.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult survival
      v = v.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seedling survival/maturation
      lambda.a = lambda.a.LOCAL, #annual fecundity
      g.eff.a  = g.eff.a.LOCAL, #annual effective germination
      s.a = s.a.LOCAL, #annual seed survival
      a.a = a.a.LOCAL, a.p = a.p.LOCAL, #foliar pathogen effect
      alpha.aa = alpha.aa.LOCAL, alpha.ap = alpha.ap.LOCAL, alpha.pp = alpha.pp.LOCAL, alpha.pa = alpha.pa.LOCAL, #competition reduction on fecundity
      beta.pNp = beta.pNp.LOCAL, beta.pAp = beta.pAp.LOCAL, beta.pa = beta.pa.LOCAL #competition reduction on perennial seed survival/maturation
    )
    out_mat <- DemographicModel(parms_temp, initial.vec, n_generation)
    1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
    # sum(out_mat[Time, 2:N]) 
  })



  #### One parameter changed at once
  #### (3) Sensitivity of s -- perennial seed survival
  s_sens <- sapply(1:n_sim, function(x){
    parms_temp <- list(
      lambda.p = lambda.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult fecundity
      g.eff.p = g.eff.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial effective germination 
      s.p = s.p.LOCAL * (1 - pchange_FOCAL_par), #perennial seed survival
      xi = xi.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult survival
      v = v.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seedling survival/maturation
      lambda.a = lambda.a.LOCAL, #annual fecundity
      g.eff.a  = g.eff.a.LOCAL, #annual effective germination
      s.a = s.a.LOCAL, #annual seed survival
      a.a = a.a.LOCAL, a.p = a.p.LOCAL, #foliar pathogen effect
      alpha.aa = alpha.aa.LOCAL, alpha.ap = alpha.ap.LOCAL, alpha.pp = alpha.pp.LOCAL, alpha.pa = alpha.pa.LOCAL, #competition reduction on fecundity
      beta.pNp = beta.pNp.LOCAL, beta.pAp = beta.pAp.LOCAL, beta.pa = beta.pa.LOCAL #competition reduction on perennial seed survival/maturation
    )
    out_mat <- DemographicModel(parms_temp, initial.vec, n_generation)
    1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
    # sum(out_mat[Time, 2:N]) 
  })


  #### One parameter changed at once
  #### (4) Sensitivity of v -- perennial seedling survival/maturation
  v_sens <- sapply(1:n_sim, function(x){
    parms_temp <- list(
      lambda.p = lambda.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult fecundity
      g.eff.p = g.eff.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial effective germination 
      s.p = s.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seed survival
      xi = xi.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult survival
      v = v.LOCAL * (1 - pchange_FOCAL_par), #perennial seedling survival/maturation
      lambda.a = lambda.a.LOCAL, #annual fecundity
      g.eff.a  = g.eff.a.LOCAL, #annual effective germination
      s.a = s.a.LOCAL, #annual seed survival
      a.a = a.a.LOCAL, a.p = a.p.LOCAL, #foliar pathogen effect
      alpha.aa = alpha.aa.LOCAL, alpha.ap = alpha.ap.LOCAL, alpha.pp = alpha.pp.LOCAL, alpha.pa = alpha.pa.LOCAL, #competition reduction on fecundity
      beta.pNp = beta.pNp.LOCAL, beta.pAp = beta.pAp.LOCAL, beta.pa = beta.pa.LOCAL #competition reduction on perennial seed survival/maturation
    )
    out_mat <- DemographicModel(parms_temp, initial.vec, n_generation)
    1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
    # sum(out_mat[Time, 2:N]) 
  })
  
  
  #### One parameter changed at once
  #### (5) Sensitivity of xi -- perennial adult survival
  xi_sens <- sapply(1:n_sim, function(x){
    parms_temp <- list(
      lambda.p = lambda.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial adult fecundity
      g.eff.p = g.eff.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial effective germination 
      s.p = s.p.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seed survival
      xi = xi.LOCAL * (1 - pchange_FOCAL_par), #perennial adult survival
      v = v.LOCAL * (1 - runif(n = 1, min = -pchange_nonfocal_par_lower, max = pchange_nonfocal_par_upper)), #perennial seedling survival/maturation
      lambda.a = lambda.a.LOCAL, #annual fecundity
      g.eff.a  = g.eff.a.LOCAL, #annual effective germination
      s.a = s.a.LOCAL, #annual seed survival
      a.a = a.a.LOCAL, a.p = a.p.LOCAL, #foliar pathogen effect
      alpha.aa = alpha.aa.LOCAL, alpha.ap = alpha.ap.LOCAL, alpha.pp = alpha.pp.LOCAL, alpha.pa = alpha.pa.LOCAL, #competition reduction on fecundity
      beta.pNp = beta.pNp.LOCAL, beta.pAp = beta.pAp.LOCAL, beta.pa = beta.pa.LOCAL #competition reduction on perennial seed survival/maturation
    )
    out_mat <- DemographicModel(parms_temp, initial.vec, n_generation)
    1 - (out_mat[Time, N+1] / sum(out_mat[Time, 2:(N+1)]))
    # sum(out_mat[Time, 2:N]) 
  })
  
  
  #### Combine all the simulation results
  pars_sens_df <- data.frame(parameter = rep(c("g", "lambda", "s", "v", "xi"), each = n_sim),
                             Np_eql = c(g_sens, lambda_sens, s_sens, v_sens, xi_sens))
  return(list(Np_equl_default = Np_equl_default,
              Np_equl_all = Np_equl_all,
              pars_sens_df = pars_sens_df))
  
}



########################################################################################################################
########################################################################################################################
### The function for visualizing the simulation results
########################################################################################################################
########################################################################################################################
plot_sim_one_par_relabd <- function(dat, pchange_focal_par, pchange_nonfocal_par_lower, pchange_nonfocal_par_upper){
  
  pars_order <- 
    dat$pars_sens_df %>% 
    group_by(parameter) %>% 
    summarise(mean = mean(Np_eql)) %>% 
    mutate(No = 1:5) %>% 
    arrange(desc(mean))
  
  dat$pars_sens_df <- 
    dat$pars_sens_df %>% 
    mutate(parameter = factor(parameter, level = pars_order$parameter, ordered = T))
  
  x_labels <- c(expression(italic(g[p])),
                expression(italic(lambda[p])),
                expression(italic(s[p])), 
                expression(italic(v)), 
                expression(italic(xi)))
  
  ggplot(data = dat$pars_sens_df) + 
    geom_hline(yintercept = dat$Np_equl_default, 
               color = "grey60", linetype = "dashed", size = 1) + 
    geom_hline(yintercept = dat$Np_equl_all, 
               color = "blue", linetype = "dashed", size = 1) + 
    geom_point(aes(x = parameter, 
                   y = Np_eql), 
               position = position_jitter(width = 0.1), 
               color = "grey50", alpha = 0.50) + 
    stat_summary(aes(x = parameter, 
                     y = Np_eql), 
                 fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange", color = "red") + 
    labs(x = "Focal parameter", 
         y = expression(Relative~abundance~of~perennials), 
         title = glue::glue("Focal parameter: {pchange_focal_par}% \n Non-focal parameters: {-pchange_nonfocal_par_lower}%~{pchange_nonfocal_par_upper}%")) +
    scale_x_discrete(labels = x_labels[pars_order$No]) +
    scale_y_continuous(expand = c(0.05, 0), limits = c(-0.05, 1.05)) +
    theme_classic() + 
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 14),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.x = element_text(size = 15, margin = margin(t = 10)),
          axis.title.y = element_text(size = 15, margin = margin(r = 8)))
   
}



########################################################################################################################
########################################################################################################################
#### Create simulation and plot
########################################################################################################################
########################################################################################################################
#### Create simulation GGPLOT file
Simulation <- 
  Sensitivity(parameters = parms,
              initial = initial.vec,
              pchange_FOCAL_par = 0.20,
              pchange_nonfocal_par_lower = 0.0,
              pchange_nonfocal_par_upper = 0.05,
              n_generation = 200, 
              n_sim = 100) %>% 
  plot_sim_one_par_relabd(pchange_focal_par = 20, 
                          pchange_nonfocal_par_lower = 0, 
                          pchange_nonfocal_par_upper = 5)


#### Save figure
dev.off()
pdf(file="FigureBox2_AnnualPerennialSensitivity.pdf", width = 4, height = 3.5)
Simulation
dev.off()
