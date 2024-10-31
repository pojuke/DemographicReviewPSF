########################################################################################################################
########################################################################################################################
#### Figure 2 in Ke et al. -- Time will tell: the temporal and demographic contexts of plant--soil microbe interactions
#### Code to simulate different annual plant models to motivate the study  
########################################################################################################################
########################################################################################################################



########################################################################################################################
########################################################################################################################
#### Install the required R packages
########################################################################################################################
########################################################################################################################
#### load the package into active memory
require('tidyverse')
library('cowplot')
library('patchwork')
library("readxl")
library("scatterpie")
library("GGally")
library("ggExtra")
library("grid")
library("magick")



########################################################################################################################
########################################################################################################################
#### Figure theme setting
########################################################################################################################
########################################################################################################################
theme_plots <- function() {
  theme_classic() +
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(),
          strip.text.x = element_text(hjust = 0),
          axis.title = element_text(size = 13), 
          strip.text = element_text(size = 12),
          plot.subtitle = element_text(size = 12),
          plot.title = element_text(size = 13), 
          legend.text.align = 0)
}
theme_set(theme_plots())



########################################################################################################################
########################################################################################################################
#### The Beverton-Holt annual plant model 
########################################################################################################################
########################################################################################################################
#### The model with CONSTANT parameters
BevertonHolt <- function(parameters, frame){
  
  Nsim = parameters$N
  Gsim = parameters$g
  Ssim = parameters$s
  Lsim = parameters$lambda
  Asim = parameters$A
  
  for(step in 2:dim(frame)[1]){
    x = frame[step-1, 2:(Nsim+1)]
    frame[step, 2:(Nsim+1)] = (1 - Gsim) * Ssim * x + Gsim * x * (Lsim / (1 + Asim %*% (Gsim * x)))
  }
  return(frame)
  
}


#### The model with TIME-VARYING parameters 
#### The time-varying parameter is a big list, with elements also bing a list for parameters
BevertonHolt_Temporal <- function(parameters.megalist, frame){
  
  for(step in 2:dim(frame)[1]){
    
    Nsim = parameters.megalist[[step-1]]$N
    Gsim = parameters.megalist[[step-1]]$g
    Ssim = parameters.megalist[[step-1]]$s
    Lsim = parameters.megalist[[step-1]]$lambda
    Asim = parameters.megalist[[step-1]]$A
    
    x = frame[step-1, 2:(Nsim+1)]
    frame[step, 2:(Nsim+1)] = (1 - Gsim) * Ssim * x + Gsim * x * (Lsim / (1 + Asim %*% (Gsim * x)))
  }
  return(frame)
  
}



########################################################################################################################
########################################################################################################################
#### Other house-keeping functions/objects/parameters
########################################################################################################################
########################################################################################################################
#### Function to create a data.frame for the temporally-varying parameter
Create_Replace <- function(parameters.megalist, initial.value, slope, extend.time){
  total.time = length(parameters.megalist) - 1
  temp = data.frame(Time = 0:total.time) %>%
    mutate(Value = initial.value + slope * Time) %>%
    mutate(Value = ifelse(Time <= extend.time, Value, Value[Time == extend.time])) 
  return(temp)
}


#### A specific function to make LAMBDA (fecundity) a temporally-varying parameter
Create_ParameterTrend_lambda <- function(parameters.megalist, initial.value, slope, extend.time){
  Replace = Create_Replace(parameters.megalist, initial.value, slope, extend.time)
  for(step in 2:dim(Replace)[1]){
    parameters.megalist[[step]] = parameters.megalist[[step]] %>%
      list_modify(lambda = c(Replace$Value[step], 
                             unname(unlist(parameters.megalist[[1]]["lambda"]))[2]))
  }
  return(parameters.megalist)
}


#### A specific function to make S (seed survival) a temporally-varying parameter
Create_ParameterTrend_s <- function(parameters.megalist, initial.value, slope, extend.time){
  Replace = Create_Replace(parameters.megalist, initial.value, slope, extend.time)
  for(step in 2:dim(Replace)[1]){
    parameters.megalist[[step]] = parameters.megalist[[step]] %>%
      list_modify(s = c(Replace$Value[step], 
                        unname(unlist(parameters.megalist[[1]]["s"]))[2]))
  }
  return(parameters.megalist)
}


#### Default parameters
#### Obtained from the species pair Festuca microstachys (N1) vs. Hordeum murinum (N2) in Van Dyke et al. (2022)
N <- 2
g <- c(0.752, 0.667)
s <- c(0.134, 0.045)
lambda <- c(2129.950, 736.667)
A <- matrix(c(0.588, 1.411,
              0.109, 0.948), N, N, byrow = TRUE)


#### Simulation and parameter setting
Time <- 30
Init.1 <- 10
Init.2 <- 10
pchange_focal_par <- 0.1   
parms <- list(N = N, 
              g = g, 
              s = s, 
              lambda = lambda, 
              A = A)
parms.megalist <- rep(list(parms), Time+1)


#### Initial data frame setup
Frame <- matrix(c(c(0:Time), 
                 rep(c(Init.1, Init.2), each = (Time+1))), 
                nrow = Time+1, ncol = N+1, byrow = F,
                dimnames = list(as.character(0:Time), 
                                c("Time", paste("N", c(1:N), sep = ""))))



########################################################################################################################
########################################################################################################################
#### Simulate/plot different demographic impact of CONSTANT microbes
########################################################################################################################
########################################################################################################################
Data_Constant <- 
  tibble(Par = list(parms,
                    parms %>%
                      list_modify(lambda = c(lambda[1] * (1 - pchange_focal_par), lambda[2])),
                    parms %>%
                      list_modify(s = c(s[1] * (1 - pchange_focal_par), s[2]))),
         Initial = list(Frame, Frame, Frame), 
         Scenario = c("No Microbe", "Vary lambda", "Vary s")) %>% 
  mutate(Sim = map2(Par, Initial, BevertonHolt)) %>% 
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  unnest(Result) %>%
  select(!c(Par, Initial, Scenario, Sim)) %>%
  rename(Scenario = "...1") %>%
  as.data.frame() %>%
  mutate(Rel_N1 = N1 / (N1 + N2)) %>%
  gather(key = Variable, value = Density, -c(Scenario, Time)) 


MotivatePlot_Constant <- 
  Data_Constant %>% 
  filter(Variable != "Rel_N1") %>%
  ggplot(aes(x = Time, 
             y = Density, 
             color = Scenario, 
             linetype = Variable)) +
  geom_line(linewidth = 1) + 
  scale_color_manual(labels = c("No Microbe" = "No pathogenic impact",
                                "Vary lambda" = expression(paste("Pathogens decrease  ", lambda[1])), 
                                "Vary s" = expression(paste("Pathogens decrease  ", s[1]))), 
                     values = c("#999999", "#56B4E9", "#E69F00"))



########################################################################################################################
########################################################################################################################
#### Simulate/plot different demographic impact of TIME-VARYING microbes
#### Scenario: the constant microbial effect (-10%) is the beginning of conditioning and goes on for 8yr (ending at -80%)
########################################################################################################################
########################################################################################################################
Data_Temporal <- 
  tibble(Par = list(parms.megalist,
                    parms.megalist %>%
                      Create_ParameterTrend_lambda(lambda[1], -212, 8),   
                    parms.megalist %>%
                      Create_ParameterTrend_s(s[1], -0.013, 8)),          
         Initial = list(Frame, Frame, Frame), 
         Scenario = c("No Microbe", "Vary lambda", "Vary s")) %>% 
  mutate(Sim = map2(Par, Initial, BevertonHolt_Temporal)) %>% 
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  unnest(Result) %>%
  select(!c(Par, Initial, Scenario, Sim)) %>%
  rename(Scenario = "...1") %>%
  as.data.frame() %>%
  mutate(Rel_N1 = N1 / (N1 + N2)) %>%
  gather(key = Variable, value = Density, -c(Scenario, Time)) 


MotivatePlot_Temporal <- 
  Data_Temporal %>% 
  filter(Variable != "Rel_N1") %>%
  ggplot(aes(x = Time, 
             y = Density, 
             color = Scenario, 
             linetype = Variable)) +
  geom_line(linewidth = 1) + 
  scale_color_manual(labels = c("No Microbe" = "No pathogenic impact",
                                "Vary lambda" = expression(paste("Pathogens decrease  ", lambda[1])), 
                                "Vary s" = expression(paste("Pathogens decrease  ", s[1]))), 
                     values = c("#999999", "#56B4E9", "#E69F00"))



########################################################################################################################
########################################################################################################################
#### Combine temporal and demographic scenarios into the same figure 
########################################################################################################################
########################################################################################################################
MotivatePlot_Combine <-
  rbind(mutate(Data_Constant, ScenarioTemp = "Constant pathogenic effect"), 
        mutate(Data_Temporal, ScenarioTemp = "Varying pathogenic effect")) %>%
  as.data.frame() %>%
  filter(Variable != "Rel_N1") %>%
  ggplot(aes(x = Time, 
             y = Density, 
             color = Scenario, 
             linetype = Variable)) +
  geom_line(linewidth = 1.0) +
  scale_y_continuous(name = "Species abundance", 
                     limits = c(0, 3000),
                     expand = c(0.005, 0)) + 
  geom_vline(xintercept = -Inf, size = 1) + theme(axis.line.y = element_blank()) + 
  facet_wrap(~ ScenarioTemp, 
             label = as_labeller(c("Constant pathogenic effect" = "Constant pathogenic effect", 
                                   "Varying pathogenic effect" = "Varying pathogenic effect"))) + 
  scale_linetype_discrete(name = "Species", 
                          labels = c("N1" = expression(N[1]), 
                                     "N2" = expression(N[2]))) + 
  scale_color_manual(name = "Scenario",
                     labels = c("No Microbe" = "No pathogenic impact",
                                "Vary lambda" = expression(paste("Pathogens decrease  ", lambda[1])),
                                "Vary s" = expression(paste("Pathogens decrease  ", s[1]))),
                     values = c("#999999", "#E69F00", "#56B4E9")) + 
  guides(colour = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) 
  


########################################################################################################################
########################################################################################################################
#### Create simulation and plot
########################################################################################################################
########################################################################################################################
#### Save figure
dev.off()
pdf(file="Figure2_MotivatingAnnualPlant.pdf", width = 6.5, height = 4)
MotivatePlot_Combine + 
  theme(legend.justification = "left", 
        legend.position = 'bottom', 
        legend.box = "vertical", 
        legend.box.just = "center", 
        legend.margin = margin(-5, 0, 0, 0))
dev.off()



# ########################################################################################################################
# ########################################################################################################################
# #### R-script to combine simulation with conceptual figure from BioRender
# #### Just to make Figure 2 fully reproducible
# ########################################################################################################################
# ########################################################################################################################
# #### Read conceptual figure PDF into R and convert to raster PNG
# Conceptual <- image_read_pdf("../Data/ModelFramework_AnnualModel.pdf") %>%
#   grid::rasterGrob()
# 
# 
# #### Create simulation GGPLOT file
# Simulation <- 
#   MotivatePlot_Combine + 
#   theme(legend.justification = "left", 
#         legend.position = 'bottom', 
#         legend.box = "vertical", 
#         legend.box.just = "center", 
#         legend.margin = margin(-5, 0, 0, 0))
# 
# 
# #### Fake empty canvas to paste raster PNG onto a GGPLOT
# Doodle <- 
#   ggplot(tibble(x = c(0, 8))) +
#   aes(x = x, y = x) +
#   geom_point() +
#   scale_x_continuous(limits = c(0.95, 7.05), expand = c(0, 0)) + 
#   scale_y_continuous(limits = c(1.40, 4.60), expand = c(0, 0)) + 
#   ggplot2::annotation_custom(Conceptual, 1, 7, 1, 5) + 
#   theme_void() + 
#   theme(plot.tag = element_text(face="bold"))
# 
# 
# #### Use patchwork to combine two GGPLOT objects
# dev.off()
# pdf(file="Figure2_MotivatingAnnualPlant.pdf", width = 8, height = 9)
# Doodle / Simulation + 
#   plot_layout(design = "\nAAAAAAAAA\n#BBBBBBB#",
#               heights = c(2, 1)) + 
#   plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") + 
#   theme(plot.tag = element_text(face="bold"))
# dev.off()
