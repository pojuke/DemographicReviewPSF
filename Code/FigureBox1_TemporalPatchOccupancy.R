########################################################################################################################
########################################################################################################################
#### Figure Box 1 in Ke et al. -- Time will tell: the temporal and demographic contexts of plant--soil microbe interactions
#### Demostrate how a patch occupancy model can be used to study the temporal dimension of PSF
#### Parameter values are obtained from Teste et al. (2017) Science
########################################################################################################################
########################################################################################################################



########################################################################################################################
########################################################################################################################
#### Install the required R packages
########################################################################################################################
########################################################################################################################
#### load the package into active memory
require('deSolve')
require('rootSolve')
require('tidyverse')
require('viridis')
require('scales')
require('gridExtra')
require('cowplot')
library('magick')
library('patchwork')
library('readxl')



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
#### The patch occupancy model based on Ke & Levine (2021) Am Nat
########################################################################################################################
########################################################################################################################
#### The model in matrix form
DecayPSF <- function(Time, State, Pars){
  
  with(as.list(c(State, Pars)),{
    
    #### Extract parameters in vector and matrix form
    r <- Pars$r 
    m <- Pars$m 
    d <- Pars$d 
    sigma <- Pars$sigma 
    n <- length(Pars$r) 
    
    #### Extract state variables from the full-length vector
    p00 <- State[1]
    pii <- State[2:(n+1)]
    p0i <- State[(n+2):(2*n+1)]
    
    #### Calculate derivatives using matrix notation
    dp00 <- - p00 * (r %*% pii) + d %*% p0i
    dpii <- r * pii * p00 - m * pii + r * pii * (sigma %*% p0i)
    dp0i <- m * pii - d * p0i - p0i * (t(sigma) %*% (r * pii))
    return(list(c(dp00, dpii, dp0i)))
    
  })
}    


#### Simulation function that includes numerical details: deSolve (time series)
SimulateODE <- function(parameters, initial){
  times <- seq(0, 10000, by = 1)
  state <- initial
  parms <- parameters
  pop_size <- ode(func = DecayPSF, times = times, y = state, parms = parms)
  return(pop_size)
}


#### Simulation function that includes numerical details: rootSolve (final equilibrium)
SimulateODE_Endpoint <- function(parameters, initial){
  # times <- seq(0, 6000, by = 1)
  state <- initial
  parms <- parameters
  pop_size <- rootSolve::runsteady(y = state, time = c(0, Inf), func = DecayPSF, parms = parms, 
                                   mf = 10, stol = 1e-8, ynames = TRUE)$y
  return(pop_size)
}



########################################################################################################################
########################################################################################################################
#### Set parameters
########################################################################################################################
########################################################################################################################
#### Sixteen species from Teste et al. (2017) Science 
#### Read in raw data
RawData <- as.data.frame(read_excel("../Data/Teste_et_al_2017_Science_Biomass.xlsx", sheet = 1, col_names = T))


#### The raw data includes the growth response of 16 species to soils of 7 root strategies (10 replicates each)
#### Calculate the average total biomass of each species to each soil 
Data_Average <- 
  RawData %>% 
  mutate(RootWt_g = as.numeric(RootWt_g),
         ShootWt_g = as.numeric(ShootWt_g)) %>% 
  filter(!is.na(RootWt_g)) %>% 
  filter(!is.na(ShootWt_g)) %>% 
  mutate(TotalWt_g = ShootWt_g + RootWt_g) %>% 
  group_by(Species, InoculaType, Strategy) %>%
  summarise(Biomass = mean(TotalWt_g)) 


#### As the data only includes species response to soil of different root strategies
#### The \sigma_{ij} matrix is based on species i's response to the soil conditioned by the root strategy of plant j 
#### Create the vector of species and their root strategies
Species_name <- c("Vd", "Xp", "Cq", "Cs", "Ea", "Ep", "Et", "Ahu", "Apu", "Gc", "Jf", "Vj", "Ba", "Bm", "Hi", "Hr")
Species_strategy <- c("AMF", "AMF", "EMF", "EMF", "EMF", "EMF", "EMF", "NFIX", "NFIX", "NFIX", "NFIX", "NFIX", "PMIN", "PMIN", "PMIN", "PMIN")
  

#### Fill in matrix based on species i's response to plant j's strategy's soil
#### Note that the diagonal is different since it's the conspecific soil response
Teste <- matrix(0, 16, 16)
for(i in 1:16){
  for(j in 1:16){
    Temp <- Data_Average[Data_Average$Species == Species_name[i], ]
    if(i == j){
      Teste[i, j] <- Temp[Temp$InoculaType == "Conspecific", ]$Biomass
    } else{
      Teste[i, j] <- Temp[Temp$InoculaType == Species_strategy[j], ]$Biomass
    }
  }
}


#### Standardized based on each species maximum growth rate
Teste_shifted <- t(apply(Teste, 1, function(x){x / max(x)}))


#### Aggregate parameters into vectors
set.seed(12345) 
N <- 16
r.vec <- runif(N, 0.2, 0.25)
m.vec <- rep(0.05, N)
d.vec <- rep(0.1, N)
sigma.mat <- Teste_shifted
init <- c(p00 = 0.2, 
          pAA = 0.05, pBB = 0.05, pCC = 0.05, pDD = 0.05,  pEE = 0.05, pFF = 0.05, pGG = 0.05, pHH = 0.05, pII = 0.05, pJJ = 0.05, pKK = 0.05, pLL = 0.05,  pMM = 0.05, pNN = 0.05, pOO = 0.05, pPP = 0.05,     
          p0A = 0.0, p0B = 0.0, p0C = 0.0, p0D = 0.0,  p0E = 0.0, p0F = 0.0, p0G = 0.0, p0H = 0.0, p0I = 0.0, p0J = 0.0, p0K = 0.0, p0L = 0.0,  p0M = 0.0, p0N = 0.0, p0O = 0.0, p0P = 0.0)



########################################################################################################################
########################################################################################################################
#### Start time series simulation
#### This simulation uses deSolve
########################################################################################################################
########################################################################################################################
#### Create data using functional programming
TimeSeries <-
  crossing(d.value = c(0.01, 0.99)) %>%
  mutate(r = list(r.vec), 
         m = list(m.vec),
         sigma = map(rep(1, n()), ~matrix(c(sigma.mat), N, N, byrow=T))) %>%
  mutate(d = map2(d.value, N, rep)) %>%
  mutate(Par = pmap(., list)) %>%
  mutate(Initial = list(init)) %>%
  mutate(Scenario = as.character(d.value)) %>%
  mutate(Sim = map2(Par, Initial, SimulateODE)) %>%
  mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
  select(Result) %>%
  unnest(Result) %>%
  rename(Scenario = "...1") %>%
  select(1:(2+length(init)-N)) %>%
  gather(key = Variable, value = Density, -c(Scenario, time))


#### Plot data
TimeSeriesPlot <-
  TimeSeries %>%
  filter(Variable != "p00") %>%
  ggplot(aes(x = time, 
             y = Density, 
             color = Variable)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  scale_x_continuous(name = "Time") +  
  scale_y_continuous(name = "Species abundance", 
                     limits = c(0, 1),
                     expand = c(0.005, 0)) + 
  geom_vline(xintercept = -Inf, size = 1) + theme(axis.line.y = element_blank()) + 
  facet_wrap(~ Scenario, 
             label = as_labeller(c("0.01" = "Slow decay", 
                                   "0.99" = "Fast decay"))) +   
  guides(color = guide_legend("Species")) + 
  theme(legend.position = "none")



########################################################################################################################
########################################################################################################################
#### Start multiple simulations with randomly-drawn fecundity
#### This simulation uses rootSolve
########################################################################################################################
########################################################################################################################
#### Create saving space for multiple simulations using functional programming
#### This simulation will take some time. Thus, you can also directly load a previously saved RDS file
#### Data <- readRDS("../Data/DecaySimulation.rds")
Run <- 100
Data <- as.data.frame(matrix(0, N * 2, Run))


#### Generate multiple simulations using functional programming
#### Skip this code-chunk if you choose to import the previously saved RDS file
for(i in 1:Run){
  
  ## Randomly-draw fecundity vector
  r.vec.sim = runif(N, 0.2, 0.25)
  
  ## One simulation
  Temp <-
    crossing(d.value = c(0.01, 0.99)) %>%
    mutate(r = list(r.vec.sim), 
           m = list(m.vec),
           sigma = map(rep(1, n()), ~matrix(c(sigma.mat), N, N, byrow=T))) %>%
    mutate(d = map2(d.value, N, rep)) %>%
    mutate(Par = pmap(., list)) %>%
    mutate(Initial = list(init)) %>%
    mutate(Scenario = as.character(d.value)) %>%
    mutate(Sim = map2(Par, Initial, SimulateODE_Endpoint)) %>%
    mutate(Result = map2(Scenario, Sim, bind_cols)) %>%
    select(Result) %>%
    mutate(Result = map(Result, ~.x[which(lapply(str_split(names(init), "p0"), length) == 1), ])) %>%
    unnest(Result) %>%
    rename(Scenario = "...1", Density = "...2") %>%
    select(Density)
  
  ## Save into data
  Data[, i] <- Temp
}


#### Plot and visualize results -- persistence frequency
#### Persistence frequency barplot 
MultipleSimulationPlot <- 
  Data %>% 
  mutate(Species = rep(LETTERS[1:N], 2)) %>%
  mutate(Scenario = rep(c("0.01", "0.99"), each = N)) %>%
  mutate(Mean = rowSums(Data > 10^(-5))) %>%
  select(c("Species", "Scenario", "Mean")) %>%
  ggplot(aes(x = Species, y = Mean, color = Species, fill = Species)) + 
  geom_bar(stat = "identity", alpha = 0.6, color="black", show_guide=FALSE) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -Inf, size = 1) + theme(axis.line.y = element_blank()) + 
  scale_x_discrete(name = "Species") +  
  scale_y_continuous(name = "# of times persisting",
                     limits = c(0, Run), 
                     expand = c(0.005, 0)) + 
  facet_wrap(~ Scenario, 
             label = as_labeller(c("0.01" = "Slow decay", 
                                   "0.99" = "Fast decay"))) + 
  theme(legend.position = "none")

  
#### Use patchwork to combine two GGPLOT objects
dev.off()
pdf(file="FigureBox1_PatchOccupancyDecay.pdf", width = 7, height = 5.5)
TimeSeriesPlot / MultipleSimulationPlot                       
dev.off()



# ########################################################################################################################
# ########################################################################################################################
# #### R-script to combine simulation with conceptual figure from BioRender
# #### Just to make Figure Box 1 fully reproducible
# ########################################################################################################################
# ########################################################################################################################
# #### Read conceptual figure PDF into R and convert to raster PNG
# Conceptual <- image_read_pdf("../Data/ModelFramework_PatchOccupancyModel.pdf") %>%
#   grid::rasterGrob()
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
# pdf(file="FigureBox1_PatchOccupancyDecay.pdf", width = 7, height = 9)
# Doodle / TimeSeriesPlot / MultipleSimulationPlot + plot_layout(design = "\n#AAAAAAAAA#\n#BBBBBBBBB#\n#CCCCCCCCC#", 
#                                                                heights = c(2, 1, 1)) + 
#   plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") + 
#   theme(plot.tag = element_text(face="bold"))
# dev.off()


