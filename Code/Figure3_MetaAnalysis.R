########################################################################################################################
########################################################################################################################
#### Figure 3 in Ke et al. -- Time will tell: the temporal and demographic contexts of plant--soil microbe interactions
#### Code for meta-analysis based on Crawford et al. (2019) and Yan et al. (2021) data set 
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



########################################################################################################################
########################################################################################################################
#### Load data 
########################################################################################################################
########################################################################################################################
Data <- read_excel("../Data/Data_DemographicReview_CrawfordAndYan.xlsx", sheet = 1)



########################################################################################################################
########################################################################################################################
#### Create new data frame with growth form matching 
########################################################################################################################
########################################################################################################################
Data.New <- 
  Data %>%
  filter(GrowthFormA != "NA", 
         GrowthFormB != "NA") %>% 
  filter(ConditionLength != "NA", 
         ResponseLength != "NA") %>%
  mutate(GrowthFormA_Cat1 = sapply(str_split(GrowthFormA, "_"), function(x) x[1]), 
         GrowthFormA_Cat2 = sapply(str_split(GrowthFormA, "_"), function(x) x[2]), 
         GrowthFormB_Cat1 = sapply(str_split(GrowthFormB, "_"), function(x) x[1]), 
         GrowthFormB_Cat2 = sapply(str_split(GrowthFormB, "_"), function(x) x[2]),
         SameFormExact = case_when(GrowthFormA == GrowthFormB ~ GrowthFormA, 
                                   GrowthFormA != GrowthFormB ~ paste(GrowthFormA, "___", GrowthFormB, sep="")),
         SameForm_Cat1 = case_when(GrowthFormA_Cat1 == GrowthFormB_Cat1 ~ GrowthFormA_Cat1, 
                                   GrowthFormA_Cat1 != GrowthFormB_Cat1 ~ "Annual/Perennial"),
         SameForm_Cat2 = case_when(GrowthFormA_Cat2 == GrowthFormB_Cat2 ~ GrowthFormA_Cat2, 
                                   GrowthFormA_Cat2 != GrowthFormB_Cat2 ~ paste(GrowthFormA_Cat2, "___", GrowthFormB_Cat2, sep=""))) %>% 
  mutate(SameForm_Cat2 = case_when(SameForm_Cat2 %in% c("Graminoid", "Forb", "Shrub", "Tree") ~ SameForm_Cat2, 
                                   SameForm_Cat2 %in% c("Graminoid___Forb", "Forb___Graminoid") ~ "Graminoid___Forb", 
                                   SameForm_Cat2 %in% c("Graminoid___Shrub", "Shrub___Graminoid") ~ "Graminoid___Shrub", 
                                   SameForm_Cat2 %in% c("Forb___Shrub", "Shrub___Forb") ~ "Forb___Shrub"))



########################################################################################################################
########################################################################################################################
#### Plot conditioning and response length information
#### Color coded with annual vs. perennial
#### Create components: scattered pie charts + Stacked histograms as marginal distribution
########################################################################################################################
########################################################################################################################
#### Create data for plotting the scatterpie -- key is pivot_wider()
Data.ScatterPie.Cat1 <-
    Data.New %>%
    group_by(Study, ConditionLength, ResponseLength, SameForm_Cat1) %>%
    summarise(Count = n()) %>%
    mutate(Freq = round(Count / sum(Count), 3)) %>%
    mutate(NumberExperiment = sum(Count)) %>%
    select(Study, ConditionLength, ResponseLength, SameForm_Cat1, NumberExperiment, Freq) %>%
    mutate(ConditionLength = case_when(ConditionLength == "Field" ~ "40",
                                       ConditionLength != "Field" ~ ConditionLength)) %>%
    pivot_wider(names_from = SameForm_Cat1, values_from = Freq, values_fill = 0) %>%
    mutate(ConditionLength = as.numeric(ConditionLength),
           ResponseLength = as.numeric(ResponseLength),
           Radius = log10(NumberExperiment+1))


#### X axis marginal histogram
X.Margin.Cat1 <- 
  Data.New %>%
  filter(ConditionLength != 48) %>%
  select(ConditionLength, ResponseLength, SameForm_Cat1) %>% 
  mutate(ConditionLength = case_when(ConditionLength == "Field" ~ "40", 
                                     ConditionLength != "Field" ~ ConditionLength)) %>% 
  mutate(across(.cols = 1:2, .fns = as.numeric)) %>% 
  mutate_if(is.character, as.factor) %>% 
  ggplot(aes(x = ConditionLength, fill = SameForm_Cat1)) + 
  geom_histogram(color = "black", bins = 20) + 
  scale_fill_manual(values = c("Perennial" = "#009E73", 
                               "Annual" = "#E69F00", 
                               "Annual/Perennial" = "#56B4E9")) +
  scale_x_continuous(expand = c(0.04, 1), 
                     breaks = c(seq(0, 20, by = 4), 26)) +
  theme_void() +
  theme(legend.position = "none")


#### Y axis marginal histogram
Y.Margin.Cat1 <- 
  Data.New %>%
  mutate(SameForm_Cat1 = fct_recode(SameForm_Cat1, "Annual-Perennial" = "Annual/Perennial")) %>% 
  filter(ConditionLength != 48) %>%
  select(ConditionLength, ResponseLength, SameForm_Cat1) %>% 
  mutate(ConditionLength = case_when(ConditionLength == "Field" ~ "40", 
                                     ConditionLength != "Field" ~ ConditionLength)) %>% 
  mutate(across(.cols = 1:2, .fns = as.numeric)) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(ResponseLength_Jitter = jitter(ResponseLength, 20)) %>% 
  ggplot(aes(x = ResponseLength, fill = SameForm_Cat1)) + 
  geom_histogram(color = "black", bins = 15) + 
  scale_fill_manual(name = "Life history", 
                    values = c("Perennial" = "#009E73", 
                               "Annual" = "#E69F00", 
                               "Annual-Perennial" = "#56B4E9")) + 
  scale_x_continuous(expand = c(0.05, 2)) +
  theme_void() + 
  coord_flip() + 
  theme(legend.position = c(0.75, 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11.5),
        plot.margin = margin(r = 26))


########################################################################################################################
########################################################################################################################
#### Layering of scattered pie charts to create the final figure
#### This approach handles overlapping pie charts better
########################################################################################################################
########################################################################################################################
#### Create the data without pivot_wider()
jitter_amt_field <- 1.5
df <-
  Data.New %>%
  filter(ConditionLength != 48) %>%
  group_by(Study, ConditionLength, ResponseLength, SameForm_Cat1) %>%
  summarise(Count = n()) %>%
  mutate(Freq = round(Count / sum(Count), 3)) %>%
  mutate(NumberExperiment = sum(Count)) %>%
  select(Study, ConditionLength, ResponseLength, SameForm_Cat1, Count, NumberExperiment) %>%
  mutate(ConditionLength = case_when(ConditionLength == "Field" ~ "40",
                                     ConditionLength != "Field" ~ ConditionLength)) %>%
  mutate(ConditionLength = as.numeric(ConditionLength),
         ResponseLength = as.numeric(ResponseLength)) %>% 
  mutate(ConditionLength = case_when(ConditionLength == 40 ~ jitter(ConditionLength, jitter_amt_field),
                                     TRUE ~ ConditionLength))


#### Draw a pie chart with each row and save information into columns
jitter_amt <- 0
df.grobs <-
  df %>%
  filter(Study %in% unlist(Data.ScatterPie.Cat1[(Data.ScatterPie.Cat1$Perennial == 1.0 | Data.ScatterPie.Cat1$Annual == 1.0), "Study"], use.names = F)) %>%
  group_by(Study, ConditionLength, ResponseLength, NumberExperiment) %>%
  do(subplots = ggplot(., aes(1, Count, fill = SameForm_Cat1)) +
       geom_col(position = "fill", alpha = 0.5) +
       coord_polar(theta = "y") +
       scale_fill_manual(values = c("Perennial" = "#009E73",
                                    "Annual" = "#E69F00",
                                    "Annual/Perennial" = "#56B4E9")) +
       theme_void()+ guides(fill = "none")) %>%
  mutate(ResponseLength_jitter = jitter(ResponseLength, jitter_amt)) %>% 
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                           x = ConditionLength - 1.3, y = ResponseLength_jitter - 1.3,
                                           xmax = ConditionLength + 1.3, ymax = ResponseLength_jitter + 1.3)))


df.grobs.mix1 <-
  df %>%
  filter(Study %in% unlist(Data.ScatterPie.Cat1[(Data.ScatterPie.Cat1$"Annual/Perennial" == 1.0), "Study"], use.names = F)) %>%
  group_by(Study, ConditionLength, ResponseLength, NumberExperiment) %>%
  do(subplots = ggplot(., aes(1, Count, fill = SameForm_Cat1)) +
       geom_col(position = "fill", alpha = 1.0, colour = "black") +
       coord_polar(theta = "y") +
       scale_fill_manual(values = c("Perennial" = "#009E73",
                                    "Annual" = "#E69F00",
                                    "Annual/Perennial" = "#56B4E9")) +
       theme_void()+ guides(fill = "none")) %>%
  mutate(ResponseLength_jitter = jitter(ResponseLength, jitter_amt)) %>% 
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                           x = ConditionLength - 1.3, y = ResponseLength_jitter - 1.3,
                                           xmax = ConditionLength + 1.3, ymax = ResponseLength_jitter + 1.3)))


df.grobs.mix2 <-
  df %>%
  filter(Study %in% unlist(Data.ScatterPie.Cat1[(Data.ScatterPie.Cat1$"Annual/Perennial" == 1.0), "Study"], use.names = F)) %>%
  group_by(Study, ConditionLength, ResponseLength, NumberExperiment) %>%
  do(subplots = ggplot(., aes(1, Count, fill = SameForm_Cat1)) +
       geom_col(position = "fill", alpha = 1.0, colour = "transparent") +
       coord_polar(theta = "y") +
       scale_fill_manual(values = c("Perennial" = "#009E73",
                                    "Annual" = "#E69F00",
                                    "Annual/Perennial" = "#56B4E9")) +
       theme_void()+ guides(fill = "none")) %>%
  mutate(ResponseLength_jitter = jitter(ResponseLength, jitter_amt)) %>% 
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                           x = ConditionLength - 1.23, y = ResponseLength_jitter - 1.23,
                                           xmax = ConditionLength + 1.23, ymax = ResponseLength_jitter + 1.23)))


df.grobs.mix3 <-
  df %>%
  filter(Study %in% unlist(Data.ScatterPie.Cat1[!(Data.ScatterPie.Cat1$Perennial == 1.0 | Data.ScatterPie.Cat1$Annual == 1.0 | Data.ScatterPie.Cat1$"Annual/Perennial" == 1.0), "Study"], use.names = F)) %>%
  mutate(ResponseLength_Jitter = jitter(ResponseLength, 3)) %>% 
  group_by(Study, ConditionLength, ResponseLength, NumberExperiment) %>%
  do(subplots = ggplot(., aes(1, Count, fill = SameForm_Cat1)) +
       geom_col(position = "fill", alpha = 1.0, colour = "black") +
       coord_polar(theta = "y") +
       scale_fill_manual(values = c("Perennial" = "#009E73",
                                    "Annual" = "#E69F00",
                                    "Annual/Perennial" = "#56B4E9")) +
       theme_void()+ guides(fill = "none")) %>%
  mutate(ResponseLength_jitter = jitter(ResponseLength, jitter_amt)) %>% 
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                           x = ConditionLength - 1.3, y = ResponseLength_jitter - 1.3,
                                           xmax = ConditionLength + 1.3, ymax = ResponseLength_jitter + 1.3)))


#### Plot pie chart, need to adjust expand() to make plots match limits
ScatterPie.Alternative.Cat1 <-
  bind_rows(df.grobs, df.grobs.mix1, df.grobs.mix2, df.grobs.mix3) %>% 
  {ggplot(data = ., aes(ConditionLength, ResponseLength)) +
      .$subgrobs +
      geom_text(aes(label = "")) +
      scale_y_continuous(name = "Response length (month)",
                         expand = c(0.13, 1.7), 
                         breaks = seq(0, 25, by = 5),
                         labels = as.character(seq(0, 25, by = 5))) +
      scale_x_continuous(name = "Conditioning length (month)",
                         expand = c(0.05, 1.5), 
                         breaks = c(seq(0, 30, by = 5), 40.5),
                         labels = c(as.character(seq(0, 30, by = 5)), "Field")) +
      theme_classic() +
      theme(legend.position = "none",
            axis.title.x = element_text(size = 18, margin = margin(t = 10)),
            axis.title.y = element_text(size = 18, margin = margin(r = 10)),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14))}


#### Create a rectangle grob
rect <- rectGrob(
  width = unit(0.39, "in"),
  height = unit(5.5, "in"),
  gp = gpar(fill = "white", col = "white"))


#### Combine plots
Final.Figure <- 
  X.Margin.Cat1 + plot_spacer() + ScatterPie.Alternative.Cat1 + Y.Margin.Cat1 + 
  plot_layout(ncol = 2, 
              nrow = 2, 
              widths = c(4, 1),
              heights = c(1, 4))
Final.Figure <- ggdraw(Final.Figure) +
  draw_grob(rect, x = 0.12, y = -0.04)



########################################################################################################################
########################################################################################################################
#### Save the figure 
########################################################################################################################
########################################################################################################################
#### As PDF
pdf(file = "Figure3_CurrentStatus.pdf", width = 8, height = 7)
Final.Figure
dev.off()



########################################################################################################################
########################################################################################################################
#### Data summary
#### Removing plants with unknown growth form (only 10 experiments)
########################################################################################################################
########################################################################################################################
#### 1. Number of studies
Data %>%
  filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>%
  select(Study) %>% 
  unique() %>% 
  nrow()


#### 2. Number of experimental pairs
Data %>%
  filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>%
  select(`Pairwise comparison`) %>% 
  unique() %>% 
  nrow()


#### 3. Median conditioning length of studies
Data %>% 
  filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>%
  filter(ConditionLength != "Field" & !is.na(ConditionLength)) %>% 
  select(Study, ConditionLength) %>% 
  distinct() %>% 
  .$ConditionLength %>% 
  as.numeric() %>% 
  median()
  # mean()


#### 4. Number of studies with field-conditioned soil
Data %>% 
  filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>% 
  filter(ConditionLength == "Field") %>% 
  select(Study, ConditionLength) %>% 
  distinct() %>%
  .$Study %>% 
  unique() %>% 
  length()
  

#### 5. Median response length of studies
Data %>% 
  filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>% 
  filter(ResponseLength != "NA") %>% 
  select(Study, ResponseLength) %>% 
  distinct() %>% 
  .$ResponseLength %>% 
  as.numeric() %>% 
  median()
  # mean()


#### 6. Number of studies with annual–perennial pairs and those with more than one type of pairs
####    Such studies, despite focusing on only annual-annual and perennial-perennial, apply the same time frame across growth forms
NumberStudies <- 
  Data %>% 
  filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>% 
  .$Study %>% 
  unique()

Table <- data.frame(Study = NumberStudies, formtype = "Multiple")
for(i in 1:length(NumberStudies)){
  Temp <- 
    Data %>%
    filter(Study == NumberStudies[i]) %>%
    filter(GrowthFormA != "NA" & GrowthFormB != "NA") %>%
    mutate(formtype = str_c(str_extract(GrowthFormA, pattern = "[A-Za-z]*"),
                            str_extract(GrowthFormB, pattern = "[A-Za-z]*"),
                            sep = "-")) 
  Table[i, 2] <- ifelse(length(unique(Temp$formtype)) == 1, 
                        unique(Temp$formtype), 
                        Table[i, 2])
}

Table %>% 
  filter(formtype == "Annual-Perennial" | formtype == "Perennial-Annual" | formtype == "Multiple") %>%
  .$Study %>% 
  unique() %>% 
  length()
