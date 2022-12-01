###########################################
#Load libraries
###########################################
library(tidyNano)
library(tidyverse)
library(ggpmisc)

#calculate standard error
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

###########################################
#Import particle data
###########################################
#define working directory with the files
path.NTA <- "Data/Fig_2-MVs_size_distribution/"

#list all Experiment files
ParticleData_files <- list.files(path=path.NTA, pattern = "-ExperimentSummary.csv",
                                 full.names = TRUE, recursive = T)


#generate a list of all experiments 
particleData_list <- lapply(ParticleData_files, function(x) {
  dat <- nanoimport(x) %>% 
    nanotidy(sep_var = c("Strain", "stamp", "Dilution")) %>% #splits the filename according to the structure
    mutate(Dilution = 1) %>% 
    separate(col = "stamp", into = c("Fraction","Date", "Time_stamp"), sep ="\\s") %>% 
    mutate(Tech_run = paste(Date, Time_stamp, sep ="_"),
           Run_n = length(unique(Tech_run)),
           Tech_run = factor(Tech_run, levels = unique(Tech_run), #generates replicates by letters according to the time stamp
                             labels = c("A","B","C","D","E")[1:Run_n]))
  
  return(dat)
})

#aggregate all the experiments into a single dataframe
tidy_data <- bind_rows(particleData_list) %>%
  mutate(Size=round(particle_size, digits = -1),#rounds the size up/down to the closest 10th
         Strain = as.character(Strain))

rm(particleData_list)

###########################################
#Bin the particle data
###########################################
#calculate mean for each particle size
tidy_data_mean <- tidy_data %>% 
  mutate(Strain = case_when(Strain =="ATCC" ~ "ATCC27126",
                            Strain =="MIT" ~ "MIT1002",
                            TRUE ~ Strain)) %>%
  group_by(Strain, particle_size) %>% 
  summarise(Count_mean = mean(Count),
            Count_sd = sd(Count)) %>% 
  mutate(Strain = factor(Strain, levels = c("AD45","ATCC27126","BGP6","BS11","HOT1A3","MIT1002")))

#plot size distribution  
tidy_data_mean %>% 
  ggplot(aes(x = particle_size, y = Count_mean)) +
  geom_errorbar(aes(ymin = Count_mean-Count_sd, ymax = Count_mean+Count_sd), colour = "gray50", alpha=0.1)+
  geom_line(size =1)+
  stat_peaks(colour = "black", span = 31, geom = "text", vjust = -2, label.fmt ="%.f",
             ignore_threshold = 0.2)+
  facet_wrap(. ~ Strain, scales = "free_y")+  #different y axis
  scale_y_continuous(labels=scales::scientific_format())+
  xlim(0, 450)+
  labs(y="Concentration (particle/ml)", x = "Size (nm)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom")
