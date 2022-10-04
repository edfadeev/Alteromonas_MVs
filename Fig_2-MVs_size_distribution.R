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
path.NTA <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/Manuscripts/Alteromonas_vesicle_production/Data/Fig_2-MVs_size_distribution/"

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

