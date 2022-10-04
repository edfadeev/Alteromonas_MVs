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
path.NTA <- "/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/Manuscripts/Alteromonas_vesicle_production/Data/Fig_3-Production_rates/"

###########################################
#Import particle data
###########################################
#list all Experiment files
ParticleData_files <- list.files(path=path.NTA, pattern = "-ExperimentSummary.csv",
                                 full.names = TRUE, recursive = T)


#generate a list of all experiments 
particleData_list <- lapply(ParticleData_files, function(x) {
  dat <- nanoimport(x) %>% 
    nanotidy(sep_var = c("Strain", "Time", "Bio_rep", "stamp", "Dilution")) %>% #splits the filename according to the structure
    mutate(Dilution = 1) #%>% 
  tidyr::separate(col = "stamp", into = c("Method","Date", "Time_stamp"), sep ="\\s")# %>% #splits the last part of the filename
  mutate(Tech_run = paste(Date, Time, sep ="_"),
         Run_n = length(unique(Tech_run)),
         Tech_run = factor(Tech_run, levels = unique(Tech_run), #generates replicates by letters according to the time stamp
                           labels = c("A","B","C","D","E")[1:Run_n]),
         True_count= Count*Dilution) # correct concentrations using dilution factor
  
  return(dat)
})

#aggregate all the experiments into a single dataframe
tidy_data <- bind_rows(particleData_list) %>%
  dplyr::filter(particle_size < 450 & particle_size > 20) %>% #filters all the particles above 250 um
  mutate_if(is.factor, as.character)%>% 
  mutate(Time = case_when(Time == "T0" ~ "0h",
                          TRUE ~ Time)) %>% 
  mutate(Strain = case_when(Strain == "AD" ~"AD45",
                            Strain == "ATCC" ~"ATCC27126",
                            Strain == "BS" ~"BS11",
                            Strain == "MIT" ~"MIT1002",
                            Strain == "1A3" ~"HOT1A3",
                            TRUE ~ Strain))

rm(particleData_list)