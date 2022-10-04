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

#define working directory with the files
path.NTA <- "Data/Fig_3-Production_rates/"

###########################################
#Import cell counts
###########################################
cell_counts_gc<- read.csv(paste0(path.NTA,"cell_counts_summary.csv"), row.names = 1) %>% 
  select(Strain, Time, Bio_rep, Concentration) %>% 
  mutate(Type = "Cells",
         Time = as.numeric(Time)) %>% 
  mutate(Strain = case_when(Strain == "AD" ~"AD45",
                            Strain == "ATCC" ~"ATCC27126",
                            Strain == "BS" ~"BS11",
                            Strain == "MIT" ~"MIT1002",
                            Strain == "1A3" ~"HOT1A3",
                            TRUE ~ Strain)) 

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
    mutate(Dilution = 1) %>% 
  tidyr::separate(col = "stamp", into = c("Method","Date", "Time_stamp"), sep ="\\s") %>% #splits the last part of the filename
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
######################################################
### Plot production of MVs over time by replicate  ###
######################################################
#Sum number of MVs per sample 
Total_MVs_bio_rep<- tidy_data %>% 
                      group_by(particle_size, Time, Bio_rep, Strain) %>% 
                      summarise(Tech_mean = mean(True_count)) %>% 
                      group_by(Time, Bio_rep, Strain) %>% 
                      summarize(Concentration = sum(Tech_mean)) %>% 
                      mutate(Time = as.numeric(gsub("h","", Time)),
                             Type ="MVs",
                             Concentration = case_when(Time == 0 ~ Concentration*0.002, # the starter samples were diluted 100 times
                                                              TRUE ~ Concentration))


#merge cell and MVs counts and calculate mean of biological replicates
MVs_CC_merged <- Total_MVs_bio_rep %>%
  rbind(cell_counts_gc) %>% 
  group_by(Time,Strain, Type) %>% 
  summarize(Concentration.mean = mean(Concentration),
            Concentration.se = se(Concentration))

#plot cell concentration and MVs
concentrations_plot <- MVs_CC_merged %>% 
                        mutate(Time = as.numeric(as.character(Time))) %>% 
                        ggplot(aes(x=Time, y= Concentration.mean, group = interaction(Strain, Type), colour = Strain))+
                        geom_line(aes(linetype = Type), size = 1, alpha = 0.3)+
                        geom_point(aes(shape = Type),size = 3.5, colour = "black")+
                        geom_point(aes(shape = Type),size = 3)+
                        labs(y="Concentration (#/ml)", x = "Time (hours)")+
                        geom_errorbar(aes(ymin = Concentration.mean-Concentration.se, ymax = Concentration.mean+Concentration.se),  width = 2)+
                        facet_wrap(.~Strain, nrow  = 3, scales = "free_y")+
                        scale_y_log10()+
                        theme_bw() + 
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                              text=element_text(size=14),legend.position = "bottom")

##########################################
### Calculate MVs production rates     ###
##########################################
#copy the T0 values for each replicate
NTA_CC_merged <- read.csv("NTA_CC_merged_EF.csv", sep = ";")


#merge cell and MVs counts by biological replicates
MVs_CC_bio_rep <- Total_MVs_bio_rep %>%
  rbind(cell_counts_gc) %>% 
  group_by(Time,Strain, Bio_rep, Type) %>% 
  summarize(Concentration.mean = mean(Concentration))

#define the time points
timePoints <- c(0,24,72)

#produce list of production rates by time points
prod.rate_list <- lapply(1:2, function(x) {
                      #define time frame
                      t.start <- timePoints[x]
                      t.end <- timePoints[x+1]
                      
                      #calculate change in MVs abundance
                      MVs_delta_calc<- MVs_CC_bio_rep %>% 
                        dplyr::filter(Time %in% c(t.start,t.end), Type == "MVs") %>% 
                        group_by(Strain, Bio_rep) %>%
                        mutate(DeltaParticles = Concentration.mean[Time == t.end] - Concentration.mean[Time == t.start])
                      
                      #calculate number of generations
                      generation_calc<- MVs_CC_bio_rep %>% 
                        dplyr::filter(Time %in% c(t.start,t.end), Type == "Cells") %>% 
                        group_by(Strain, Bio_rep) %>%
                        mutate(growth_rate = log(Concentration.mean[Time == t.end]/ Concentration.mean[Time == t.start])/(t.end-t.start),
                               generation = (t.end-t.start)*growth_rate/log(2))
                      
                      #calculate production rates
                      prod.rate.calc <- left_join(generation_calc[,c("Strain", "Bio_rep","Time","generation","Concentration.mean")], 
                                                  MVs_delta_calc[,c("Strain", "Bio_rep","Time","DeltaParticles")], by=c("Strain","Time"))  %>% 
                        dplyr::filter(Time == t.start) %>% # to get the initial cell counts
                        mutate(V = DeltaParticles, 
                               x = (2^(generation)-1)* Concentration.mean,
                               prod.rate = V / x,
                               Period = paste0(t.start,"-",t.end))
                      
                      return(prod.rate.calc)
                    })

#aggregate into a dataframe
prod.rate_df <- bind_rows(prod.rate_list)

#plot prod.rate
prod.rate_plot<- prod.rate_df %>% 
  mutate(Strain = factor(Strain, levels = c("AD45","ATCC27126","BGP6","BS11","HOT1A3","MIT1002"))) %>% 
  ggplot(aes(x=Strain, y = round(prod.rate,0), group = Period, label = round(prod.rate,0)))+
  geom_col(width = 1, color = "black")+
  geom_text(nudge_y = 1, size = 5 )+
  geom_errorbar(aes(ymin = prod.rate-prod.rate.se, ymax = prod.rate+prod.rate.se),  width = 0.2)+
  labs(y="MVs production rate (# per cell per generation)", x = "Strain")+
  facet_wrap(Period~., ncol = 1, scales="free_y")+
  #scale_y_log10()+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

##########################################
### Produce merged figure
##########################################

plot_grid(concentrations_plot, prod.rate_plot)
