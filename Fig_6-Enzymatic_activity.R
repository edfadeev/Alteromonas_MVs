require(rio)
require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpmisc)
require(purrr)
require(reshape2)

###########################################
#Create standard curves for MUF and MCA
###########################################
#import data
curves_raw <- read.csv("data/EEA_cal_curves.csv",
                       header=T, row.names=1, sep=",")

#calculate mean per each concentration and subtract the blank
std_curves <- curves_raw %>%
  group_by(t, sample, fluorphore, concentration_nM) %>%
  summarize(fl_mean = mean(fluorescence)) %>% #calculate mean
  group_by(t) %>% 
  mutate(fl_corr = fl_mean - fl_mean[sample == "asw"]) %>% #subtract blank
  filter(sample != "asw") %>% 
  ungroup()

#plot regressions
std_curves %>% 
  ggplot(aes(concentration_nM, fl_corr, colour = sample)) + 
  geom_point(size =3)+
  geom_smooth(method = "lm")+
  scale_y_continuous(labels=scales::scientific_format())+
  stat_poly_eq(formula = y~x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,
               colour = "black") +
  facet_wrap(fluorphore~t, scales = "free_y")+
  theme_bw(14)

#produce linear models for each timepoint
std_curves.lm <-std_curves %>% 
  select(t, fluorphore,concentration_nM, fl_corr) %>% 
  group_by(fluorphore, t) %>% 
  do(fit.lm=lm(concentration_nM ~ fl_corr, .))

###########################################
#Calculate concentrations based on calibration curves
###########################################
#import data
raw_measurements <- read.csv("data/EEA_measurements.csv",
                             header=T, row.names=1, sep=",")

#calculate the mean fluorescence and substruct the blank
data_corr <- raw_measurements %>%
  group_by(t, time_h, sample, substrate) %>%
  summarize(fluorescence_mean = mean(fluorescence)) %>% 
  group_by(t, time_h, substrate) %>%
  mutate(fl_corr = fluorescence_mean - fluorescence_mean[sample == "NC"] ) %>%   #substracts the blank
  filter(sample != "NC") %>% 
  separate(col="substrate", into=c("fluorphore","substrate"),
           sep="-", extra ="merge", remove=FALSE) %>% 
  ungroup()


#calculate substrate conc. based on the linear models of the calibration curves
sample_data.conc <- data_corr %>% 
  group_by(fluorphore, time_h, t) %>% #media
  nest() %>% 
  inner_join(std_curves.lm) %>% 
  mutate(subs_conc_nM = map2(fit.lm, data, predict)) %>% 
  unnest(c(data,subs_conc_nM)) %>% 
  separate(col="sample", into=c("fraction","replicate"), sep="_", remove=FALSE)

#############################
# Calculate enzymatic activity
#############################
sample_data.act <- sample_data.conc %>% 
  filter(t %in% c("t0","t1")) %>% 
  group_by(fraction, replicate, substrate) %>% 
  arrange(time_h) %>%
  mutate(activity.nM_h = (subs_conc_nM - lag(subs_conc_nM))/ (time_h-lag(time_h))) %>% 
  mutate(activity.nM_h = case_when(fraction =="100kDa-concentrate" ~ 2*activity.nM_h, # the MVs were diluted 2X in ASW
                                   #activity.nM_h < 0 ~ 0,
                                   TRUE ~ activity.nM_h)) %>% 
  filter(t =="t1", 
         activity.nM_h> 0) %>% 
  group_by(fraction, substrate) %>% 
  summarize(act.mean = mean(activity.nM_h), act.se = se(activity.nM_h))



#plot enzyme activity vs. substrate concentration
sample_data.act %>% 
  filter(fraction != "total-fraction") %>% 
  ggplot(aes(x = fraction, act.mean, fill = fraction))+ #fill=fraction
  geom_col()+
  geom_errorbar(aes(ymin=act.mean-act.se, ymax=act.mean+act.se),
                width=.2) +
  labs(y = bquote('EEA (nmol '~L^-1~~h^-1~")")) +
  facet_wrap(substrate~., scales = "free_y", ncol = 2)+
  theme_bw(base_size=14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")





sample_data.act.prop <- sample_data.act %>% 
  group_by(replicate, substrate,t) %>% 
  mutate(activity.prop.022 = activity.nM_h[fraction=="0.2um-filtrate"]/activity.nM_h[fraction=="total-fraction"],
         activity.prop.MVs = activity.nM_h[fraction=="100kDa-concentrate"]/activity.nM_h[fraction=="total-fraction"],
         activity.prop.free = activity.nM_h[fraction=="100kDa-outflow"]/activity.nM_h[fraction=="total-fraction"]) %>% 
  select(t, replicate, substrate, activity.prop.022, activity.prop.MVs, activity.prop.free) %>% 
  melt() %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(value > 0) %>%
  group_by(variable, substrate,t) %>% 
  summarize(prop.mean = mean(value), prop.se = se(value))





sample_data.act.prop %>% 
  ggplot(aes(x=variable, y = prop.mean)) +
  geom_col()+
  geom_errorbar(aes(ymin=prop.mean-prop.se, ymax=prop.mean+prop.se),
                width=.2) +
  facet_wrap(substrate~t, scales = "free_y", ncol = 2)+
  theme_bw(base_size=14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")






# New facet label names for substrate variable
enzyme.labs <- c("LAPase", "AGase", "BGase", "APase")
names(enzyme.labs) <- c("MCA-Leucine", "MUF-alpha-glucoside", "MUF-beta-glucoside", "MUF-Phosphate")



ggsave("EEA_fraction.pdf",
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)
