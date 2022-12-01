###########################################
#Load libraries
###########################################
require(ggplot2)
require(phyloseq)
require(dplyr)
require(tidyr)

source("scripts/write_fasta.R")
source("scripts/nsaf_trans.R")


cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", 
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


########################################
#import phyloseq
########################################
#import the phyloseq objects 
strains <- c("Amac_ATCC27126","Amac_AD45","Amac_HOT1A3","Amac_BS11",
             "Amac_MIT1002","Amac_BGP6")
for (str in strains) {
  assign(paste0(str,"_phy0"),readRDS(paste0("Data/PD2phy/", str,"_phyloseq.rds")))
}

########################################
#export top 5% of proteins
########################################
top_prot_list<- lapply(strains, function(x) {
  prot_nsaf <- get(paste0(x,"_phy0")) %>% 
    merge_samples("Fraction", fun = "sum") %>% 
    add_nsaf(.,"prot_length") %>% 
    psmelt()%>% 
    filter(Abundance > 0) %>%
    mutate(Strain = x) %>% 
    dplyr::rename(gene_caller_id = OTU) %>% 
    select(-c("contig","Fraction","Replicate", "gene_callers_id","start","stop","partial"))
  
  top_prot_MP<- prot_nsaf %>% 
    filter(Sample == "MP") %>%
    slice_max(Abundance, prop = 0.05) 
  
  writeFasta(top_prot_MP[,c("gene_caller_id","aa_sequence")],
              paste0("Data/Localization/CELLO2GO/",x,"_MP_top_proteins.fasta"))
  
  
  top_prot_EV<- prot_nsaf %>% 
    filter(Sample == "EV") %>%
    slice_max(Abundance, prop = 0.05) 
  
  writeFasta(top_prot_EV[,c(c("gene_caller_id","aa_sequence"))],
          paste0("Data/Localization/CELLO2GO/",x,"_EV_top_proteins.fasta"))
  
  
  top_prot<- rbind(top_prot_MP,top_prot_EV)
  
  return(top_prot)
})

#overview of number of proteins
top_prot_df <- bind_rows(top_prot_list) 

#upload the fasta manually to http://cello.life.nctu.edu.tw/cello2go/

########################################
# Analyze localization in Cello2go
# http://cello.life.nctu.edu.tw/cello2go/
########################################
#export the top 5% of all proteins in each fraction as fasta files 
#(no need to run)
#source("./scripts/top_prot2fasta.R")

#plot results
cello2go_list <- lapply(strains, function(x) {
  dat_EV <- read.csv(paste0("Data/Localization/CELLO2GO/",x, "_EV_out.txt"),
                     h=F, sep =" ") %>% 
    filter(V1 =="LA") %>% 
    separate(V4, into =c("Localization","Amount"), sep =":") %>% 
    mutate(Amount = gsub(";","",Amount),
           Strain = x,
           Fraction = "EV") %>% 
    select(c("Strain","Fraction","Localization","Amount")) 
  
  
  dat_MP <- read.csv(paste0("Data/Localization/CELLO2GO/",x, "_MP_out.txt"),
                     h=F, sep =" ") %>% 
    filter(V1 =="LA") %>% 
    separate(V4, into =c("Localization","Amount"), sep =":") %>% 
    mutate(Amount = gsub(";","",Amount),
           Strain = x,
           Fraction = "MP") %>% 
    select(c("Strain","Fraction","Localization","Amount"))
  
  dat <- rbind(dat_MP,dat_EV) %>% 
    group_by(Fraction) %>% 
    mutate(Amount = as.numeric(Amount),
           Proportion = 100*(Amount / sum(Amount)))
  
  return(dat)
} )

#overview of number of proteins
cello2go_predictions_df <- bind_rows(cello2go_list) %>% 
                            mutate(Source ="cello2go")

########################################
# localization predictions from psortdb
# https://db.psort.org/
########################################
# import predictions (available only for 4 strains)
psortdb_list <- lapply(c("Amac_ATCC27126","Amac_AD45",
                         "Amac_HOT1A3","Amac_BS11"), 
                       function(str) {
                            psortdb_raw <- read.table(paste0("Data/Localization/psortdb/",str,"-psortdb.tsv"),
                                                      sep ="\t", h= T, fill = TRUE) %>% 
                                           mutate(Strain = str,
                                                  NCBI_PGAP_accession = gsub("ref\\|","",SeqID)) %>% 
                                           select(c("Strain","NCBI_PGAP_accession", "Final_Localization",         
                                                    "Final_Localization_2", "Final_Localization_Comments",
                                                    "Secondary_Localization", "Final_Score"))
                              return(psortdb_raw)
                                } )

#localization per protein
top_prot_psortdb<- top_prot_df %>% 
  filter(Strain %in% c("Amac_ATCC27126","Amac_AD45","Amac_HOT1A3","Amac_BS11")) %>% 
  left_join(bind_rows(psortdb_list), by = c("Strain","NCBI_PGAP_accession"))

#summarize
top_prot_psortdb_sum<- top_prot_psortdb %>% 
  group_by(Strain, Sample, Final_Localization) %>% 
  dplyr::rename(c("Localization" ="Final_Localization",
                  "Fraction" = "Sample")) %>% 
  mutate(Localization = case_when(Localization== "Cytoplasmic Membrane"~ "Innermembrane",
                                  Localization== "Outer Membrane" ~ "Outermembrane",
                                  is.na(Localization)==TRUE ~"Unknown",
                                  TRUE ~ Localization)) %>% 
  summarize(Amount = n()) %>% 
  mutate(Proportion = round(100*(Amount / sum(Amount)),1),
         Source ="psortdb")

########################################
# Plot both localization predictions
########################################
prot_loc_df <- rbind(top_prot_psortdb_sum, 
                     cello2go_predictions_df) %>% 
                mutate(Localization = as.factor(Localization))

ggplot(prot_loc_df, 
  aes(x= interaction(Fraction,Strain), y= Proportion, group = Strain, fill = Localization))+
  geom_col()+
  scale_fill_manual(values = rev(cbbPalette))+
  facet_wrap(~Source, scales = "free_x")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x =element_text(angle = 90),
        legend.position = "bottom")

ggsave(paste0("Figures/","Fig_4-localizations.pdf"),
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)
