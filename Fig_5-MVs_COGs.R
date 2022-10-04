###########################################
#Load libraries
###########################################
require(ggplot2)
require(phyloseq)
require(dplyr)
require(tidyr)
require(purrr)

source("./scripts/nsaf_trans.R") #NSAF trasformation function

########################################
#import phyloseq
########################################
#import the output of Proteome Discoverer into phyloseq objects (run only once)
#source("./scripts/PD2phyloseq.R")

#import the phyloseq objects 
strains <- c("Amac_ATCC27126","Amac_AD45","Amac_HOT1A3","Amac_BS11",
             "Amac_MIT1002","Amac_BGP6")
for (str in strains) {
  assign(paste0(str,"_phy0"),readRDS(paste0("Data/PD2phy/", str,"_phyloseq.rds")))
}

#extract all annotations for each strain
genes_description <- lapply(strains,  function(x) {
  desc<- as.data.frame(tax_table(get(paste0(x,"_phy0")))) %>% 
    select("gene_callers_id", contains("_function"),contains("_accession"),"aa_sequence") %>% 
    mutate(Strain = x) %>% 
    #tidyr::separate_rows(c(COG20_FUNCTION_accession, COG20_FUNCTION_function), sep = '!!!') %>% #break double annotations into individual rows
    distinct()
  
  return(desc)
})

genes_description_df <- bind_rows(genes_description)

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
             paste0(wd,"Localization/CELLO2GO/",x,"_EV_top_proteins.fasta"))
  
  
  top_prot<- rbind(top_prot_MP,top_prot_EV)
  
  return(top_prot)
})

names(top_prot_list) <- strains 

#overview of number of proteins
top_prot_df <- bind_rows(top_prot_list) 

########################################
#explore EV proteins on COG level
########################################
#extract COGs descriptions
COGs_description <-  genes_description_df %>% select(contains("COG"))%>% 
  distinct()

top_prot_COG<- top_prot_df %>% 
  replace(is.na(.), "Unk") %>% 
  group_by(Strain,Sample, COG20_FUNCTION_accession) %>%
  summarize(N_genes= n()) %>% 
  mutate(Proportion = round(100*(N_genes / sum(N_genes)),1)) %>% 
  filter(COG20_FUNCTION_accession != "Unk") %>% 
  left_join(COGs_description, by = "COG20_FUNCTION_accession")

prot_COG_sum<- top_prot_COG %>% 
  group_by(Sample, COG20_FUNCTION_accession,COG20_FUNCTION_function) %>% 
  summarize(N_strains = n())

#extract COGs list of EV proteins for each strain
EVs_COGs_list <- sapply(strains, simplify = FALSE, USE.NAMES =FALSE, function(x) {
  COGs<- top_prot_list[[x]] %>%
    dplyr::select(COG20_FUNCTION_accession, Abundance) %>%
    replace(is.na(.), "Unk") %>% #replace NAs - needed for pres_abs function
    group_by(COG20_FUNCTION_accession) %>%
    summarize(Total.abundance =sum(Abundance))
  return(COGs)
})

#add strains as labels for each vector in the list
names(EVs_COGs_list) <- strains 


#merge COGs counts with description
EVs_COGs_df <- EVs_COGs_list %>% 
  purrr::reduce(full_join, by = "COG20_FUNCTION_accession") %>% 
  setNames(c("COG20_FUNCTION_accession",strains)) %>% 
  replace(is.na(.), 0) %>% 
  left_join(COGs_description, by = "COG20_FUNCTION_accession") %>% 
  separate_rows(c(COG20_CATEGORY_accession, COG20_CATEGORY_function), sep = '!!!')#break double annotations into individual rows


#summarize the number of proteins in each COG category
EVs_COGs_sums<- EVs_COGs_df %>% 
  melt() %>% 
  group_by(variable, COG20_CATEGORY_accession, COG20_CATEGORY_function) %>% 
  summarize(Total_proteins = sum(value)) %>% 
  filter(Total_proteins>0) %>% 
  group_by(variable) %>% 
  mutate(Proportion = 100*(Total_proteins / sum(Total_proteins))) 


#select the most prominent COG categories
selected_pathways <- c("Energy_production_and_conversion",
                       "Amino_acid_transport_and_metabolism",
                       "Nucleotide_transport_and_metabolism",
                       "Carbohydrate_transport_and_metabolism",
                       "Cell_wall/membrane/envelope_biogenesis",
                       "Coenzyme_transport_and_metabolism", 
                       "Signal_transduction_mechanisms",
                       "Translation,_ribosomal_structure_and_biogenesis",
                       "Intracellular_trafficking,_secretion,_and_vesicular_transport",
                       "Inorganic_ion_transport_and_metabolism",
                       "Posttranslational_modification,_protein_turnover,_chaperones")

#plot the proportion of proteins out of total in each category
EVs_COGs_sums %>% 
  #filter(COG20_CATEGORY_function %in% selected_pathways) %>% 
  ggplot(aes(x=COG20_CATEGORY_function, y= Proportion, fill = variable))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = tol21rainbow)+
  #facet_wrap(variable~.)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

ggsave(paste0(wd,"R_data/figures/","EVs_COGs_sums.pdf"),
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

