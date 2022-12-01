###########################################
#Load libraries
###########################################
require(DESeq2)
require(ggplot2)
require(phyloseq)
require(dplyr)
require(tidyr)
require(vegan)

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

########################################
#explore the protein dataset
########################################
#generate dataframe with relative abundance of each protein
prot_nsaf_list <- lapply(strains, function(x) {
  prot_nsaf <- get(paste0(x,"_phy0")) %>% 
    merge_samples("Fraction", fun = "sum") %>% 
    add_nsaf(.,"prot_length") %>% 
    psmelt()%>% 
    filter(Abundance > 0) %>%
    mutate(Strain = x) %>% 
    dplyr::rename(gene_caller_id = OTU) %>% 
    select(-c("contig","Fraction","Replicate", "gene_callers_id","start","stop","partial"))
  return(prot_nsaf)
})  

names(prot_nsaf_list) <- strains 
prot_nsaf_df <- bind_rows(prot_nsaf_list) 

#overview of number of proteins
n_proteins_by_strain<- prot_nsaf_df %>% 
  group_by(Strain, Sample) %>% 
  summarize(Total = n(),
            Mean.abund = mean(Abundance),
            SD.abund = sd(Abundance))

#Number of proteins per genome
n_prot_genome <- data.frame(Strain = strains,
                            Proteins = c(3770,3817,3839,3677,3845,3869))

#calculate proportion of detection
prot_proportion <- n_proteins_by_strain %>% 
  left_join(n_prot_genome) %>% 
  mutate(Proportion = 100*Total/Proteins)


prot_proportion_mean <- prot_proportion %>% 
  mutate(Fraction = gsub("_1|_2","",Sample)) %>% 
  group_by(Fraction) %>% 
  summarize(Mean.cov =mean(Proportion),
            SD.cov = sd(Proportion))

########################################
#Test whether the fractions are significantly different in each strain
########################################
prot_sig_list <- lapply(strains, function(x) {
  prot_phy0 <- get(paste0(x,"_phy0")) 
  
  counts <- as(otu_table(prot_phy0),"matrix")
  
  counts <- as.data.frame(counts) %>% 
    mutate_if(is.numeric, as.integer) %>% 
    replace(is.na(.), 0)
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sample_data(prot_phy0),
                                design = ~ Fraction)
  
  
  #conduct variance stabilization of the dataset
  metaP.DEseq.vsd <- vst(dds,  fitType="local", blind=FALSE)
  
  #test differences
  df <- colData(metaP.DEseq.vsd)[,c("Replicate","Fraction")]
  sample_distance <- dist(t(assay(metaP.DEseq.vsd)))
  adonis_all <- adonis2(sample_distance ~ Replicate+Fraction, df)
  adonis_all
  
  #posthoc to check which ponds are different
  groups <- df[["Fraction"]]
  mod <- betadisper(sample_distance, groups)
  perm <- permutest(mod)
  
  return(perm)
})  

names(prot_sig_list) <- strains 

########################################
#Explore how many MV proteins shared between strains
########################################
#import the list of gene clusters
gene_cluster <- read.csv("Data/anvio/Alteromonas-Pangenome_gene_clusters_summary.txt", sep="\t", h= T) %>% 
  select("genome_name", "gene_callers_id", "gene_cluster_id")


MV_prot_overlap_list <- lapply(strains, function(x) {
  gene_cluster_strain <- gene_cluster %>%  filter(genome_name == x)
  
  prot_EV <- prot_nsaf_list[[x]] %>%
    select(c("gene_caller_id","Sample","Strain","Abundance")) %>% 
    dplyr::rename(gene_callers_id = gene_caller_id) %>% 
    mutate(gene_callers_id = as.integer(gene_callers_id)) %>% 
    filter(Sample =="EV", Abundance >0.00) %>% 
    left_join(gene_cluster_strain, by = "gene_callers_id") %>% 
    select(gene_cluster_id,Strain)
  
    return(prot_EV)
      })

MV_prot_overlap_df <- MV_prot_overlap_list %>%  reduce(full_join, by = c("gene_cluster_id")) %>% 
                      rename_with(~ c("gene_cluster_id", strains)) %>% unique() %>% 
                      mutate_at(vars(strains), funs(case_when(is.na(.)~0, TRUE ~ 1))) %>% 
                      group_by(gene_cluster_id) %>% 
                      mutate(Total = sum(across(strains)))

#total protein clusters shared between all strains
MV_prot_overlap_df %>% filter(Total == 6) %>% dim()



########################################
#Explore the most abundant proteins
########################################
#calculate overlaps of proteins between fractions
prot_overlap_list <- lapply(strains, function(x) {
  prot_MP <- prot_nsaf_list[[x]] %>%
    select(c("gene_caller_id","Sample","Abundance")) %>% 
    filter(Sample =="MP", Abundance >0.00) %>% 
    #filter(Sample == "MP") %>%
    slice_max(Abundance, prop = 0.05) 
  
  prot_EV <- prot_nsaf_list[[x]] %>%
    select(c("gene_caller_id","Sample","Abundance")) %>% 
    filter(Sample =="EV", Abundance >0.00) %>% 
    #filter(Sample == "EV") %>%
    slice_max(Abundance, prop = 0.05) 
  
  overview <- data.frame(Strain = x,
                         shared = length(intersect(prot_MP$gene_caller_id,prot_EV$gene_caller_id)),
                         unique_MP = length(setdiff(prot_MP$gene_caller_id,prot_EV$gene_caller_id)),
                         unique_EV = length(setdiff(prot_EV$gene_caller_id,prot_MP$gene_caller_id))) %>% 
    mutate(Prop_u_MP = unique_MP/(unique_MP+shared),
           Prop_u_EV= unique_EV/(unique_EV+shared))
  return(overview)
})

#transform into dataframe
prot_overlap_df <- bind_rows(prot_overlap_list) %>% 
  left_join(n_prot_genome)

prot_overlap_df %>% summarize(Prop_u_MP.mean = mean(Prop_u_MP), Prop_u_MP.sd = sd(Prop_u_MP), 
                              Prop_u_EV.mean = mean(Prop_u_EV), Prop_u_EV.sd = sd(Prop_u_EV))


