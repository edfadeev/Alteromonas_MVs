require(tidyr)
require(dplyr)
require(purrr)
require(phyloseq)

#define working directory with the files
wd <- "Data/PD_results/"

#set variables
strains <- c("Amac_ATCC27126","Amac_AD45","Amac_HOT1A3","Amac_BS11",
             "Amac_MIT1002","Amac_BGP6")

source_dbs <- c("NCBI_PGAP","COG20_FUNCTION","COG20_CATEGORY","COG20_PATHWAY",
                "GO","Pfam","InterPro","MetaCyc", 
                "KOfam","KEGG_Class","KEGG_Module") 


################################################################################
#generate phyloseq objects for each strain and export it
################################################################################
for (str in strains) {
  #############################################
  #import annotation of the genomes from anvio
  #############################################
  #gene_type <- unique_genes %>%  filter(Strain == str)
  
  #genes list
  gene_annotations_df <- read.csv(paste0(wd,"Ref_annotations/",str,"-genes.txt"),
                                  sep="\t", h= T) #%>% 
                          #left_join(gene_type, by = "gene_callers_id")
    
  
  
  source_files <- list.files(path = paste0(wd,"Ref_annotations/"), 
                              pattern = paste0("*",str,"-.*functions\\.txt"),
                              full.names = TRUE)
  
  annotations_list <- lapply(source_files, function(x) {
                              source_db <- unlist(strsplit(basename(x), split ="-"))[2]
                              dat <- read.csv(x,sep="\t", 
                                      col.names = c("gene_callers_id","db",
                                                    "accession","function.","e_value")) %>% 
                                      #mutate_all(as.character) %>% 
                                select(gene_callers_id, accession, function.) %>% 
                                group_by(gene_callers_id)%>%
                                summarise_each(funs(paste(unique(.), collapse='|')), matches('^\\D+$')) %>% 
                                plyr::rename(replace= c("accession"=paste(source_db, "accession", sep ="_"), 
                                                        "function."=paste(source_db, "function", sep ="_"))) 
                              return(dat)
                              } )

  #merge all dataframes
  annotations_df <- annotations_list %>% 
                    purrr::reduce(full_join, by = "gene_callers_id") %>% 
                    select(gene_callers_id,contains(source_dbs)) %>% 
                    replace(is.na(.), "Unk")

  gene_meta_df <- gene_annotations_df %>% 
                  select(gene_callers_id, contig, start, stop, partial, aa_sequence) %>% 
                  merge(annotations_df, by = "gene_callers_id") %>% 
                  mutate(prot_length = nchar(aa_sequence)) %>% 
                  mutate(rows = gene_callers_id) %>% 
                  tibble::column_to_rownames('rows')
    
  #gene ids as taxonomy table
  genome_annotation<- tax_table(as.matrix(gene_meta_df))
  #############################################
  #import protein profiles from Protein Discoverer
  #############################################
  sample_names <- c(MP_1 = "S_1", MP_2 = "S_2", EV_1 = "S_3", EV_2 = "S_4")
  
  metaP_raw <- read.csv(paste0(wd,str,"-Proteins.txt"),
                        sep="\t", h= T)%>% 
    dplyr::rename(gene_caller_id = Accession) %>% 
    mutate_if(is.numeric, funs(replace_na(., 0))) %>% 
    filter(Number.of.PSMs >=2 , Number.of.Unique.Peptides>=1)%>% 
    rename_with(~gsub("Abundance.F","S_",.), everything()) %>% 
    rename_with(~gsub(".Sample","",.), everything()) %>% 
    dplyr::rename(any_of(sample_names)) %>%
    select(any_of(c("gene_caller_id", "MP_1", "MP_2", "EV_1", "EV_2")))
    
  
  totalP <- metaP_raw %>% 
    replace(is.na(.), 0) %>%
    tibble::column_to_rownames('gene_caller_id') %>%
    otu_table(taxa_are_rows = TRUE)
  
            

  #############################################
  #generate phyloseq object
  #############################################
  sam_data <- data.frame(Strain = str,
              Fraction = gsub("_.*","",colnames(totalP)),
              Replicate = gsub(".*_","",colnames(totalP)),
              row.names = colnames(totalP)) %>% 
              sample_data()

  
  phyloseq <- phyloseq(totalP, genome_annotation, 
                       sam_data)
  
  #remove empty observation
  phyloseq_clean <- prune_taxa(taxa_sums(phyloseq)>0, phyloseq)
  
  #export phyloseq object
  saveRDS(phyloseq_clean, paste0("Data/PD2phy/", str,"_phyloseq.rds"))
}

#clean workspace
rm(list = ls())
gc()
