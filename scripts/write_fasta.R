#write out fasta file from a dataframe 
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"gene_caller_id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"aa_sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}