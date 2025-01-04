#Before going into R a match list needs to be created this is simply done by the following bash script
#you´ll need blast+ for this and an ASV.fa file of your ASV.seqs
makeblastdb -in ASVs.fa -parse_seqids -dbtype nucl
blastn -db ASVs.fa -outfmt '6 qseqid sseqid pident' -out match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query ASVs.fa

#This is the R part of the lulu process
library(lulu)
## Import data
otutab <- read.csv("ASV_counts_deconprev05.tsv",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table("match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
tax <- read.csv("ASVs_taxonomy_deconprev05.tsv",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)



curated_result <- lulu(otutab, matchlist)



### Compared output
curated_result$curated_count

percentage.removed <- round((as.numeric(length(rownames(otutab)))-curated_result$curated_count)/as.numeric(length(rownames(otutab)))*100,digits = 2)



Curated <- curated_result$curated_table

tax_curated <- tax[match(rownames(Curated),rownames(tax)),]


write.table(tax_curated, "Curated_Tax.csv", sep = ',', row.names = FALSE)                            
write.table(curated_result$curated_table, "Curated_Table.txt", sep = ',')  
