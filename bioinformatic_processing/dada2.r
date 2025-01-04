#set dependencies
library(dada2)
packageVersion("dada2") 
list.files() # make sure what we think is here is actually here


#r set create data and check QC
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
path <- "trimmed"
# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern="_R1_001_val_1.fq.gz", full.names = TRUE))

# and one with the reverse
fnRs <- sort(list.files(path, pattern="_R2_001_val_2.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1_001_val_1"), `[`, 1)

#Filtering
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filtered_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,trimLeft=10, truncLen=c(290,200), compress=TRUE, multithread=TRUE)
### Store filtered reads as R dataset
saveRDS(filtered_out, "filtered_out.RData")
#r filtering based on error rates
#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtFs, multithread=TRUE)
err_reverse_reads <- learnErrors(filtRs, multithread=TRUE)

#dereplication
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

#Inferring of ASVs

#Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)

#Merge reads
#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, verbose=TRUE)
### Store filtered amplicons as R dataset
saveRDS(merged_amplicons, "merged_amplicons.RData")

###Generate ASV table
#Generating a count table
seqtab <- makeSequenceTable(merged_amplicons)
table(nchar(getSequences(seqtab)))
#to see the breakdown of the size of these amplicons. On the top row is the size of the merged reads, and on the botton is the frequency

#Chimera identification
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.RData")

### Create a Summary
# this is one quick way to look at sequences that have been lost, to know whether they held a lot in terms of abundance
sum(seqtab.nochim)/sum(seqtab)

#Overview of counts throughout:quick way to pull out how many reads were dropped at various points of the pipeline
# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
write.table(summary_tab, "summary_reads_table_CRISPR_full_new_Trunclen200_180_O2.txt", sep="\t", quote=F)


#Assigning taxonomy(before this step, download the database, at: https://zenodo.org/record/1172783#.XzvzQJMzbmE)
#                  There are different DADA2-formatted databases available in DADA2 website

#Assign Taxanomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", tryRC=T)
saveRDS(taxa, "taxa.RData")


#Extracting the standard goods from DADA2
#Write out tables for further processing}
#giving to seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
  }

#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
