# I am confident about the majority of the script. However, the DESeq2 analysis at the end will need looking at.
# I am not entirely clear whether the multifactorial analysis is applied correctly yet. I will continue working on
# it but figured you might want to work on your own stuff in the meantime. You may wish to restrict the analysis to
# the samples of interest only in order to simplify things.

# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Tips from:
# https://scilifelab.github.io/courses/rnaseq/labs/kallisto
# https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
setwd("/Volumes/Maxtor/Rita/trna_12082018/tRNA_2017")
library(tximport)
library(readr)
library(DESeq2)
library(biomaRt)

# Prepare 2 tables
# 1) Describing each sample
# 2) Transcript to gene conversion

# 1) Description table:
alignment_dir <- "dicer_dep_expression.own_data/"
samples <- list.files(alignment_dir,"*070817*")
samples
alignment_directories <- paste(alignment_dir,samples,sep="")
samples <- sub("^(\\w+)\\..*","\\1",samples,perl=TRUE)
samples
desc_table <- data.frame(path=alignment_directories,sample=samples,stringsAsFactors=FALSE)

sample_pairings <- c(rep("HEK_induced",3),rep("HEK_control",3),rep("Sh_induced",3),rep("Sh_control",3))
names(sample_pairings) <- c("CoA","CoB","CoC","CoX","CoY","CoZ","ShA","ShB","ShC","ShX","ShY","ShZ")
desc_table$cell_type <- sample_pairings[desc_table$sample]

# CHECKED

# 2) Transcript to gene conversion
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
assign_transcripts_to_genes <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = ensembl)



# CHECKED

# Read in the Kallisto data
# NOTE: Error in hd5 import
# Error in make.names(col.names, unique = TRUE) : 
# invalid multibyte string at '<89>HDF'

files <- paste(desc_table$path,"abundance.tsv",sep="/")
names(files) <- desc_table$sample
if(!all(file.exists(files))){stop("Missing an abundance file")}
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kallisto

# CHECKED

# Convert to gene counts
txi.kallisto.gene <- summarizeToGene(txi.kallisto,assign_transcripts_to_genes[,c("ensembl_transcript_id","ensembl_gene_id")],ignoreTxVersion=TRUE)
rownames(txi.kallisto.gene$counts)
rowSums(txi.kallisto.gene$counts)
write.table(rowSums(txi.kallisto.gene$counts),file="allgenes_kallistocounts",quote=F,col.names=F)
#Prepare for DESeq2
rownames(desc_table) <- desc_table$sample
if(! all(colnames(txi.kallisto.gene)%in%rownames(desc_table))){stop("Missing some samples in the Kallisto sample file")}

# CHECKED


### FROM HERE PROCEED WITH EXTREME CAUTION:
### NEED TO LEARN MORE CONCERNING INTERACTION TERMS

desc_table$cell_type <- as.factor(desc_table$cell_type)
desc_table$cell <- as.factor(c(rep("HEK293T",6),rep("2B",6)))
desc_table$dox <- as.factor(rep(c(rep("yes",3),rep("no",3)),2))
desc_table<-desc_table[desc_table$dox=="yes",]# See ?results for a description of what is going on here
names(txi.kallisto.gene)

txi.kallisto.gene$length<-txi.kallisto.gene$length[,c("CoA","CoB","CoC","ShA","ShB","ShC")]
txi.kallisto.gene$counts<-txi.kallisto.gene$counts[,c("CoA","CoB","CoC","ShA","ShB","ShC")]
txi.kallisto.gene$abundance<-txi.kallisto.gene$abundance[,c("CoA","CoB","CoC","ShA","ShB","ShC")]
txi.kallisto.gene
desc_table
dds <- DESeqDataSetFromTximport(txi.kallisto.gene, desc_table, ~cell_type) 

dds <- DESeq(dds)
resultsNames(dds)



sh_co_results <- results(dds, contrast=c("cell_type","Sh_induced","HEK_induced"))


sh_co_results
hist(sh_co_results$pval)
(sh_co_results <- sh_co_results[order(sh_co_results$padj),])
# trna_prev_vector<-c("ENST00000262776","ENST00000371655","ENST00000380486","ENST00000361566","ENST00000244519","ENST00000372040","ENST00000259526","ENST00000378943","ENST00000589390","ENST00000316403","ENST00000339364","ENST00000338702","ENST00000267842","ENST00000263816","ENST00000526355","ENST00000396625","ENST00000380861","ENST00000581347","ENST00000370819","ENST00000452846","ENST00000525985","ENST00000227507","ENST00000377113","ENST00000315901","ENST00000591639","ENST00000344051","ENST00000272223")
# trna_prev_trans_genes<-assign_transcripts_to_genes[assign_transcripts_to_genes$ensembl_transcript_id%in%vector,]
# trna_prev_trans_genes$ensembl_gene_id
# View(sh_co_results[trna_prev_trans_genes$ensembl_gene_id,])
(sh_co_results_rest<-subset(sh_co_results, padj < 10^(-3)))
seq = getSequence(id=row.names(sh_co_results_rest), type="ensembl_gene_id", seqType="gene_exon_intron", mart = ensembl)
exportFASTA(seq,file="sh_co_results_rest_sequence_test.fasta")

row.names(sh_co_results_rest)
g = getBM(c("hgnc_symbol","chromosome_name","strand","start_position","end_position","ensembl_gene_id"),"ensembl_gene_id", row.names(sh_co_results_rest) , ensembl)
show(g)
write.table(g, "~/Desktop/sh_co_results_mat_coordinates.txt",  quote=FALSE, sep="\t")

write.table(sh_co_results, "~/Desktop/sh_co_results_mat.txt",  quote=FALSE, sep="\t")
