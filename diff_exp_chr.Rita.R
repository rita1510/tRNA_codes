# I am confident about the majority of the script. However, the DESeq2 analysis at the end will need looking at.
# I am not entirely clear whether the multifactorial analysis is applied correctly yet. I will continue working on
# it but figured you might want to work on your own stuff in the meantime. You may wish to restrict the analysis to
# the samples of interest only in order to simplify things.

# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# Tips from:
# https://scilifelab.github.io/courses/rnaseq/labs/kallisto
# https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
setwd("~/Desktop/")
library(tximport)
library(readr)
library(DESeq2)
library(biomaRt)

# Prepare 2 tables
# 1) Describing each sample
# 2) Transcript to gene conversion

# 1) Description table:
alignment_dir <- "kallisto_chr/"
samples <- list.files(alignment_dir,"WTCHG_6*")
samples
alignment_directories <- paste(alignment_dir,samples,sep="")
samples <- sub("^(\\w+)\\..*","\\1",samples,perl=TRUE)
samples
desc_table <- data.frame(path=alignment_directories,sample=samples,stringsAsFactors=FALSE)

sample_pairings <- c(rep("HEK293_chr",2),rep("shDrosha",2),rep("shDicer",2),rep("shAgo2",2))
names(sample_pairings) <- c("WTCHG_6_205142","WTCHG_6_206130","WTCHG_6_207118","WTCHG_6_208106","WTCHG_6_209189","WTCHG_6_210177","WTCHG_6_211165","WTCHG_6_212153")
desc_table$cell_type <- sample_pairings[desc_table$sample]

# CHECKED

# 2) Transcript to gene conversion
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'http://may2017.archive.ensembl.org')
assign_transcripts_to_genes <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = ensembl)

gene_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = ensembl)


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

#Prepare for DESeq2
rownames(desc_table) <- desc_table$sample
if(! all(colnames(txi.kallisto.gene)%in%rownames(desc_table))){stop("Missing some samples in the Kallisto sample file")}

# CHECKED


### FROM HERE PROCEED WITH EXTREME CAUTION:
### NEED TO LEARN MORE CONCERNING INTERACTION TERMS
desc_table$name<-NULL
desc_table$condition<-NULL

desc_table$name[desc_table$sample=="WTCHG_6_205142_kallisto"]<-"HEK293_chr_1"
desc_table$condition[desc_table$sample=="WTCHG_6_205142_kallisto"]<-"HEK293"
desc_table$name[desc_table$sample=="WTCHG_6_206130_kallisto"]<-"HEK293_chr_2"
desc_table$condition[desc_table$sample=="WTCHG_6_206130_kallisto"]<-"HEK293"
desc_table$name[desc_table$sample=="WTCHG_6_207118_kallisto"]<-"shDrosha_chr_3"
desc_table$condition[desc_table$sample=="WTCHG_6_207118_kallisto"]<-"shDrosha"
desc_table$name[desc_table$sample=="WTCHG_6_208106_kallisto"]<-"shDrosha_chr_4"
desc_table$condition[desc_table$sample=="WTCHG_6_208106_kallisto"]<-"shDrosha"
desc_table$name[desc_table$sample=="WTCHG_6_209189_kallisto"]<-"shDicer_chr_5"
desc_table$condition[desc_table$sample=="WTCHG_6_209189_kallisto"]<-"shDicer"
desc_table$name[desc_table$sample=="WTCHG_6_210177_kallisto"]<-"shDicer_chr_6"
desc_table$condition[desc_table$sample=="WTCHG_6_210177_kallisto"]<-"shDicer"
desc_table$name[desc_table$sample=="WTCHG_6_211165_kallisto"]<-"shAgo2_chr_7"
desc_table$condition[desc_table$sample=="WTCHG_6_211165_kallisto"]<-"shAgo2"
desc_table$name[desc_table$sample=="WTCHG_6_212153_kallisto"]<-"shAgo2_chr_8"
desc_table$condition[desc_table$sample=="WTCHG_6_212153_kallisto"]<-"shAgo2"

desc_table$cell_type <- as.factor(rep("HEK293T",8))
names(txi.kallisto.gene)

dds <- DESeqDataSetFromTximport(txi.kallisto.gene, desc_table, ~condition) 

dds <- DESeq(dds)

resultsNames(dds)

'%ino%' <- function(x, table) {
  xSeq <- seq(along = x)
  names(xSeq) <- x
  Out <- xSeq[as.character(table)]
  Out[!is.na(Out)]
}

HEK_dicer_results <- results(dds, contrast=c("condition","shDicer","HEK293"))
dim(HEK_dicer_results)

HEK_dicer_results$symbol<-gene_symbol[gene_symbol$ensembl_gene_id%ino%rownames(HEK_dicer_results),"hgnc_symbol"]
View(HEK_dicer_results)
write.table(HEK_dicer_results[order(HEK_dicer_results$padj),],"~/Desktop/Dicer.txt",sep="\t", quote=F)


HEK_drosha_results <- results(dds, contrast=c("condition","shDrosha","HEK293"))
HEK_drosha_results$symbol<-gene_symbol[gene_symbol$ensembl_gene_id%ino%rownames(HEK_dicer_results),"hgnc_symbol"]
View(HEK_drosha_results)
write.table(HEK_drosha_results[order(HEK_drosha_results$padj),],"~/Desktop/Drosha.txt",sep="\t", quote=F)


HEK_ago2_results <- results(dds, contrast=c("condition","shAgo2","HEK293"))
HEK_ago2_results$symbol<-gene_symbol[gene_symbol$ensembl_gene_id%ino%rownames(HEK_dicer_results),"hgnc_symbol"]
View(HEK_ago2_results)
write.table(HEK_ago2_results[order(HEK_ago2_results$padj),],"~/Desktop/Ago2.txt",sep="\t", quote=F)


HEK_dicer_results$padj[is.na(HEK_dicer_results$padj)]<-1
HEK_ago2_results$padj[is.na(HEK_ago2_results$padj)]<-1
HEK_drosha_results$padj[is.na(HEK_drosha_results$padj)]<-1

dim(HEK_dicer_results[HEK_dicer_results$log2FoldChange>0 &HEK_dicer_results$padj<0.005,])
dim(HEK_ago2_results[HEK_ago2_results$log2FoldChange>0 &HEK_ago2_results$padj<0.005,])
dim(HEK_drosha_results[HEK_drosha_results$log2FoldChange>0 &HEK_drosha_results$padj<0.05,])

dim(HEK_dicer_results[HEK_dicer_results$log2FoldChange>0 &HEK_dicer_results$padj<0.005 &HEK_ago2_results$log2FoldChange>0 &HEK_ago2_results$padj<0.005,])
dim(HEK_ago2_results[HEK_ago2_results$log2FoldChange>0 &HEK_ago2_results$padj<0.005 & HEK_drosha_results$log2FoldChange>0 &HEK_drosha_results$padj<0.05 ,])
dim(HEK_dicer_results[HEK_dicer_results$log2FoldChange>0 &HEK_dicer_results$padj<0.005 & HEK_drosha_results$log2FoldChange>0 &HEK_drosha_results$padj<0.05,])

dim(HEK_dicer_results[HEK_dicer_results$log2FoldChange>0 &HEK_dicer_results$padj<0.005 &HEK_ago2_results$log2FoldChange>0 &HEK_ago2_results$padj<0.005&HEK_drosha_results$padj>0.05 ,])


HEK_dicer_results1<-HEK_dicer_results[HEK_dicer_results$log2FoldChange>0&HEK_ago2_results$log2FoldChange>0&HEK_dicer_results$padj<0.001&HEK_ago2_results$padj<0.001&HEK_drosha_results$padj>0.05,]
write.table(HEK_dicer_results1[order(HEK_dicer_results1$padj),],"~/Desktop/HEK_dicer_results1_dicer001_ago001_drosha05_upregulated.txt",sep="\t",quote=F)
HEK_dicer_results2<-HEK_dicer_results[HEK_dicer_results$log2FoldChange>0&HEK_ago2_results$log2FoldChange>0&HEK_dicer_results$padj<0.005&HEK_ago2_results$padj<0.005&HEK_drosha_results$padj>0.05,]
HEK_dicer_results3<-HEK_dicer_results2[!(rownames(HEK_dicer_results2)%in%rownames(HEK_dicer_results1)),]

write.table(HEK_dicer_results2,"~/Desktop/trna_mappings_ensembl89/HEK_dicer_results1_dicer005_ago005_drosha05_upregulated.txt",sep="\t",quote=F)

dim(HEK_ago2_results[HEK_ago2_results$padj<0.001,])
dim(HEK_ago2_results[HEK_ago2_results$padj<0.005,])
dim(HEK_dicer_results[HEK_dicer_results$padj<0.001,])
dim(HEK_dicer_results[HEK_dicer_results$padj<0.005,])

dim(HEK_dicer_results1)
is.na(HEK_drosha_results)
HEK_dicer_results$log2FoldChange>0
is.na(HEK_drosha_results$padj)
View(HEK_drosha_results[is.na(HEK_drosha_results$padj),])


HEK_dicer_results2<-HEK_dicer_results[HEK_dicer_results$log2FoldChange>0&HEK_ago2_results$log2FoldChange>0&HEK_dicer_results$padj<0.001&HEK_ago2_results$padj<0.001,]
View(HEK_dicer_results2)
write.table(HEK_dicer_results2[order(HEK_dicer_results2$padj),],"~/Desktop/HEK_dicer_results1_dicer001_ago001_droshaNA_upregulated.txt",sep="\t",quote=F)
View(HEK_dicer_results2)



seq = getSequence(id=row.names(HEK_dicer_results1), type="ensembl_gene_id", seqType="gene_exon_intron", mart = ensembl)
exportFASTA(seq,file="HEK_dicer_results1_dicer001_ago001_drosha05.fasta")
dim(seq)
length(row.names(HEK_dicer_results3))
g = getBM(c("hgnc_symbol","chromosome_name","strand","start_position","end_position","ensembl_gene_id"),"ensembl_gene_id", row.names(HEK_dicer_results2) , ensembl)
dim(g)
write.table(g, "~/Desktop/trna_mappings_ensembl89/HEK_dicer_results1_dicer005_ago005_drosha05_coord.txt",  quote=FALSE, sep="\t")



HEK_dicer_results[,HEK_dicer_results$padj<0.001]
tab<-read.delim("~/Desktop/ONC_TSG.tab", sep ="\t",header=T)
ONC<-as.character(tab$Gene[tab$Type=="ONC"])
TSG<-as.character(tab$Gene[tab$Type=="TSG"])

View(HEK_dicer_results[HEK_dicer_results$symbol%ino%ONC,])
View(HEK_drosha_results[HEK_drosha_results$symbol%ino%ONC,])
View(HEK_ago2_results[HEK_ago2_results$symbol%ino%ONC,])


View(HEK_dicer_results[HEK_dicer_results$symbol%ino%TSG,])
View(HEK_drosha_results[HEK_drosha_results$symbol%ino%TSG,])
View(HEK_ago2_results[HEK_ago2_results$symbol%ino%TSG,])

HEK_dicer_results["ENSG00000235501",]
HEK_drosha_results["ENSG00000235501",]
HEK_ago2_results["ENSG00000235501",]

sh_co_results
hist(sh_co_results$pval)
(sh_co_results <- sh_co_results[order(sh_co_results$padj),])
trna_prev_vector<-c("ENST00000262776","ENST00000371655","ENST00000380486","ENST00000361566","ENST00000244519","ENST00000372040","ENST00000259526","ENST00000378943","ENST00000589390","ENST00000316403","ENST00000339364","ENST00000338702","ENST00000267842","ENST00000263816","ENST00000526355","ENST00000396625","ENST00000380861","ENST00000581347","ENST00000370819","ENST00000452846","ENST00000525985","ENST00000227507","ENST00000377113","ENST00000315901","ENST00000591639","ENST00000344051","ENST00000272223")
trna_prev_trans_genes<-assign_transcripts_to_genes[assign_transcripts_to_genes$ensembl_transcript_id%in%vector,]
trna_prev_trans_genes$ensembl_gene_id
View(sh_co_results[trna_prev_trans_genes$ensembl_gene_id,])
(sh_co_results_rest<-subset(sh_co_results, padj < 10^(-3)))
seq = getSequence(id=row.names(sh_co_results_rest), type="ensembl_gene_id", seqType="gene_exon_intron", mart = ensembl)
exportFASTA(seq,file="sh_co_results_rest_sequence.fasta")

row.names(sh_co_results_rest)
g = getBM(c("hgnc_symbol","chromosome_name","strand","start_position","end_position","ensembl_gene_id"),"ensembl_gene_id", row.names(sh_co_results_rest) , ensembl)
show(g)
write.table(g, "~/Desktop/sh_co_results_mat_coordinates.txt",  quote=FALSE, sep="\t")

write.table(sh_co_results, "~/Desktop/sh_co_results_mat.txt",  quote=FALSE, sep="\t")
