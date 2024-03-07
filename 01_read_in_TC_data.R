### This script is used to read in the time course data from GSE145653
### The raw counts are processed with DESeq2 and normalized
### Shrunken fold changes are calculated as no samples are present for
### timepoints between 0 and 32 hours

## reading in packages that are needed
library(DESeq2)
library(ggplot2)
library(biomaRt)
library(GenomicFeatures)
#library(edgeR)


### load read counts in dataframe ####

## get a list of all files with mapped reads
files <- list.files("data/counts_TC/")
## sort files by time point not by name
files <- files[c(12:15,11,19:21,1:10,16:18)]

## merge counts per gene from each file together in one data.frame
count_df <- data.frame(gene=factor())

for(f in files){
  counts <- data.frame(read.table(paste("data/counts_TC/",f,sep="")))
  counts <- counts[c(5:dim(counts)[1]),c(1,4)]
  colnames(counts) <- c("gene", strsplit(f,split=".",fixed=TRUE)[[1]][1])

  count_df <- merge.data.frame(count_df, counts, by = "gene", all = TRUE)
}

## set gene names as row names and remove gene name column
row.names(count_df) <- count_df$gene
count_df <- count_df[,c(2:dim(count_df)[2])]

head(count_df)

### sort out count table by biotype ####

## retrieve a lookup between ensembl gene ids and mgi symbols from biomart
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol"),
                 mart=ensembl)


## remove version info from gene names
final_count_df <- as.data.frame(count_df)
new_names <- sapply(strsplit(row.names(final_count_df), ".", fixed=T), function(x) x[1])
row.names(final_count_df) <- new_names
head(final_count_df)

## add columns with info from the biomart lookup to the count table
final_count_df$biotype <- NA
final_count_df$mgi <- NA
for(d in 1:dim(final_count_df)[1]){
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])==1){
    final_count_df[d,22] <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2]
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])>1){
    mart_types <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2]
    if(length(unique(mart_types))==1){
      final_count_df[d,22] <- mart_types[1]
    }
    else{
      print(d)
    }
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])==0){
    final_count_df[d,22] <- NA
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),3])==TRUE){
    final_count_df[d,23] <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),3]
  }
  else{
    final_count_df[d,23] <- NA
  }
}

## sort out genes that do not belong to certain biotypes
final_count_df$biotype <- as.factor(final_count_df$biotype)

short_count_df <- final_count_df[which(final_count_df$biotype%in%
                                         c("protein_coding","lincRNA",
                                           "processed_transcript","antisense",
                                           "3prime_overlapping_ncRNA",
                                           "bidirectional_promoter_lncRNA",
                                           "macro_lncRNA","miRNA",
                                           "misc_RNA","lncRNA",
                                           "scaRNA","scRNA",
                                           "sense_intronic","sense_overlapping",
                                           "snoRNA","snRNA","sRNA")),]

head(short_count_df)


### get TPMs ####

## retrieve length of exons to normalize for read length
txdb_length <- makeTxDbFromGFF("data/gencode.vM11.primary_assembly.annotation.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb_length,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(x))})
exon_names <- sapply(strsplit(names(exonic.gene.sizes), ".", fixed=T), function(x) x[1])
names(exonic.gene.sizes) <- exon_names

## only use exon sizes that are in the count table after sorting out by biotype
short_exon_sizes <- exonic.gene.sizes[which(names(exonic.gene.sizes)%in%row.names(short_count_df))]

## calculate TPMs by normalizing for read length and sequencing depth
tpm_table <- as.matrix(short_count_df[,c(1:21)])
tpm_table <- tpm_table / as.numeric(short_exon_sizes)
tpm_table <- t(t(tpm_table)*1e6/colSums(tpm_table))
tpm_table <- data.frame(tpm_table)
tpm_table$biotype <- short_count_df$biotype
tpm_table$mgi <- short_count_df$mgi


### load 100 KO results & get FPKM table ####
## load in object from KO data
CRISP2 <- readRDS("data/RDS_data/CRISP2_Kdm6a_updated_091020.rds")
CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
SB2 <- readRDS("data/RDS_data/SB2_recalc_121020.RDS")

## revtrieve log2FCs and adjusted pvalues for 24 hours vs 0 hours (2i)
## for all KOs and WT
KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
KOvKO_padj <- CLIM2$singlecontrasts$fit_KOVKO$adj.P.Value
## revtrieve log2FCs and adjusted pvalues for KO vs WT
## for all KOs at 24 hours and 0 hours (2i)
KOvWT <- CLIM2$singlecontrasts$fit_N1$coefficients
KOvWT_padj <- CLIM2$singlecontrasts$fit_N1$adj.P.Value

## last column to get WT 24 hours vs 0 hours (2i)
WT_diff <- row.names(KOvKO)[which(abs(KOvKO[,74])>0.5&
                                    KOvKO_padj[,74]<0.05)]

## get FPKM table from KO data
FPKM_mat_short <- data.frame(RC9_N2=rowMeans(CRISP2$fpkm[,which(colnames(CRISP2$fpkm)%in%
                                                                  row.names(CRISP2$design_rall[which(CRISP2$design_rall$gene_x2=="wtN1"),]))]),
                             RC9_2i=rowMeans(CRISP2$fpkm[,which(colnames(CRISP2$fpkm)%in%
                                                                  row.names(CRISP2$design_rall[which(CRISP2$design_rall$gene_x2=="wt2i"),]))]))

## log2 transform FPKM tanle
FPKM_log <- log2(FPKM_mat_short)
FPKM_log <- FPKM_log[row.names(KOvKO)[which(row.names(KOvKO)%in%row.names(FPKM_log))],]

## create a table that only has genes that are also in the WT comparison
FPKM_table_log <- matrix(0,nrow=dim(KOvKO)[1],ncol=2)
row.names(FPKM_table_log) <- row.names(KOvKO)
colnames(FPKM_table_log) <- colnames(FPKM_log)

for(d in 1:dim(FPKM_table_log)[1]){
  if(row.names(FPKM_table_log)[d]%in%row.names(KOvKO)){
    FPKM_table_log[d,] <- as.numeric(FPKM_log[which(row.names(FPKM_log)==row.names(FPKM_table_log)[d]),])
  }
}

## delog the logged table again to have both options and add additional info
FPKM_table <- data.frame(2**FPKM_table_log)
FPKM_table$ensembl <- NA
for(d in 1:dim(FPKM_table)[1]){
  if(length(biomart[which(biomart$mgi_symbol==
                          row.names(FPKM_table)[d]),1])==1){
    FPKM_table[d,3] <- biomart[which(biomart$mgi_symbol==
                                       row.names(FPKM_table)[d]),1]
  }
}


### get normalized read counts from DESeq####

## run DESeq analysis to retrieve shrunken fold changes
design <- read.csv("data/samplesheet.csv")

ddsFullCountTable_mx <- DESeqDataSetFromMatrix(countData = short_count_df[,c(1:21)],
                                               colData = design,
                                               design= ~ time)

dds <- estimateSizeFactors(ddsFullCountTable_mx)
normalized_counts <- counts(dds, normalized=TRUE)

dds$time <- relevel(dds$time, ref = "0h")
dds <- DESeq(dds)

resultsNames(dds)
shrunken_foldchanges <- matrix(NA, nrow=dim(normalized_counts)[1], ncol=18)
row.names(shrunken_foldchanges) <- row.names(normalized_counts)
colnames(shrunken_foldchanges) <- c("RC9_2iL", "RC9_2i", "RC9_2h",
                                    "RC9_4h", "RC9_6h", "RC9_8h",
                                    "RC9_10h", "RC9_12h", "RC9_14h",
                                    "RC9_16h", "RC9_18h", "RC9_20h",
                                    "RC9_22h", "RC9_24h", "RC9_26h",
                                    "RC9_28h", "RC9_30h", "RC9_32h")

## calculate shrunken fold changes for each condition
for(d in 1:dim(shrunken_foldchanges)[2]){
## for 2iLIF
  if(d==1){
    shrunken_foldchanges[,d] <- lfcShrink(dds,coef=2, type = "apeglm")[,2]#,type="apeglm"
  }
## always 0 for 2i
  if(d==2){
    shrunken_foldchanges[,d] <- 0
  }
## for all other time points
  if(d>2){
    shrunken_foldchanges[,d] <- lfcShrink(dds,coef=paste("time_",(d-2)*2,"h_vs_0h",
                                                         sep = ""), type = "apeglm")[,2]#,type="apeglm"
  }
}
shrunken_foldchanges_complete <- na.omit(shrunken_foldchanges)


### PCA ####
## run standard PCA analysis
vsd <- vst(dds)
pcaData <- plotPCA(vsd,"time", returnData=TRUE,ntop=1000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$group <- factor(pcaData$group,
                        levels = levels(pcaData$group)[c(2,1,13,16:18,3:12,14,15)])
## plot first two dimensions
#Fig 5.1
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw()

## select 50 and 100 most variable genes from PC1
rv <- rowVars(assay(vsd))
select35000 <- order(rv,decreasing=TRUE)[seq_len(min(35000,length(rv)))]
pca <- prcomp(t(assay(vsd)[select35000,]))
loadings <- as.data.frame((pca$rotation))
loadings_order <- order(abs(loadings$PC1),decreasing = TRUE)

genes50 <- row.names(loadings)[loadings_order[1:50]]
genes100 <- row.names(loadings)[loadings_order[1:100]]


### log normalized data to 2i ####
## construct table to fill
log_norm_counts <- data.frame(matrix(NA, nrow=dim(normalized_counts)[1], ncol=21))
row.names(log_norm_counts) <- row.names(normalized_counts)
colnames(log_norm_counts) <- c("RC9_2iL_1", "RC9_2iL_2",
                               "RC9_2i_1", "RC9_2i_2", "RC9_2h",
                               "RC9_4h", "RC9_6h", "RC9_8h",
                               "RC9_10h", "RC9_12h", "RC9_14h",
                               "RC9_16h", "RC9_18h", "RC9_20h",
                               "RC9_22h", "RC9_24h", "RC9_26h",
                               "RC9_28h", "RC9_30h",
                               "RC9_32h_1","RC9_32h_2")

for(d in 1:dim(log_norm_counts)[1]){
  ## log2 of counts per time point divided by mean counts in 2i
  ## pseudo count of 1 added to avoid deviding by 0 or log2(0)
  log_norm_counts[d,] <- log2((normalized_counts[d,c(3:4,1:2,5:21)]+1)/(mean(normalized_counts[d,c(1:2)]+1)))
}

log_norm_counts <- as.matrix(log_norm_counts)

log_norm_counts_short <- log_norm_counts[row.names(shrunken_foldchanges_complete),]


### export rds files in folder for other scripts ####

saveRDS(log_norm_counts_short,"RDS/log_norm_counts_short.rds")
saveRDS(shrunken_foldchanges_complete,"RDS/shrunken_foldchanges_complete.rds")
write.csv(shrunken_foldchanges_complete,"output_tables/shrunkenFCs_TC.csv",quote = F)
saveRDS(tpm_table,"RDS/tpm_table.rds")
write.csv(tpm_table,"output_tables/TPMs_TC.csv",quote = F)
saveRDS(FPKM_table,"RDS/FPKM_table.rds")
saveRDS(short_count_df,"RDS/short_count_df.rds")
write.csv(short_count_df,"output_tables/counts_TC.csv",quote = F)
