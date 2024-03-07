library(biomaRt)
library(DESeq2)
library(ggplot2)
library(GenomicFeatures)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
#library(edgeR)


setwd("/cellfile/datapublic/ftitztei/Rscripts/TCpackage/")
### load read counts in dataframe ####

files <- list.files("data/counts_medium_change/")

count_df <- data.frame(gene=factor())


for(f in files){

  counts <- data.frame(read.table(paste("data/counts_medium_change/",f,sep="")))
  counts <- counts[c(5:dim(counts)[1]),c(1,4)]
  colnames(counts) <- c("gene", strsplit(f,split="_",fixed=TRUE)[[1]][2])

  count_df <- merge.data.frame(count_df, counts, by = "gene", all = TRUE)

}

row.names(count_df) <- count_df$gene
count_df <- count_df[,c(2:dim(count_df)[2])]

head(count_df)


### sort out count table by biotype ####

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol"),
                 mart=ensembl)


final_count_df <- as.data.frame(count_df)
new_names <- sapply(strsplit(row.names(final_count_df), ".", fixed=T), function(x) x[1])
row.names(final_count_df) <- new_names
head(final_count_df)

final_count_df$biotype <- NA
final_count_df$mgi <- NA
for(d in 1:dim(final_count_df)[1]){
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])==1){
    final_count_df[d,28] <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2]
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])>1){
    mart_types <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2]
    if(length(unique(mart_types))==1){
      final_count_df[d,28] <- mart_types[1]
    }
    else{
      print(d)
    }
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])==0){
    final_count_df[d,28] <- NA
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),3])==TRUE){
    final_count_df[d,29] <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),3]
  }
  else{
    final_count_df[d,29] <- NA
  }
}

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

txdb_length <- makeTxDbFromGFF("data/gencode.vM11.primary_assembly.annotation.gtf",format="gtf")
exons.list.per.gene <- exonsBy(txdb_length,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(x))})
exon_names <- sapply(strsplit(names(exonic.gene.sizes), ".", fixed=T), function(x) x[1])
names(exonic.gene.sizes) <- exon_names

short_exon_sizes <- exonic.gene.sizes[which(names(exonic.gene.sizes)%in%row.names(short_count_df))]

tpm_table <- as.matrix(short_count_df[,c(1:27)])
tpm_table <- tpm_table / as.numeric(short_exon_sizes)
tpm_table <- t(t(tpm_table)*1e6/colSums(tpm_table))
tpm_table <- data.frame(tpm_table)
tpm_table$biotype <- short_count_df$biotype
tpm_table$mgi <- short_count_df$mgi

tpm_summary <- data.frame(gene=row.names(tpm_table),
                          mean_TPM_all=NA,
                          mean_TPM_2i=NA)
for(d in 1:nrow(tpm_table)){
  tpm_summary[d,2] <- rowMeans(tpm_table[d,c(1:27)])
  tpm_summary[d,3] <- rowMeans(tpm_table[d,c(3:5)])
}

### load 100 KO results & get FPKM table ####
CRISP2 <- readRDS("data/RDS_data/CRISP2_Kdm6a_updated_091020.rds")
CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
SB2 <- readRDS("data/RDS_data/SB2_recalc_121020.RDS")

KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
KOvKO_padj <- CLIM2$singlecontrasts$fit_KOVKO$adj.P.Value
KOvWT <- CLIM2$singlecontrasts$fit_N1$coefficients
KOvWT_padj <- CLIM2$singlecontrasts$fit_N1$adj.P.Value


WT_diff <- row.names(KOvKO)[which(abs(KOvKO[,74])>0.5&
                                    KOvKO_padj[,74]<0.05)]


FPKM_mat_short <- data.frame(RC9_N2=rowMeans(CRISP2$fpkm[,which(colnames(CRISP2$fpkm)%in%
                                                                  row.names(CRISP2$design_rall[which(CRISP2$design_rall$gene_x2=="wtN1"),]))]),
                             RC9_2i=rowMeans(CRISP2$fpkm[,which(colnames(CRISP2$fpkm)%in%
                                                                  row.names(CRISP2$design_rall[which(CRISP2$design_rall$gene_x2=="wt2i"),]))]))

FPKM_log <- log2(FPKM_mat_short)
FPKM_log <- FPKM_log[row.names(KOvKO)[which(row.names(KOvKO)%in%row.names(FPKM_log))],]


FPKM_table_log <- matrix(0,nrow=dim(KOvKO)[1],ncol=2)
row.names(FPKM_table_log) <- row.names(KOvKO)
colnames(FPKM_table_log) <- colnames(FPKM_log)

for(d in 1:dim(FPKM_table_log)[1]){
  if(row.names(FPKM_table_log)[d]%in%row.names(KOvKO)){
    FPKM_table_log[d,] <- as.numeric(FPKM_log[which(row.names(FPKM_log)==row.names(FPKM_table_log)[d]),])
  }
}

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
design <- read.csv("data/samplesheet_medium_change.csv")
design$time <- as.factor(design$time)
design$medium <- factor(as.factor(design$medium),
                        levels = levels(as.factor(design$medium))[c(1,3,2,4)])
design$condition <- factor(as.factor(design$condition),
                           levels = levels(as.factor(design$condition))[c(1:3,6,7,4,5,8,9)])
design <- design[c(26,27,1:25),]

ddsFullCountTable_mx <- DESeqDataSetFromMatrix(countData = short_count_df[,c(1:27)],
                                               colData = design,
                                               design= ~ condition)

dds <- estimateSizeFactors(ddsFullCountTable_mx)
normalized_counts <- counts(dds, normalized=TRUE)

dds$condition <- relevel(dds$condition, ref = "2i_0")
dds <- DESeq(dds)

resultsNames(dds)


log2FC_N4 <- log2((rowMeans(normalized_counts[,12:14]+0.1))/(rowMeans(normalized_counts[,3:5]+0.1)))
log2FC_N8 <- log2((rowMeans(normalized_counts[,15:17]+0.1))/(rowMeans(normalized_counts[,3:5]+0.1)))

summary(abs(log2FC_N4))
summary(abs(log2FC_N8))

sum(abs(log2FC_N4)>=0.5)
sum(abs(log2FC_N8)>=0.5)

### get results from DESeq2 ####
all_contrasts <- data.frame(ref=c(rep("2i_0",8),
                                      "2i_4","2i_4","2i_4","2i_4",
                                      "2i_8","2i_8","2i_8",
                                      "N2_4","N2_4","N2_4",
                                      "N2_8","N2_8",
                                      "CH_4","CH_8",
                                      "CH_4","PD_4"),
                            group=c("2i_4","2i_8",
                                    "N2_4","N2_8",
                                    "CH_4","CH_8",
                                    "PD_4","PD_8",
                                    "2i_8","N2_4","CH_4","PD_4",
                                    "N2_8","CH_8","PD_8",
                                    "N2_8","CH_4","PD_4",
                                    "CH_8","PD_8",
                                    "PD_4","PD_8",
                                    "CH_8","PD_8"))
results_deseq <- list()

for(c in 1:nrow(all_contrasts)){
  results_deseq <- c(results_deseq,list(results(dds, contrast = c("condition",as.character(all_contrasts[c,2]),
                                                  as.character(all_contrasts[c,1])))))
}
names(results_deseq) <- paste(all_contrasts$group,all_contrasts$ref, sep="_vs_")

shrunken_lfcs <-list()
for(c in 1:(length(resultsNames(dds))-1)){
  print(resultsNames(dds)[c+1])
  shrunken_lfcs <- c(shrunken_lfcs,list(lfcShrink(dds, coef = resultsNames(dds)[c+1],
                                                  res=results_deseq[[c]],
                                                  type = "apeglm")))
}
names(shrunken_lfcs) <- resultsNames(dds)[2:9]


medium_change_results <- list(log2FCs=matrix(data=NA,nrow = nrow(as.data.frame(results_deseq[[1]])),
                                             ncol = length(names(results_deseq))),
                              padj=matrix(data=NA,nrow = nrow(as.data.frame(results_deseq[[1]])),
                                             ncol = length(names(results_deseq))))

rownames(medium_change_results$log2FCs) <- row.names(as.data.frame(results_deseq[[1]]))
rownames(medium_change_results$padj) <- row.names(as.data.frame(results_deseq[[1]]))

colnames(medium_change_results$log2FCs) <- names(results_deseq)
colnames(medium_change_results$padj) <- names(results_deseq)

for(d in 1:length(results_deseq)){
  medium_change_results$log2FCs[,d] <- as.data.frame(results_deseq[[d]])[,2]
  medium_change_results$padj[,d] <- as.data.frame(results_deseq[[d]])[,6]
}

medium_change_shrunken <- matrix(data=NA,nrow = nrow(as.data.frame(shrunken_lfcs[[1]])),
                                ncol = length(names(shrunken_lfcs)))
rownames(medium_change_shrunken) <- row.names(as.data.frame(shrunken_lfcs[[1]]))
colnames(medium_change_shrunken) <- names(shrunken_lfcs)

for(d in 1:length(shrunken_lfcs)){
  medium_change_shrunken[,d] <- as.data.frame(shrunken_lfcs[[d]])[,2]
}

medium_change_results$log2FCs <- medium_change_results$log2FCs[rowSums(is.na(medium_change_results$padj))<
                                                                 ncol(medium_change_results$padj),]

medium_change_results$padj <- medium_change_results$padj[rowSums(is.na(medium_change_results$padj))<
                                                           ncol(medium_change_results$padj),]

medium_change_shrunken <- medium_change_shrunken[rowSums(is.na(medium_change_shrunken))<
                                                   ncol(medium_change_shrunken),]


summary(abs(medium_change_results$log2FCs[,3]))
summary(abs(medium_change_results$log2FCs[,4]))

sum(abs(medium_change_results$log2FCs[,3])>=0.5)
sum(abs(medium_change_results$log2FCs[,4])>=0.5)

diff_genes <- list(twoi_4_vs_2i_0=NA,
                   twoi_8_vs_2i_0=NA,
                   N2_4_vs_2i_0=NA,
                   N2_8_vs_2i_0=NA,
                   CH_4_vs_2i_0=NA,
                   CH_8_vs_2i_0=NA,
                   PD_4_vs_2i_0=NA,
                   PD_8_vs_2i_0=NA,
                   twoi_8_vs_2i_4=NA,
                   N2_4_vs_2i_4=NA,
                   CH_4_vs_2i_4=NA,
                   PD_4_vs_2i_4=NA,
                   N2_8_vs_2i_8=NA,
                   CH_8_vs_2i_8=NA,
                   PD_8_vs_2i_8=NA,
                   N2_8_vs_N2_4=NA,
                   CH_4_vs_N2_4=NA,
                   PD_4_vs_N2_4=NA,
                   CH_8_vs_N2_8=NA,
                   PD_8_vs_N2_8=NA,
                   PD_4_vs_CH_4=NA,
                   PD_8_vs_CH_8=NA,
                   CH_8_vs_CH_4=NA,
                   PD_8_vs_PD_4=NA)

## 2i 4 vs 2i 0 ## 44
diff_genes$twoi_4_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,1])>=0.5),1] <= 0.05))

## 2i 8 vs 2i 0 ## 32
diff_genes$twoi_8_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,2])>=0.5),2] <= 0.05))

## N2 4 vs 2i 0 ## 2316
diff_genes$N2_4_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,3])>=0.5),3] <= 0.05))

## N2 8 vs 2i 0 ## 2536
diff_genes$N2_8_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,4])>=0.5),4] <= 0.05))

## CH 4 vs 2i 0 ## 1741
diff_genes$CH_4_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,5])>=0.5),5] <= 0.05))

## CH 8 vs 2i 0 ## 1426
diff_genes$CH_8_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,6])>=0.5),6] <= 0.05))

## PD 4 vs 2i 0 ## 1081
diff_genes$PD_4_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,7])>=0.5),7] <= 0.05))

## PD 8 vs 2i 0 ## 1633
diff_genes$PD_8_vs_2i_0 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,8])>=0.5),8] <= 0.05))

## 2i 8 vs 2i 4 ## 27
diff_genes$twoi_8_vs_2i_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,9])>=0.5),9] <= 0.05))

## N2 4 vs 2i 4 ## 2013
diff_genes$N2_4_vs_2i_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,10])>=0.5),10] <= 0.05))

## CH 4 vs 2i 4 ## 1493
diff_genes$CH_4_vs_2i_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,11])>=0.5),11] <= 0.05))

## PD 4 vs 2i 4 ## 841
diff_genes$PD_4_vs_2i_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,12])>=0.5),12] <= 0.05))

## N2 8 vs 2i 8 ## 2242
diff_genes$N2_8_vs_2i_8 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,13])>=0.5),13] <= 0.05))

## CH 8 vs 2i 8 ## 1228
diff_genes$CH_8_vs_2i_8 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,14])>=0.5),14] <= 0.05))

## PD 8 vs 2i 8 ## 1413
diff_genes$PD_8_vs_2i_8 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,15])>=0.5),15] <= 0.05))

## N2 8 vs N2 4 ## 1350
diff_genes$N2_8_vs_N2_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,16])>=0.5),16] <= 0.05))

## CH 4 vs N2 4 ## 1009
diff_genes$CH_4_vs_N2_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,17])>=0.5),17] <= 0.05))

## PD 4 vs N2 4 ## 1858
diff_genes$PD_4_vs_N2_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,18])>=0.5),18] <= 0.05))

## CH 8 vs N2 8 ## 1612
diff_genes$CH_8_vs_N2_8 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,19])>=0.5),19] <= 0.05))

## PD 8 vs N2 8 ## 1999
diff_genes$PD_8_vs_N2_8 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,20])>=0.5),20] <= 0.05))

## PD 4 vs CH 4 ## 2567
diff_genes$PD_4_vs_CH_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,21])>=0.5),21] <= 0.05))

## PD 8 vs CH 8 ## 2979
diff_genes$PD_8_vs_CH_8 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,22])>=0.5),22] <= 0.05))

## CH 8 vs CH 4 ## 1067
diff_genes$CH_8_vs_CH_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,23])>=0.5),23] <= 0.05))

## PD 8 vs PD 4 ## 343
diff_genes$PD_8_vs_PD_4 <- names(which(medium_change_results$padj[which(abs(medium_change_results$log2FCs[,24])>=0.5),24] <= 0.05))




### PCA ####
pca <- prcomp(normalized_counts,scale. = T, center = T)
gg_pca <- data.frame(sample=row.names(pca$rotation),
                     PC1=pca$rotation[,1],
                     PC2=pca$rotation[,2],
                     PC3=pca$rotation[,3],
                     time=NA,
                     medium=NA)
for(r in 1:dim(gg_pca)[1]){
  gg_pca[r,5] <- as.character(design[which(design$sample==gg_pca[r,1]),3])
  gg_pca[r,6] <- as.character(design[which(design$sample==gg_pca[r,1]),2])
}
var_explained <- pca$sdev**2/sum(pca$sdev**2)

ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, alpha=time, color=medium, fill=medium)) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  scale_color_manual(values=c("black","orangered3","steelblue4","orange3")) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))

ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC3,size=4, alpha=time, color=medium, fill=medium)) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  scale_color_manual(values=c("black","orangered3","steelblue4","orange3")) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%"))

ggplot(gg_pca)+
  geom_point(aes(x=PC2,y=PC3,size=4, alpha=time, color=medium, fill=medium)) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  scale_color_manual(values=c("black","orangered3","steelblue4","orange3")) +
  theme_bw() +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%"))


## changes vs 2i 0
#Fig 5.7
upset(fromList(diff_genes[c(1,2,3,4)]),order.by = "freq")
upset(fromList(diff_genes[c(5:8)]),order.by = "freq")

## PD & CH vs N2
upset(fromList(diff_genes[c(17,18)]),order.by = "freq")
upset(fromList(diff_genes[c(19,20)]),order.by = "freq")

## reverse stress
upset(fromList(diff_genes[c(9,16,23,24)]),order.by = "freq")


Heatmap(medium_change_results$log2FCs,
        show_row_dend = F, show_row_names = F,
        cluster_rows = T,
        name="log2FC")


vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("medium","time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData) +
geom_point(aes(x=PC1,y=PC2,size=4, alpha=time, color=medium, fill=medium)) +
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  scale_color_manual(values=c("black","steelblue4","orangered3","orange3")) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(percentVar[1],1),"%"),
       y=paste0("PC2: ",round(percentVar[2],1),"%"))




### compare foldchanges between groups that are significant vs non significant ####

gg_scatter_N8 <- as.data.frame(medium_change_results$log2FCs[diff_genes[["N2_8_vs_2i_0"]],c("N2_4_vs_2i_0",
                                                                               "N2_8_vs_2i_0")])



ggplot(gg_scatter_N8,aes(x=N2_8_vs_2i_0,y=N2_4_vs_2i_0))+
  geom_point() +
  coord_fixed() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="red") +
  labs(x="log2FC N2B27 8h vs 2i 0h",
       y="log2FC N2B27 4h vs 2i 0h") +
  ggtitle("comparing log2FCs for significant in N2B27 8h vs 2i 0h comparison") +
  theme(plot.title=element_text(hjust=0.5))

### how do gene groups behave in gene groups ####
CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")

KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
KOvKO_padj <- CLIM2$singlecontrasts$fit_KOVKO$adj.P.Value


WT_diff <- row.names(KOvKO)[which(abs(KOvKO[,74])>=0.5&
                                    KOvKO_padj[,74] <= 0.05)]
WT_diff_ensembl <- biomart[which(biomart$mgi_symbol%in%c(WT_diff)),1]


NAGs_frame <- read.table("data/results_marker_multiple_regression_updated.tsv",
                         sep = "\t",colClasses = c("character","numeric","numeric",
                                                   "numeric","numeric","numeric",
                                                   "numeric","numeric","numeric",
                                                   "numeric","numeric"),
                         skip = 1)
colnames(NAGs_frame) <- c("gene", "Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                          "Klf4", "Prdm14", "Zfp42", "rsq", "nondet", "nondet_z")
NAGs <- NAGs_frame[which(NAGs_frame$rsq >= 0.65),1]
NAGs_ensembl <- biomart[which(biomart$mgi_symbol%in%c(NAGs)),1]
upNAGs <- NAGs[which(NAGs%in%row.names(KOvKO[which(KOvKO[,74]>0),]))]
upNAGs_ensembl <- biomart[which(biomart$mgi_symbol%in%c(upNAGs)),1]
downNAGs <- NAGs[which(NAGs%in%row.names(KOvKO[which(KOvKO[,74]<0),]))]
downNAGs_ensembl <- biomart[which(biomart$mgi_symbol%in%c(downNAGs)),1]

naive_ensembl <- biomart[which(biomart$mgi_symbol%in%c("Nanog", "Prdm14",
                                                       "Klf4", "Esrrb", "Tfcp2l1", "Tbx3")),1]

formative_ensembl <- biomart[which(biomart$mgi_symbol%in%c("Dnmt3a", "Dnmt3b",
                                                       "Fgf5", "Otx2", "Pou3f1")),1]

colnames(medium_change_shrunken) <- c("2i_4_vs_2i_0", "2i_8_vs_2i_0", "N2_4_vs_2i_0",
                                      "N2_8_vs_2i_0", "CH_4_vs_2i_0", "CH_8_vs_2i_0",
                                      "PD_4_vs_2i_0", "PD_8_vs_2i_0")

Heatmap(medium_change_shrunken[formative_ensembl,],
        show_row_dend = F, show_row_names = F,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c("blue","white","red")),
        name="log2FC formative")

Heatmap(medium_change_shrunken[naive_ensembl,],
        show_row_dend = F, show_row_names = F,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c("blue","white","red")),
        name="log2FC naive")

Heatmap(medium_change_shrunken[which(row.names(medium_change_shrunken)%in%NAGs_ensembl),],
        show_row_dend = F, show_row_names = F,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c("blue","white","red")),
        name="log2FC NAGs")

Heatmap(medium_change_shrunken[which(row.names(medium_change_shrunken)%in%WT_diff_ensembl),],
        show_row_dend = F, show_row_names = F,
        cluster_rows = T,
        col=colorRamp2(c(-2,0,2),c("blue","white","red")),
        name="log2FC WT diff")

### test for equivalence ####
results_deseq_eq <- list()

for(c in 1:nrow(all_contrasts)){
  results_deseq_eq <- c(results_deseq_eq,list(results(dds, contrast = c("condition",as.character(all_contrasts[c,2]),
                                                                  as.character(all_contrasts[c,1])),
                                                      lfcThreshold = 0.5, altHypothesis = "lessAbs")))
}
names(results_deseq_eq) <- paste(all_contrasts$group,all_contrasts$ref, sep="_vs_")


medium_change_results_eq <- list(log2FCs=matrix(data=NA,nrow = nrow(as.data.frame(results_deseq_eq[[1]])),
                                             ncol = length(names(results_deseq_eq))),
                              padj=matrix(data=NA,nrow = nrow(as.data.frame(results_deseq_eq[[1]])),
                                          ncol = length(names(results_deseq_eq))))

rownames(medium_change_results_eq$log2FCs) <- row.names(as.data.frame(results_deseq_eq[[1]]))
rownames(medium_change_results_eq$padj) <- row.names(as.data.frame(results_deseq_eq[[1]]))

colnames(medium_change_results_eq$log2FCs) <- names(results_deseq_eq)
colnames(medium_change_results_eq$padj) <- names(results_deseq_eq)

for(d in 1:length(results_deseq_eq)){
  medium_change_results_eq$log2FCs[,d] <- as.data.frame(results_deseq_eq[[d]])[,2]
  medium_change_results_eq$padj[,d] <- as.data.frame(results_deseq_eq[[d]])[,6]
}


medium_change_results_eq$log2FCs <- medium_change_results_eq$log2FCs[rowSums(is.na(medium_change_results_eq$padj))<
                                                                 ncol(medium_change_results_eq$padj),]

medium_change_results_eq$padj <- medium_change_results_eq$padj[rowSums(is.na(medium_change_results_eq$padj))<
                                                           ncol(medium_change_results_eq$padj),]



eq_genes <- list(twoi_4_vs_2i_0=NA,
                   twoi_8_vs_2i_0=NA,
                   N2_4_vs_2i_0=NA,
                   N2_8_vs_2i_0=NA,
                   CH_4_vs_2i_0=NA,
                   CH_8_vs_2i_0=NA,
                   PD_4_vs_2i_0=NA,
                   PD_8_vs_2i_0=NA,
                   twoi_8_vs_2i_4=NA,
                   N2_4_vs_2i_4=NA,
                   CH_4_vs_2i_4=NA,
                   PD_4_vs_2i_4=NA,
                   N2_8_vs_2i_8=NA,
                   CH_8_vs_2i_8=NA,
                   PD_8_vs_2i_8=NA,
                   N2_8_vs_N2_4=NA,
                   CH_4_vs_N2_4=NA,
                   PD_4_vs_N2_4=NA,
                   CH_8_vs_N2_8=NA,
                   PD_8_vs_N2_8=NA,
                   PD_4_vs_CH_4=NA,
                   PD_8_vs_CH_8=NA,
                   CH_8_vs_CH_4=NA,
                   PD_8_vs_PD_4=NA)

## 2i 4 vs 2i 0 ## 12029
eq_genes$twoi_4_vs_2i_0 <- names(which(medium_change_results_eq$padj[,1] <= 0.05))

## 2i 8 vs 2i 0 ## 12025
eq_genes$twoi_8_vs_2i_0 <- names(which(medium_change_results_eq$padj[,2] <= 0.05))

## N2 4 vs 2i 0 ## 9835
eq_genes$N2_4_vs_2i_0 <- names(which(medium_change_results_eq$padj[,3] <= 0.05))

## N2 8 vs 2i 0 ## 9716
eq_genes$N2_8_vs_2i_0 <- names(which(medium_change_results_eq$padj[,4] <= 0.05))

## CH 4 vs 2i 0 ## 10271
eq_genes$CH_4_vs_2i_0 <- names(which(medium_change_results_eq$padj[,5] <= 0.05))

## CH 8 vs 2i 0 ## 10533
eq_genes$CH_8_vs_2i_0 <- names(which(medium_change_results_eq$padj[,6] <= 0.05))

## PD 4 vs 2i 0 ## 10856
eq_genes$PD_4_vs_2i_0 <- names(which(medium_change_results_eq$padj[,7] <= 0.05))

## PD 8 vs 2i 0 ## 10407
eq_genes$PD_8_vs_2i_0 <- names(which(medium_change_results_eq$padj[,8] <= 0.05))

## 2i 8 vs 2i 4 ## 11888
eq_genes$twoi_8_vs_2i_4 <- names(which(medium_change_results_eq$padj[,9] <= 0.05))

## N2 4 vs 2i 4 ## 10084
eq_genes$N2_4_vs_2i_4 <- names(which(medium_change_results_eq$padj[,10] <= 0.05))

## CH 4 vs 2i 4 ## 10488
eq_genes$CH_4_vs_2i_4 <- names(which(medium_change_results_eq$padj[,11] <= 0.05))

## PD 4 vs 2i 4 ## 11073
eq_genes$PD_4_vs_2i_4 <- names(which(medium_change_results_eq$padj[,12] <= 0.05))

## N2 8 vs 2i 8 ## 9862
eq_genes$N2_8_vs_2i_8 <- names(which(medium_change_results_eq$padj[,13] <= 0.05))

## CH 8 vs 2i 8 ## 10606
eq_genes$CH_8_vs_2i_8 <- names(which(medium_change_results_eq$padj[,14] <= 0.05))

## PD 8 vs 2i 8 ## 10492
eq_genes$PD_8_vs_2i_8 <- names(which(medium_change_results_eq$padj[,15] <= 0.05))

## N2 8 vs N2 4 ## 10481
eq_genes$N2_8_vs_N2_4 <- names(which(medium_change_results_eq$padj[,16] <= 0.05))

## CH 4 vs N2 4 ## 10847
eq_genes$CH_4_vs_N2_4 <- names(which(medium_change_results_eq$padj[,17] <= 0.05))

## PD 4 vs N2 4 ## 10139
eq_genes$PD_4_vs_N2_4 <- names(which(medium_change_results_eq$padj[,18] <= 0.05))

## CH 8 vs N2 8 ## 10417
eq_genes$CH_8_vs_N2_8 <- names(which(medium_change_results_eq$padj[,19] <= 0.05))

## PD 8 vs N2 8 ## 10071
eq_genes$PD_8_vs_N2_8 <- names(which(medium_change_results_eq$padj[,20] <= 0.05))

## PD 4 vs CH 4 ## 9559
eq_genes$PD_4_vs_CH_4 <- names(which(medium_change_results_eq$padj[,21] <= 0.05))

## PD 8 vs CH 8 ## 9320
eq_genes$PD_8_vs_CH_8 <- names(which(medium_change_results_eq$padj[,22] <= 0.05))

## CH 8 vs CH 4 ## 10762
eq_genes$CH_8_vs_CH_4 <- names(which(medium_change_results_eq$padj[,23] <= 0.05))

## PD 8 vs PD 4 ## 11399
eq_genes$PD_8_vs_PD_4 <- names(which(medium_change_results_eq$padj[,24] <= 0.05))

sum(tpm_summary[which(tpm_summary$mean_TPM_2i>=10),1]%in%eq_genes$N2_4_vs_2i_0)

### export rds files in folder for other scripts ####

saveRDS(tpm_table,"RDS/tpms_MC.rds")
saveRDS(medium_change_shrunken,"RDS/log2FCs_shrunken_MC.rds")
saveRDS(medium_change_results,"RDS/MC_DESeq2_results.rds")
saveRDS(diff_genes,"RDS/diff_genes_MC.rds")

saveRDS(medium_change_results_eq,"RDS/MC_DESeq2_results_eq.rds")
saveRDS(eq_genes,"RDS/eq_genes_MC.rds")


### create export tables ####
out_count <- short_count_df
colnames(out_count)[1:27] <- paste(as.character(design$condition),
                                   as.character(design$replicate),
                                  sep="_")
out_count <- out_count[,c(3:20,1:2,21:29)]

write.csv(out_count, "output_tables/MC_counts.csv", quote = F)
write.csv(medium_change_results$log2FCs, "output_tables/MC_log2FCs.csv", quote = F)
write.csv(medium_change_results$padj, "output_tables/MC_padj.csv", quote = F)


### naive formative ####
naive_formative_log2FCs <- medium_change_results[[1]][biomart[which(biomart$mgi_symbol%in%
                                                                      c("Nanog", "Tfcp2l1",
                                                                        "Esrrb", "Tbx3",
                                                                        "Klf4", "Prdm14",
                                                                        "Fgf5", "Pou3f1",
                                                                        "Otx2","Dnmt3a",
                                                                        "Dnmt3b")),1],]

naive_formative_log2FCs <- as.data.frame(naive_formative_log2FCs[c(1,9,2,8,7,3,10,6,4,5,11),])
row.names(naive_formative_log2FCs) <- c("Nanog", "Tfcp2l1",
                                        "Esrrb", "Tbx3",
                                        "Klf4", "Prdm14",
                                        "Fgf5", "Pou3f1",
                                        "Otx2","Dnmt3a",
                                        "Dnmt3b")

naive_formative_padj <- medium_change_results[[2]][biomart[which(biomart$mgi_symbol%in%
                                                                      c("Nanog", "Tfcp2l1",
                                                                        "Esrrb", "Tbx3",
                                                                        "Klf4", "Prdm14",
                                                                        "Fgf5", "Pou3f1",
                                                                        "Otx2","Dnmt3a",
                                                                        "Dnmt3b")),1],]

naive_formative_padj <- as.data.frame(naive_formative_padj[c(1,9,2,8,7,3,10,6,4,5,11),])
row.names(naive_formative_padj) <- c("Nanog", "Tfcp2l1",
                                        "Esrrb", "Tbx3",
                                        "Klf4", "Prdm14",
                                        "Fgf5", "Pou3f1",
                                        "Otx2","Dnmt3a",
                                        "Dnmt3b")

star_mat <- naive_formative_padj
star_mat[1:nrow(star_mat),1:ncol(star_mat)] <- ""
for(r in 1:nrow(star_mat)){
  for(c in 1:ncol(star_mat)){
    if(abs(naive_formative_log2FCs[r,c]) >= 0.5 &
       naive_formative_padj[r,c] <= 0.05){
      star_mat[r,c] <- "*"
    }
  }
}


Heatmap(naive_formative_log2FCs[,1:8],
        show_row_dend = T, show_row_names = T,
        cluster_rows = T,
        col=colorRamp2(c(-3,0,3),c("blue","white","red")),
        name="log2FC vs 2i",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(star_mat[i, j]=="*")
          grid.text(star_mat[i, j], x, y)
        })


