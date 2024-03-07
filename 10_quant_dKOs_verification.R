library(biomaRt)
library(DESeq2)
library(limma)
library(ggplot2)
library(ComplexHeatmap)

### read in counts ####
files <- list.files("/cellfile/datapublic/ftitztei/quant_seq_prep_SC/counts/")

count_df <- data.frame(gene=factor())


for(f in files){

  counts <- data.frame(read.table(paste("/cellfile/datapublic/ftitztei/quant_seq_prep_SC/counts/",f,sep="")))
  counts <- counts[,]
  colnames(counts) <- c("gene", strsplit(f,split="_",fixed=TRUE)[[1]][2])

  count_df <- merge.data.frame(count_df, counts, by = "gene", all = TRUE)

}

row.names(count_df) <- count_df$gene
count_df <- count_df[c(6:dim(count_df)[1]),c(2:dim(count_df)[2])]

head(count_df)

### sort out by biotype ####
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
    final_count_df[d,78] <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2]
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])>1){
    mart_types <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2]
    if(length(unique(mart_types))==1){
      final_count_df[d,78] <- mart_types[1]
    }
    else{
      print(d)
    }
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),2])==0){
    final_count_df[d,78] <- NA
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),3])==TRUE){
    final_count_df[d,79] <- biomart[which(biomart$ensembl_gene_id==rownames(final_count_df)[d]),3]
  }
  else{
    final_count_df[d,79] <- NA
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
short_count_df$`114222_2i` <- short_count_df$`114222`
short_count_df <- short_count_df[,c(1:77,80,78,79)]
saveRDS(short_count_df, "RDS/short_counts_quantseq.RDS")


cpm_df <- short_count_df[,1:78]
cpm_df <- as.data.frame(sapply(cpm_df,FUN=as.numeric))
row.names(cpm_df) <- row.names(short_count_df[,1:78])
for(d in 1:dim(cpm_df)[2]){
  cpm_df[,d] <- cpm_df[,d]/(colSums(cpm_df)[d]/1000000)
}

### differenial expression needed for P10 dKO ####
## do one deseq design combining genotype and time to extract FCs
## due to incomplete design (P10 only in one dKO)

design_P10 <- read.csv("data/samplesheet_quantseq_KO1_KO2_dKO.csv")

epressed_names_full <- row.names(cpm_df[which(rowMax(as.matrix(cpm_df))>=0.1),])

ddsFullCountTable_P10 <- DESeqDataSetFromMatrix(countData = short_count_df[epressed_names_full,
                                                                           as.character(design_P10$CCG)],
                                                 colData = design_P10,
                                                 design= ~ condition)

dds_P10 <- estimateSizeFactors(ddsFullCountTable_P10)
normalized_counts_P10 <- counts(dds_P10, normalized=TRUE)

dds_P10$condition <- relevel(dds_P10$condition,"WT_2i")
dds_P10_deseq <- DESeq(dds_P10)
resultsNames(dds_P10_deseq)

## WT
diff_WT_N16vsWT_2i <- data.frame(results(dds_P10_deseq, name="condition_WT_N16_vs_WT_2i"))
diff_WT_N16vsWT_2i$mgi <- final_count_df[row.names(diff_WT_N16vsWT_2i),79]
sig_diff_WT_N16vsWT_2i <- diff_WT_N16vsWT_2i[which(diff_WT_N16vsWT_2i$padj<0.05 & abs(diff_WT_N16vsWT_2i$log2FoldChange)>0.5),]

diff_KO1_KO2_2ivsWT_2i <- data.frame(results(dds_P10_deseq, name="condition_KO1_KO2_2i_vs_WT_2i"))
diff_KO1_KO2_2ivsWT_2i$mgi <- final_count_df[row.names(diff_KO1_KO2_2ivsWT_2i),79]
sig_KO1_KO2_2ivsWT_2i <- diff_KO1_KO2_2ivsWT_2i[which(diff_KO1_KO2_2ivsWT_2i$padj<0.05 &
                                                                      abs(diff_KO1_KO2_2ivsWT_2i$log2FoldChange)>0.5),]

diff_KO1_KO2_N16vsKO1_KO2_2i <- data.frame(results(dds_P10_deseq, contrast = c("condition","KO1_KO2_N16",
                                                                                           "KO1_KO2_2i")))
diff_KO1_KO2_N16vsKO1_KO2_2i$mgi <- final_count_df[row.names(diff_KO1_KO2_N16vsKO1_KO2_2i),79]
sig_KO1_KO2_N16vsKO1_KO2_2i <- diff_KO1_KO2_N16vsKO1_KO2_2i[which(diff_KO1_KO2_N16vsKO1_KO2_2i$padj<0.05 &
                                                                      abs(diff_KO1_KO2_N16vsKO1_KO2_2i$log2FoldChange)>0.5),]

diff_KO1_KO2_P10vsKO1_KO2_2i <- data.frame(results(dds_P10_deseq, contrast = c("condition","KO1_KO2_P10_N16",
                                                                                           "KO1_KO2_2i")))
diff_KO1_KO2_P10vsKO1_KO2_2i$mgi <- final_count_df[row.names(diff_KO1_KO2_P10vsKO1_KO2_2i),79]
sig_KO1_KO2_P10vsKO1_KO2_2i <- diff_KO1_KO2_P10vsKO1_KO2_2i[which(diff_KO1_KO2_P10vsKO1_KO2_2i$padj<0.05 &
                                                                                            abs(diff_KO1_KO2_P10vsKO1_KO2_2i$log2FoldChange)>0.5),]

diff_KO1_KO2_P10vsKO1_KO2_N16 <- data.frame(results(dds_P10_deseq, contrast = c("condition","KO1_KO2_P10_N16",
                                                                                           "KO1_KO2_N16")))
diff_KO1_KO2_P10vsKO1_KO2_N16$mgi <- final_count_df[row.names(diff_KO1_KO2_P10vsKO1_KO2_N16),79]
sig_KO1_KO2_P10vsKO1_KO2_N16 <- diff_KO1_KO2_P10vsKO1_KO2_N16[which(diff_KO1_KO2_P10vsKO1_KO2_N16$padj<0.05 &
                                                                                            abs(diff_KO1_KO2_P10vsKO1_KO2_N16$log2FoldChange)>0.5),]

diff_KO1_KO2_N16vsWT_2i <- data.frame(results(dds_P10_deseq, name="condition_KO1_KO2_N16_vs_WT_2i"))
diff_KO1_KO2_N16vsWT_2i$mgi <- final_count_df[row.names(diff_KO1_KO2_N16vsWT_2i),79]
sig_KO1_KO2_N16vsWT_2i <- diff_KO1_KO2_N16vsWT_2i[which(diff_KO1_KO2_N16vsWT_2i$padj<0.05 &
                                                                      abs(diff_KO1_KO2_N16vsWT_2i$log2FoldChange)>0.5),]

diff_KO1_KO2_P10vsWT_2i <- data.frame(results(dds_P10_deseq, name="condition_KO1_KO2_P10_N16_vs_WT_2i"))
diff_KO1_KO2_P10vsWT_2i$mgi <- final_count_df[row.names(diff_KO1_KO2_P10vsWT_2i),79]
sig_KO1_KO2_P10vsWT_2i <- diff_KO1_KO2_P10vsWT_2i[which(diff_KO1_KO2_P10vsWT_2i$padj<0.05 &
                                                                      abs(diff_KO1_KO2_P10vsWT_2i$log2FoldChange)>0.5),]


naive_formative_frame <- data.frame(row.names = row.names(diff_WT_N16vsWT_2i[which(diff_WT_N16vsWT_2i$mgi%in%
                                                                                     c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                       "Tbx3", "Klf4", "Prdm14",
                                                                                       "Fgf5", "Pou3f1", "Otx2",
                                                                                       "Dnmt3a", "Dnmt3b")),]),
                                    FC_WT_N16_vs_WT_2i=diff_WT_N16vsWT_2i[which(diff_WT_N16vsWT_2i$mgi%in%
                                                                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                    "Tbx3", "Klf4", "Prdm14",
                                                                                    "Fgf5", "Pou3f1", "Otx2",
                                                                                    "Dnmt3a", "Dnmt3b")),2],
                                    padj_WT_N16_vs_WT_2i=diff_WT_N16vsWT_2i[which(diff_WT_N16vsWT_2i$mgi%in%
                                                                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                    "Tbx3", "Klf4", "Prdm14",
                                                                                    "Fgf5", "Pou3f1", "Otx2",
                                                                                    "Dnmt3a", "Dnmt3b")),6],
                                    FC_dKO_2i_vs_WT_2i=diff_KO1_KO2_2ivsWT_2i[which(diff_KO1_KO2_2ivsWT_2i$mgi%in%
                                                                                                                c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                                  "Tbx3", "Klf4", "Prdm14",
                                                                                                                  "Fgf5", "Pou3f1", "Otx2",
                                                                                                                  "Dnmt3a", "Dnmt3b")),2],
                                    padj_dKO_2i_vs_WT_2i=diff_KO1_KO2_2ivsWT_2i[which(diff_KO1_KO2_2ivsWT_2i$mgi%in%
                                                                                                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                                    "Tbx3", "Klf4", "Prdm14",
                                                                                                                    "Fgf5", "Pou3f1", "Otx2",
                                                                                                                    "Dnmt3a", "Dnmt3b")),6],
                                    FC_dKO_N16_vs_dKO_2i=diff_KO1_KO2_N16vsKO1_KO2_2i[which(diff_KO1_KO2_N16vsKO1_KO2_2i$mgi%in%
                                                                                                 c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                   "Tbx3", "Klf4", "Prdm14",
                                                                                                   "Fgf5", "Pou3f1", "Otx2",
                                                                                                   "Dnmt3a", "Dnmt3b")),2],
                                    padj_dKO_N16_vs_dKO_2i=diff_KO1_KO2_N16vsKO1_KO2_2i[which(diff_KO1_KO2_N16vsKO1_KO2_2i$mgi%in%
                                                                                                   c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                     "Tbx3", "Klf4", "Prdm14",
                                                                                                     "Fgf5", "Pou3f1", "Otx2",
                                                                                                     "Dnmt3a", "Dnmt3b")),6],
                                    FC_dKO_P10_vs_dKO_2i=diff_KO1_KO2_P10vsKO1_KO2_2i[which(diff_KO1_KO2_P10vsKO1_KO2_2i$mgi%in%
                                                                                                 c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                   "Tbx3", "Klf4", "Prdm14",
                                                                                                   "Fgf5", "Pou3f1", "Otx2",
                                                                                                   "Dnmt3a", "Dnmt3b")),2],
                                    padj_dKO_P10_vs_dKO_2i=diff_KO1_KO2_P10vsKO1_KO2_2i[which(diff_KO1_KO2_P10vsKO1_KO2_2i$mgi%in%
                                                                                                   c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                     "Tbx3", "Klf4", "Prdm14",
                                                                                                     "Fgf5", "Pou3f1", "Otx2",
                                                                                                     "Dnmt3a", "Dnmt3b")),6],

                                    FC_dKO_N16_vs_WT_2i=diff_KO1_KO2_N16vsWT_2i[which(diff_KO1_KO2_N16vsWT_2i$mgi%in%
                                                                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                    "Tbx3", "Klf4", "Prdm14",
                                                                                    "Fgf5", "Pou3f1", "Otx2",
                                                                                    "Dnmt3a", "Dnmt3b")),2],
                                    padj_dKO_N16_vs_WT_2i=diff_KO1_KO2_N16vsWT_2i[which(diff_KO1_KO2_N16vsWT_2i$mgi%in%
                                                                                    c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                      "Tbx3", "Klf4", "Prdm14",
                                                                                      "Fgf5", "Pou3f1", "Otx2",
                                                                                      "Dnmt3a", "Dnmt3b")),6],
                                    FC_dKO_P10_vs_WT_2i=diff_KO1_KO2_P10vsWT_2i[which(diff_KO1_KO2_P10vsWT_2i$mgi%in%
                                                                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                    "Tbx3", "Klf4", "Prdm14",
                                                                                    "Fgf5", "Pou3f1", "Otx2",
                                                                                    "Dnmt3a", "Dnmt3b")),2],
                                    padj_dKO_P10_vs_WT_2i=diff_KO1_KO2_P10vsWT_2i[which(diff_KO1_KO2_P10vsWT_2i$mgi%in%
                                                                                    c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                      "Tbx3", "Klf4", "Prdm14",
                                                                                      "Fgf5", "Pou3f1", "Otx2",
                                                                                      "Dnmt3a", "Dnmt3b")),6],

                                    FC_dKO_P10_vs_dKO_N16=diff_KO1_KO2_P10vsKO1_KO2_N16[which(diff_KO1_KO2_P10vsKO1_KO2_N16$mgi%in%
                                                                                                 c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                   "Tbx3", "Klf4", "Prdm14",
                                                                                                   "Fgf5", "Pou3f1", "Otx2",
                                                                                                   "Dnmt3a", "Dnmt3b")),2],
                                    padj_dKO_P10_vs_dKO_N16=diff_KO1_KO2_P10vsKO1_KO2_N16[which(diff_KO1_KO2_P10vsKO1_KO2_N16$mgi%in%
                                                                                                   c("Nanog", "Tfcp2l1", "Esrrb",
                                                                                                     "Tbx3", "Klf4", "Prdm14",
                                                                                                     "Fgf5", "Pou3f1", "Otx2",
                                                                                                     "Dnmt3a", "Dnmt3b")),6])

naive_formative_frame$mgi <- NA
for(d in 1:nrow(naive_formative_frame)){
  naive_formative_frame[d,15] <- final_count_df[row.names(naive_formative_frame)[d],79]
}

naive_formative_frame <- naive_formative_frame[c(2,7,5,3,1,10,9,11,6,4,8),]
row.names(naive_formative_frame) <- naive_formative_frame$mgi


star_mat <- naive_formative_frame[,c(1,3,5,7,9,11,13)]
star_mat[1:nrow(star_mat),1:ncol(star_mat)] <- ""
for(r in 1:nrow(star_mat)){
  for(c in 1:ncol(star_mat)){
    if(abs(naive_formative_frame[r,2*c-1]) >= 0.5 &
       naive_formative_frame[r,2*c] <= 0.05){
      star_mat[r,c] <- "*"
    }
  }
}


Heatmap(naive_formative_frame[,c(1,3,5,7,9,11,13)],
        show_row_dend = T, show_row_names = T,
        cluster_rows = T,
        cluster_columns = F,
        col=colorRamp2(c(-3,0,3),c("blue","white","red")),
        name="log2FC vs 2i",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(star_mat[i, j]=="*")
            grid.text(star_mat[i, j], x, y)
        })

### read in MC data to compare PD to dKO response ####
medium_change_results <- readRDS("RDS/MC_DESeq2_results.rds")
sum(row.names(diff_KO1_KO2_P10vsKO1_KO2_2i)!=
      row.names(diff_KO1_KO2_N16vsKO1_KO2_2i))

gg_dKO_MC <- cbind(diff_KO1_KO2_N16vsWT_2i[,c(2,6)],
                   diff_KO1_KO2_P10vsWT_2i[,c(2,6,7)])
colnames(gg_dKO_MC) <- c("FC_dKO_N16_vs_WT_2i",
                         "padj_dKO_N16_vs_WT_2i",
                         "FC_dKO_P10_vs_WT_2i",
                         "padj_dKO_P10_vs_WT_2i",
                         "mgi")
gg_dKO_MC <- gg_dKO_MC[which(row.names(gg_dKO_MC)%in%row.names(medium_change_results[[1]])),]
gg_dKO_MC <- gg_dKO_MC[which(!is.na(gg_dKO_MC$FC_dKO_N16_vs_WT_2i)),]
gg_dKO_MC <- gg_dKO_MC[which(!is.na(gg_dKO_MC$FC_dKO_P10_vs_WT_2i)),]

names_MC <- row.names(medium_change_results[["log2FCs"]])[which(row.names(medium_change_results[["log2FCs"]])%in%row.names(gg_dKO_MC))]

gg_dKO_MC <- cbind(gg_dKO_MC, medium_change_results[["log2FCs"]][names_MC, c(3:8)])
colnames(gg_dKO_MC)[6:11] <- paste0("FC_",colnames(gg_dKO_MC)[6:11])

gg_dKO_MC <- cbind(gg_dKO_MC, medium_change_results[["padj"]][names_MC, c(3:8)])
colnames(gg_dKO_MC)[12:17] <- paste0("padj_",colnames(gg_dKO_MC)[12:17])

gg_dKO_MC <- gg_dKO_MC[,c(1:4,6,12,7,13,8,14,9,15,10,16,11,17,5)]

ggplot(gg_dKO_MC, aes(FC_N2_8_vs_2i_0,FC_dKO_N16_vs_WT_2i))+
  geom_point(alpha=0.1)+
  xlab("log2FC MC WT N8 vs 2i")+
  ylab("log2FC dKO N16 vs WT 2i")+
  ggtitle("log2FCs across conditions") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ggplot(gg_dKO_MC, aes(FC_N2_4_vs_2i_0,FC_dKO_N16_vs_WT_2i))+
  geom_point(alpha=0.1)+
  xlab("log2FC MC WT N4 vs 2i")+
  ylab("log2FC dKO N16 vs WT 2i")+
  ggtitle("log2FCs across conditions") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ggplot(gg_dKO_MC, aes(FC_PD_8_vs_2i_0,FC_dKO_N16_vs_WT_2i))+
  geom_point(alpha=0.1)+
  xlab("log2FC MC WT PD8 vs 2i")+
  ylab("log2FC dKO N16 vs WT 2i")+
  ggtitle("log2FCs across conditions") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ggplot(gg_dKO_MC, aes(FC_CH_8_vs_2i_0,FC_dKO_N16_vs_WT_2i))+
  geom_point(alpha=0.1)+
  xlab("log2FC MC WT CH8 vs 2i")+
  ylab("log2FC dKO N16 vs WT 2i")+
  ggtitle("log2FCs across conditions") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ggplot(gg_dKO_MC, aes(FC_CH_8_vs_2i_0,FC_PD_8_vs_2i_0))+
  geom_point(alpha=0.1)+
  xlab("log2FC MC WT CH8 vs 2i")+
  ylab("log2FC MC WT PD8 vs 2i")+
  ggtitle("log2FCs across conditions") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))


gg_dKO_MC <- gg_dKO_MC[which(!gg_dKO_MC$mgi%in%names(which(table(gg_dKO_MC$mgi)>1))),]
row.names(gg_dKO_MC) <- gg_dKO_MC$mgi

naive_formative_dKO_MC <- gg_dKO_MC[which(row.names(gg_dKO_MC)%in%
                                            c("Nanog", "Tfcp2l1", "Esrrb",
                                              "Tbx3", "Klf4", "Prdm14",
                                              "Fgf5", "Pou3f1", "Otx2",
                                              "Dnmt3a", "Dnmt3b")),]#[c(2,7,5,3,1,10,9,11,6,4,8),]


star_mat_dKO_MC <- naive_formative_dKO_MC[,c(1,3,5,7,9,11,13,15)]
star_mat_dKO_MC[1:nrow(star_mat_dKO_MC),1:ncol(star_mat_dKO_MC)] <- ""
for(r in 1:nrow(star_mat_dKO_MC)){
  for(c in 1:ncol(star_mat_dKO_MC)){
    if(abs(naive_formative_dKO_MC[r,2*c-1]) >= 0.5 &
       naive_formative_dKO_MC[r,2*c] <= 0.05){
      star_mat_dKO_MC[r,c] <- "*"
    }
  }
}

## Fig 5.33
Heatmap(naive_formative_dKO_MC[,c(1,3,5,7,9,11,13,15)],
        show_row_dend = T, show_row_names = T,
        cluster_rows = T,
        cluster_columns = F,
        col=colorRamp2(c(-3,0,3),c("blue","white","red")),
        name="log2FC vs 2i",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(star_mat_dKO_MC[i, j]=="*")
            grid.text(star_mat_dKO_MC[i, j], x, y)
        })

library(corrplot)
## Fig 5.32
corrplot(cor(gg_dKO_MC[,c(1,3,5,7,9,11,13,15)]))


clusters_h <- readRDS("RDS/clusters_h.rds")
genes_upstream <- clusters_h[which(clusters_h$cluster%in%c(2,4,5)),2]

formative_ext <-  data.frame(gg_dKO_MC[which(row.names(gg_dKO_MC)%in%c("Fgf5", "Pou3f1", "Otx2",
                                                         "Dnmt3a", "Dnmt3b", genes_upstream)),
                              c(1,3,9,11,13,15, 5,7)])



formative_ext_clust <- hclust(dist(formative_ext), method="ward.D2" )

library(foreach)
library(tidyverse)

rs.withinClusterVariance = function( m, clus, sum_variances=T ) {
  require(foreach)
  withinvar = foreach( c = unique(clus), .combine = c ) %do% {
    centroid = apply(m[,which(clus==c),drop=F],1,mean)
    return( sum( apply(m[,which(clus==c),drop=F],2,function(x) x - centroid)^2 ) )
  }

  return( sum(withinvar) )
}
## test cluster sizes in reagrads to within cluster variance
## on unscales data
WCV_cluster_size = foreach(i = 6:50, .combine=c) %do% {
  rs.withinClusterVariance(t(formative_ext),
                            cutree(formative_ext_clust,i) )
}
## plot change of within cluster variance to find good k
plot(diff(WCV_cluster_size) ~ c(7:50), type="b", xlim=c(0,50),
     xlab="number of clusters (k)", main="clusters of genes",
     ylab=expression( W[k] ~ "-" ~ W[k-1]), pch=17, cex=1.5, col="black" )
abline(v=c(16,25), lwd=2, lty=2)
number_clusters <- 18

formative_ext_clusters <- data.frame(cutree(formative_ext_clust,number_clusters))
formative_ext_clusters$mgi <- row.names(formative_ext_clusters)
colnames(formative_ext_clusters)[1] <- "cluster"

formative_ext_clusters[c("Fgf5", "Pou3f1", "Otx2",
                         "Dnmt3a", "Dnmt3b"),]

Heatmap(formative_ext,
        show_row_dend = F, show_row_names = F,
        cluster_rows = F,
        cluster_columns = F,
        row_split = formative_ext_clusters$cluster,
        row_gap = unit(5, "mm"),
        col=colorRamp2(c(-3,0,3),c("blue","white","red")),
        name="log2FC vs 2i")


CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
CRISP2 <- readRDS("data/RDS_data/CRISP2_Kdm6a_updated_091020.rds")

formative_ext_clusters$WT_log2FC <- NA
for(d in 1:nrow(formative_ext_clusters)){
  formative_ext_clusters[d,3] <- CLIM2$singlecontrasts$fit_KOVKO$coefficients[formative_ext_clusters[d,2],74]
}

View(formative_ext_clusters[which(formative_ext_clusters$cluster%in%c(2,6,7,9,13) &
                                    formative_ext_clusters$WT_log2FC > 0),])

Heatmap(formative_ext[formative_ext_clusters[which(formative_ext_clusters$cluster%in%c(2,6,7,9,13) &
                                       formative_ext_clusters$WT_log2FC > 0),2],],
        show_row_dend = F, show_row_names = F,
        cluster_rows = F,
        cluster_columns = F,
        row_split = formative_ext_clusters[which(formative_ext_clusters$cluster%in%c(2,6,7,9,13) &
                                                   formative_ext_clusters$WT_log2FC > 0),1],
        row_gap = unit(5, "mm"),
        col=colorRamp2(c(-3,0,3),c("blue","white","red")),
        name="log2FC vs 2i")

## try correlation in N1

pred_form <- as.data.frame( t( CRISP2$voom_corrected_merged_N1[ c("Fgf5", "Pou3f1", "Otx2", "Dnmt3a", "Dnmt3b"),] ) )

ext_form_cor = foreach(i = 1:nrow(formative_ext),
               .combine = rbind) %dopar% {
  this = summary( lm( as.numeric(CRISP2$voom_corrected_merged_N1[row.names(formative_ext)[i],])
                      ~ ., data = pred_form ) )
  c(this$coefficients[2:nrow(this$coefficients),4], "rsq"=this$r.squared)
} %>% as.data.frame

row.names(ext_form_cor) <- row.names(formative_ext)
ext_form_cor$mgi <- row.names(ext_form_cor)
ext_form_cor$WT_log2FC <- formative_ext_clusters$WT_log2FC
View(ext_form_cor)

extended_form_KO <- ext_form_cor[which(ext_form_cor$rsq >= 0.65 &
                                         ext_form_cor$WT_log2FC > 0),7]

rev_form_KO <- ext_form_cor[which(ext_form_cor$rsq >= 0.65 &
                                         ext_form_cor$WT_log2FC < 0),7]



### export objects ####
saveRDS(extended_form_KO, "RDS/extended_form_KO.rds")
saveRDS(rev_form_KO, "RDS/rev_form_KO.rds")
save.image("workspaces/quant_dKOs_verification.RData")
