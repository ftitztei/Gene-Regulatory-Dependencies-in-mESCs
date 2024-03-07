
library(biomaRt)
library(ggplot2)
library(ComplexHeatmap)

## load functions for this script and other scripts in the analysis
source("R/functions.R")


### read in rds and objects from QC and prerocessing ####
CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
gpr_list <- readRDS("RDS/gpr_list_shrunkenFCS.rds")

FPKM_table <- readRDS("RDS/FPKM_table.rds")
tpms_TC_gpr <- readRDS("RDS/tpms_TC_gpr.rds")


## biomart lookup between ensembl and mgi symbols
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")#, host = "uswest.ensembl.org")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol"),
                 mart=ensembl)
biomart <- biomart[which(!is.na(biomart$entrezgene)),]
biomart <- biomart[complete.cases(biomart),]

## revtrieve log2FCs and adjusted pvalues for 24 hours vs 0 hours (2i)
## for all KOs and WT
KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
KOvKO_padj <- CLIM2$singlecontrasts$fit_KOVKO$adj.P.Value
## revtrieve log2FCs and adjusted pvalues for KO vs WT
## for all KOs at 24 hours and 0 hours (2i)
KOvWT <- CLIM2$singlecontrasts$fit_N1$coefficients
KOvWT_padj <- CLIM2$singlecontrasts$fit_N1$adj.P.Value

### find possible process process dependencies ####
## get log2 FCs after Gaussian Process Regression
log_norm_heatmap_gpr <- gpr_list[[1]]
## set everything in relation to 2i
log_norm_heatmap_gpr <- log_norm_heatmap_gpr-log_norm_heatmap_gpr[,2]
## genes that change during the TC are those that have a log2 FC of at least 0.5
## between the maximum expression and the minimum expression
log_norm_heatmap_gpr_changed <- log_norm_heatmap_gpr[which(rowMax(log_norm_heatmap_gpr[,2:18])-
                                                             rowMin(log_norm_heatmap_gpr[,2:18])>=0.5),]
changed_genes_ensembl <- row.names(log_norm_heatmap_gpr_changed)
changed_genes_mgi <- biomart[which(biomart$ensembl_gene_id%in%changed_genes_ensembl),3]


DEGs <- row.names(KOvKO_padj[which(KOvKO_padj[,74] <= 0.05 &
                                     abs(KOvKO[,74]) >= 0.5),])

lookup <- data.frame(names=c("RC9_2iL",
                             "RC9_2i", "RC9_2h",
                             "RC9_4h", "RC9_6h", "RC9_8h",
                             "RC9_10h", "RC9_12h", "RC9_14h",
                             "RC9_16h", "RC9_18h", "RC9_20h",
                             "RC9_22h", "RC9_24h", "RC9_26h",
                             "RC9_28h", "RC9_30h",
                             "RC9_32h"),
                     time=c(-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32))

gpr_TPMs <- readRDS("RDS/tpms_TC_gpr.rds")
gpr_TPMs_short <- data.frame(RC9_2i=rowMeans(gpr_TPMs[,1:2]),
                             RC9_2h=gpr_TPMs[,5],
                             RC9_4h=gpr_TPMs[,6],
                             RC9_6h=gpr_TPMs[,7],
                             RC9_8h=gpr_TPMs[,8],
                             RC9_10h=gpr_TPMs[,9],
                             RC9_12h=gpr_TPMs[,10],
                             RC9_14h=gpr_TPMs[,11],
                             RC9_16h=gpr_TPMs[,12],
                             RC9_18h=gpr_TPMs[,13],
                             RC9_20h=gpr_TPMs[,14],
                             RC9_22h=gpr_TPMs[,15],
                             RC9_24h=gpr_TPMs[,16],
                             RC9_26h=gpr_TPMs[,17],
                             RC9_28h=gpr_TPMs[,18],
                             RC9_30h=gpr_TPMs[,19],
                             RC9_32h=rowMeans(gpr_TPMs[,20:21]),
                             mgi=gpr_TPMs[,23],
                             min_TP=NA,
                             max_TP=NA,
                             type=NA)

### normalize data & compare to TC ####
compare_gpr_WT_diff <- data.frame(WT_diff = KOvKO[which(row.names(KOvKO)%in%DEGs),74],
                                  row.names = row.names(KOvKO[which(row.names(KOvKO)%in%
                                                                      DEGs),]))
compare_gpr_WT_diff$ensembl <- NA
for(d in 1:dim(compare_gpr_WT_diff)[1]){
  if(length(biomart[which(biomart$mgi_symbol==
                          row.names(compare_gpr_WT_diff)[d]),1])==1){
    compare_gpr_WT_diff[d,2] <- biomart[which(biomart$mgi_symbol==
                                                row.names(compare_gpr_WT_diff)[d]),1]
  }
}
compare_gpr_WT_diff <- merge(compare_gpr_WT_diff,data.frame(gpr=log_norm_heatmap_gpr_changed[,14]),
                             all.x=T,all.y=F,by.x="ensembl",by.y="row.names")
compare_gpr_WT_diff$tresh_TC <- 0
compare_gpr_WT_diff$tresh_100KO <- 0
for(d in 1:dim(compare_gpr_WT_diff)[1]){
  if(compare_gpr_WT_diff$ensembl[d]%in%row.names(tpms_TC_gpr)){
    if(tpms_TC_gpr[which(row.names(tpms_TC_gpr)==compare_gpr_WT_diff$ensembl[d]),16]>=1){
      compare_gpr_WT_diff[d,4] <- 1
    }
  }
  if(compare_gpr_WT_diff$ensembl[d]%in%FPKM_table$ensembl){
    if(!is.na(compare_gpr_WT_diff$ensembl[d])){
      if(!is.na(max(FPKM_table[which(FPKM_table$ensembl==compare_gpr_WT_diff$ensembl[d]),1:2])>=1)){
        if(max(FPKM_table[which(FPKM_table$ensembl==compare_gpr_WT_diff$ensembl[d]),1:2])>=1){
          compare_gpr_WT_diff[d,5] <- 1
        }
      }
    }
  }
}
head(compare_gpr_WT_diff)

compare_gpr_WT_diff_short <- compare_gpr_WT_diff
#compare_gpr_WT_diff_short <- compare_gpr_WT_diff[which(rowSums(compare_gpr_WT_diff[,4:5])==2),]
#compare_gpr_WT_diff_short <- compare_gpr_WT_diff_short[which(sign(compare_gpr_WT_diff_short$WT_diff)==
#                                                               sign(compare_gpr_WT_diff_short$gpr)),]
compare_gpr_WT_diff_short <- compare_gpr_WT_diff_short[which(!is.na(compare_gpr_WT_diff_short$gpr)),]
lm(compare_gpr_WT_diff_short$WT_diff~0+compare_gpr_WT_diff_short$gpr)

gene_profiles_WT_diff <- log_norm_heatmap_gpr_changed[which(row.names(log_norm_heatmap_gpr_changed)%in%
                                                              compare_gpr_WT_diff_short$ensembl),]

for(d in 1:dim(gene_profiles_WT_diff)[1]){
  gene_profiles_WT_diff[d,] <- gene_profiles_WT_diff[d,]* 1.312
    #mean(compare_gpr_WT_diff_short[which(compare_gpr_WT_diff_short$ensembl==row.names(gene_profiles_WT_diff)),2]/
    #   compare_gpr_WT_diff_short[which(compare_gpr_WT_diff_short$ensembl==row.names(gene_profiles_WT_diff)),3])
}

gene_profiles_WT_diff_24 <- gene_profiles_WT_diff[,1:18]
for(d in 1:dim(gene_profiles_WT_diff)[1]){
  gene_profiles_WT_diff_24[d,] <- gene_profiles_WT_diff[d,1:18] - gene_profiles_WT_diff[d,14]
}


for(d in 1:dim(gene_profiles_WT_diff_24)[1]){
  row.names(gene_profiles_WT_diff_24)[d] <- biomart[which(biomart$ensembl_gene_id==row.names(gene_profiles_WT_diff_24)[d]),3]
}


KO_wt_mat_all <- KOvWT[which(row.names(KOvWT)%in%row.names(gene_profiles_WT_diff_24)),]
KO_wt_mat_all <- as.matrix(KO_wt_mat_all[match(row.names(gene_profiles_WT_diff_24),row.names(KO_wt_mat_all)),])

distances_KOs_all <- distanceKOvsTP(KO_wt_mat_all,gene_profiles_WT_diff_24)

rel_distances_KOs_all <- getRelativeDistance(distances_KOs_all)


compare_gpr_WT_diff_short_naive <- compare_gpr_WT_diff_short[which(compare_gpr_WT_diff_short$ensembl%in%
                                                                     c("ENSMUSG00000012396","ENSMUSG00000021255",
                                                                       "ENSMUSG00000042414","ENSMUSG00000026380",
                                                                       "ENSMUSG00000003032","ENSMUSG00000018604")),]
lm(compare_gpr_WT_diff_short_naive$WT_diff~0+compare_gpr_WT_diff_short_naive$gpr)

gene_profiles_naive <- log_norm_heatmap_gpr_changed[which(row.names(log_norm_heatmap_gpr_changed)%in%
                                                            compare_gpr_WT_diff_short_naive$ensembl),]

for(d in 1:dim(gene_profiles_naive)[1]){
  gene_profiles_naive[d,] <- gene_profiles_naive[d,]* 1.454
}

gene_profiles_naive_24 <- gene_profiles_naive[,1:18]
for(d in 1:dim(gene_profiles_naive)[1]){
  gene_profiles_naive_24[d,] <- gene_profiles_naive[d,1:18] - gene_profiles_naive[d,14]
}


for(d in 1:dim(gene_profiles_naive_24)[1]){
  row.names(gene_profiles_naive_24)[d] <- biomart[which(biomart$ensembl_gene_id==row.names(gene_profiles_naive_24)[d]),3]
}


KO_wt_mat_naive <- KOvWT[which(row.names(KOvWT)%in%row.names(gene_profiles_naive_24)),]
KO_wt_mat_naive <- as.matrix(KO_wt_mat_naive[match(row.names(gene_profiles_naive_24),row.names(KO_wt_mat_naive)),])

distances_KOs_naive <- distanceKOvsTP(KO_wt_mat_naive,gene_profiles_naive_24)

rel_distances_KOs_naive <- getRelativeDistance(distances_KOs_naive)



compare_gpr_WT_diff_short_formative <- compare_gpr_WT_diff_short[which(compare_gpr_WT_diff_short$ensembl%in%
                                                                     c("ENSMUSG00000021848", "ENSMUSG00000020661",
                                                                       "ENSMUSG00000090125", "ENSMUSG00000029337",
                                                                       "ENSMUSG00000027478")),]
lm(compare_gpr_WT_diff_short_formative$WT_diff~0+compare_gpr_WT_diff_short_formative$gpr)

gene_profiles_formative <- log_norm_heatmap_gpr_changed[which(row.names(log_norm_heatmap_gpr_changed)%in%
                                                            compare_gpr_WT_diff_short_formative$ensembl),]

for(d in 1:dim(gene_profiles_formative)[1]){
  gene_profiles_formative[d,] <- gene_profiles_formative[d,]* 1.035
}

gene_profiles_formative_24 <- gene_profiles_formative[,1:18]
for(d in 1:dim(gene_profiles_formative)[1]){
  gene_profiles_formative_24[d,] <- gene_profiles_formative[d,1:18] - gene_profiles_formative[d,14]
}


for(d in 1:dim(gene_profiles_formative_24)[1]){
  row.names(gene_profiles_formative_24)[d] <- biomart[which(biomart$ensembl_gene_id==row.names(gene_profiles_formative_24)[d]),3]
}


KO_wt_mat_formative <- KOvWT[which(row.names(KOvWT)%in%row.names(gene_profiles_formative_24)),]
KO_wt_mat_formative <- as.matrix(KO_wt_mat_formative[match(row.names(gene_profiles_formative_24),row.names(KO_wt_mat_formative)),])

distances_KOs_formative <- distanceKOvsTP(KO_wt_mat_formative,gene_profiles_formative_24)

rel_distances_KOs_formative <- getRelativeDistance(distances_KOs_formative)


Heatmap(rel_distances_KOs_all,show_column_dend = FALSE, show_row_dend = FALSE,
        col = colorRamp2(c(1,0.1,0.1,0), c("white",'black','darkred','red')),
        show_row_names = T,show_column_names=TRUE,
        cluster_columns = FALSE,cluster_rows = TRUE, name="distance",
        gap=unit(5,"mm"),
        column_title = "timing of KOs\nbased on all DEGs\n")

Heatmap(rel_distances_KOs_naive,show_column_dend = FALSE, show_row_dend = FALSE,
        col = colorRamp2(c(1,0.1,0.1,0), c("white",'black','darkred','red')),
        show_row_names = T,show_column_names=TRUE,
        cluster_columns = FALSE,cluster_rows = TRUE, name="distance",
        gap=unit(5,"mm"),
        column_title = "timing of KOs\nbased on naive marker genes\n")

Heatmap(rel_distances_KOs_formative,show_column_dend = FALSE, show_row_dend = FALSE,
        col = colorRamp2(c(1,0.1,0.1,0), c("white",'black','darkred','red')),
        show_row_names = T,show_column_names=TRUE,
        cluster_columns = FALSE,cluster_rows = TRUE, name="distance",
        gap=unit(5,"mm"),
        column_title = "timing of KOs\nbased on formative marker genes\n")

## get best timings
KOs_rank <- as.data.frame(matrix(NA,nrow=73,ncol=12))
colnames(KOs_rank) <- c("timing naive marker","min euclidean distance naive marker",
                        "confidence interval naive marker","rank naive marker",
                        "timing formative marker","min euclidean distance formative marker",
                        "confidence interval formative marker","rank formative marker",
                        "timing DEGs","min euclidean distance DEGs",
                        "confidence interval DEGs","rank DEGs")

for(d in 1:dim(KOs_rank)[1]){
  row.names(KOs_rank)[d] <- row.names(distances_KOs_naive)[d]
  
  row <- which(distances_KOs_naive[d,]==min(distances_KOs_naive[d,]))
  row_name <- row.names(KOs_rank)[d]
  
  best_time_inter <- getGlobalTiming(gene_profiles_naive_24,distances_KOs_naive,
                                     row,KO_wt_mat_naive,row_name,
                                     duplicates=F,lookup=lookup)
  KOs_rank[d,1] <- best_time_inter[1]
  KOs_rank[d,2] <- best_time_inter[2]
  
  confidence_interval <- lookup[which(lookup$names%in%colnames(rel_distances_KOs_naive)[which(rel_distances_KOs_naive[d,]<0.1)]),2]
  if(length(unique(confidence_interval))==
     length(seq(from=min(confidence_interval), to=max(confidence_interval),by=2))){
    KOs_rank[d,3] <- paste(min(confidence_interval), "to", max(confidence_interval))
  }
  else{
    KOs_rank[d,3] <- paste(min(confidence_interval), "to", max(confidence_interval))
    print(paste(row.names(KOs_rank)[d],"has a gap in the confidence interval for naive markers"))
    print(rel_distances_KOs_naive[d,which(rel_distances_KOs_naive[d,]<0.1)])
  }
  
  formative_d <- which(row.names(distances_KOs_formative)==row.names(KOs_rank)[d])
  
  row <- which(distances_KOs_formative[formative_d,]==min(distances_KOs_formative[formative_d,]))
  
  best_time_inter <- getGlobalTiming(gene_profiles_formative_24,distances_KOs_formative,
                                     row,KO_wt_mat_formative,row_name,
                                     duplicates=F,lookup=lookup)
  KOs_rank[d,5] <- best_time_inter[1]
  KOs_rank[d,6] <- best_time_inter[2]
  confidence_interval <- lookup[which(lookup$names%in%colnames(rel_distances_KOs_formative)[which(rel_distances_KOs_formative[formative_d,]<0.1)]),2]
  if(length(unique(confidence_interval))==
     length(seq(from=min(confidence_interval), to=max(confidence_interval),by=2))){
    KOs_rank[d,7] <- paste(min(confidence_interval), "to", max(confidence_interval))
  }
  else{
    KOs_rank[d,7] <- paste(min(confidence_interval), "to", max(confidence_interval))
    print(paste(row.names(KOs_rank)[d],"has a gap in the confidence interval for formative markers"))
    print(rel_distances_KOs_formative[formative_d,which(rel_distances_KOs_formative[formative_d,]<0.1)])
  }
  
  all_d <- which(row.names(distances_KOs_all)==row.names(KOs_rank)[d])
  
  row <- which(distances_KOs_all[all_d,]==min(distances_KOs_all[all_d,]))
  
  best_time_inter <- getGlobalTiming(gene_profiles_WT_diff_24,distances_KOs_all,
                                     row,KO_wt_mat_all,row_name,
                                     duplicates=F,lookup=lookup)
  KOs_rank[d,9] <- best_time_inter[1]
  KOs_rank[d,10] <- best_time_inter[2]
  confidence_interval <- lookup[which(lookup$names%in%colnames(rel_distances_KOs_all)[which(rel_distances_KOs_all[all_d,]<0.1)]),2]
  if(length(unique(confidence_interval))==
     length(seq(from=min(confidence_interval), to=max(confidence_interval),by=2))){
    KOs_rank[d,11] <- paste(min(confidence_interval), "to", max(confidence_interval))
  }
  else{
    KOs_rank[d,11] <- paste(min(confidence_interval), "to", max(confidence_interval))
    print(paste(row.names(KOs_rank)[d],"has a gap in the confidence interval for DEGs"))
    print(paste(confidence_interval, ":",seq(from=min(confidence_interval), to=max(confidence_interval),by=2)))
  }
}

KOs_rank[,4] <- rank(KOs_rank$`timing naive marker`)
KOs_rank[,8] <- rank(KOs_rank$`timing formative marker`)
KOs_rank[,12] <- rank(KOs_rank$`timing DEGs`)
KOs_rank <- KOs_rank[order(KOs_rank$`rank naive marker`),]


rel_distances_KOs_all_ordered <- rel_distances_KOs_all[row.names(KOs_rank),]
rel_distances_KOs_naive_ordered <- rel_distances_KOs_naive[row.names(KOs_rank),]
rel_distances_KOs_formative_ordered <- rel_distances_KOs_formative[row.names(KOs_rank),]
#Fig 5.9 C
Heatmap(rel_distances_KOs_all_ordered,show_column_dend = FALSE, show_row_dend = FALSE,
        col = colorRamp2(c(1,0.1,0.1,0), c("white",'black','darkred','red')),
        show_row_names = T,show_column_names=TRUE,
        cluster_columns = FALSE,cluster_rows = F, name="distance",
        gap=unit(5,"mm"),
        column_title = "timing of KOs\nbased on all DEGs\n")
#Fig 5.9 A
Heatmap(rel_distances_KOs_naive_ordered,show_column_dend = FALSE, show_row_dend = FALSE,
        col = colorRamp2(c(1,0.1,0.1,0), c("white",'black','darkred','red')),
        show_row_names = T,show_column_names=TRUE,
        cluster_columns = FALSE,cluster_rows = F, name="distance",
        gap=unit(5,"mm"),
        column_title = "timing of KOs\nbased on naive marker genes\n")
#Fig 5.9 B
Heatmap(rel_distances_KOs_formative_ordered,show_column_dend = FALSE, show_row_dend = FALSE,
        col = colorRamp2(c(1,0.1,0.1,0), c("white",'black','darkred','red')),
        show_row_names = T,show_column_names=TRUE,
        cluster_columns = FALSE,cluster_rows = F, name="distance",
        gap=unit(5,"mm"),
        column_title = "timing of KOs\nbased on formative marker genes\n")


###
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart_2 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol"),
                 mart=ensembl)


lookup_2 <- biomart_2[which(biomart_2$mgi_symbol%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                      "Klf4", "Prdm14")),]

foldchanges_KOs <- CLIM2$singlecontrasts$fit_N1$coefficients[which(row.names(CLIM2$singlecontrasts$fit_N1$coefficients)%in%
                                                                     c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                       "Klf4", "Prdm14")),]

gpr_mean <- gpr_list[[1]][which(row.names(gpr_list[[1]])%in%lookup_2$ensembl_gene_id),2:18]
gpr_upper <- gpr_list[[2]][which(row.names(gpr_list[[2]])%in%lookup_2$ensembl_gene_id),]
gpr_lower <- gpr_list[[3]][which(row.names(gpr_list[[3]])%in%lookup_2$ensembl_gene_id),]


naive_timing <- data.frame(matrix(NA,nrow = 0,ncol=3))
colnames(naive_timing) <- c("time","log2FC","Gene")

for(d in 1:dim(gpr_mean)[1]){
  if(row.names(gpr_mean)[d]%in%lookup_2[which(lookup_2$mgi_symbol%in%
                                            c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                              "Klf4", "Prdm14")),1]){
    naive_timing <- rbind(naive_timing,data.frame(time=c(0,2,4,6,8,10,12,14,16,
                                                         18,20,22,24,26,28,30,32),
                                                  log2FC=gpr_mean[d,]-gpr_mean[d,13],
                                                  Gene=lookup_2[which(lookup_2$ensembl_gene_id==row.names(gpr_mean)[d]),3]))
  }
}

scatter_KOs <- data.frame(matrix(NA,nrow = 0,ncol=3))
colnames(scatter_KOs) <- c("time","log2FC","Gene")

scatter_KOs <- rbind(scatter_KOs,data.frame(time=17.5,
                                            log2FC=foldchanges_KOs[which(row.names(foldchanges_KOs)%in%
                                                                           c("Esrrb", "Nanog", "Tfcp2l1", "Tbx3",
                                                                             "Prdm14","Klf4"))
                                                                   ,"N1_Pten"],
                                            Gene=names(foldchanges_KOs[which(row.names(foldchanges_KOs)%in%
                                                                               c("Esrrb", "Nanog", "Tfcp2l1", "Tbx3",
                                                                                 "Prdm14","Klf4"))
                                                                       ,"N1_Pten"])))

scatter_KOs <- rbind(scatter_KOs,data.frame(time=27.5,
                                            log2FC=foldchanges_KOs[which(row.names(foldchanges_KOs)%in%
                                                                           c("Esrrb", "Nanog", "Tfcp2l1", "Tbx3",
                                                                             "Prdm14","Klf4"))
                                                                   ,"N1_Myc"],
                                            Gene=names(foldchanges_KOs[which(row.names(foldchanges_KOs)%in%
                                                                               c("Esrrb", "Nanog", "Tfcp2l1", "Tbx3",
                                                                                 "Prdm14","Klf4"))
                                                                       ,"N1_Myc"])))


#Fig 5.8
ggplot(naive_timing,aes(x=time,y=log2FC*1.454,col=Gene)) +
  geom_line(size=1) +
  #stat_summary(fun.y=mean,size=1,color="grey50",geom="line",linetype="dashed") +
  geom_point(data = scatter_KOs,aes(x=scatter_KOs$time,y=scatter_KOs$log2FC,
                                    col=scatter_KOs$Gene),
             size=3) + xlab("hours after 2i removal") + ylab("log2FC vs 24hours")+
  theme_bw()
###

sum(KOs_rank$`timing naive marker` >=22 & KOs_rank$`timing naive marker` <=26)
sum(KOs_rank$`timing formative marker` >=22 & KOs_rank$`timing formative marker` <=26)
sum(KOs_rank$`timing DEGs` >=22 & KOs_rank$`timing DEGs` <=26)
#Fig 5.B
ggplot(KOs_rank, aes(x=24-`timing naive marker`,y=24-`timing DEGs`)) +
  geom_point(shape=21, size=4) +
  geom_abline(slope = 1, intercept = 0)+
  xlab("Delay based on naive marker")+
  ylab("Delay based on DEGs")+
  xlim(c(-5,25))+
  ylim(c(-5,25))+
  theme_bw() + coord_fixed()
#Fig 5.A
ggplot(KOs_rank, aes(x=24-`timing naive marker`,y=24-`timing formative marker`)) +
  geom_point(shape=21, size=4) +
  geom_abline(slope = 1, intercept = 0)+
  xlab("Delay based on naive marker")+
  ylab("Delay based on formative marker")+
  xlim(c(-10,25))+
  ylim(c(-10,25))+
  theme_bw() + coord_fixed()

out_table <- KOs_rank[,c(1,1,5,5,9,9)]
out_table[,c(2,4,6)] <- 24-out_table[,c(2,4,6)]
colnames(out_table) <- c("timing naive marker", "delay naive marker",
                         "timing formative marker", "delay formative marker",
                         "timing DEGs", "delay DEGs")
row.names(out_table) <- sapply(row.names(out_table), function(x) strsplit(x,split = "_")[[1]][2])
write.csv(out_table, "output_tables/KO_ranks.csv",
          quote = F,row.names = T)

