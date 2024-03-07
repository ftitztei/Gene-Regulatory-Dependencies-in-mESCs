## This script aims to apply a similar concept to the single-cell data
## as it was applied to the KO data. Instead of completion percentage,
## the ratios of cells that show late expression of one gene and
## early expression of another gene & the number of cells that
## show early expression of one gene and late expression of the other.
## These possible dependencies are then tested for consistency with
## dependencies from the KO data.


library(ggplot2)
library(sva)
library(biomaRt)
library(GenomicFeatures)
library(gridExtra)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)

library(scran)
library(slingshot)
library(RColorBrewer)
library(psupertime)
library(SeuratObject)

library(org.Mm.eg.db)

## load functions for this script and other scripts in the analysis
source("R/functions.R")

### read in rds and objects from QC and prerocessing ####
SC_lookup <- readRDS("RDS/count_sct_Christa_lookup.rds")
SC_RC9 <- readRDS("RDS/SC_RC9_Christa.rds")
sc_data <- as.matrix(readRDS("RDS/count_sct_Christa.rds"))
colnames(sc_data) <- paste(SC_RC9@meta.data$type,
                                             1:length(SC_RC9@meta.data$type),sep="_")
sc_data_0to24 <- sc_data[,which(SC_RC9@meta.data$type%in%c("0h","6h","12h","24h"))]
gain_loss_ggplot <- readRDS("RDS/gain_loss_ggplot_sct_Christa.rds")

non_neg_cells <- which(SC_RC9@meta.data$type!="Negative")
T0_cells <- which(SC_RC9@meta.data$type=="0h")
T24_cells <- which(SC_RC9@meta.data$type=="24h")
T48_cells <- which(SC_RC9@meta.data$type=="48h")

time_points_0to24 <- sapply(colnames(sc_data_0to24),
                            function(x) strsplit(x,split="_")[[1]][1])
T0_cells <- which(time_points_0to24=="0h")
T24_cells <- which(time_points_0to24=="24h")

## biomart lookup between ensembl and mgi symbols
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")#, host = "uswest.ensembl.org")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol"),
                 mart=ensembl)
biomart <- biomart[which(!is.na(biomart$entrezgene)),]
biomart <- biomart[complete.cases(biomart),]


### pseudotime ####
sc_exp_in <- sc_data_0to24
row.names(sc_exp_in) <- SC_lookup$mgi
sc_exp_obj <- SingleCellExperiment(list(count=sc_exp_in))
logcounts(sc_exp_obj) <- log2(sc_exp_in)
pseudo_obj <- psupertime(x = sc_exp_obj,
                         y = factor(time_points_0to24,
                                    levels = c("0h", "6h", "12h", "24h")),
                         sel_genes = "all")#,
# gene_list = c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
#               "Klf4", "Prdm14","Fgf5", "Pou3f1", "Otx2",
#               "Dnmt3a", "Dnmt3b"))

plot_train_results(pseudo_obj)
plot_labels_over_psupertime(pseudo_obj, label_name='Time')
plot_identified_gene_coefficients(pseudo_obj)
plot_identified_genes_over_psupertime(pseudo_obj, label_name='Time')

pseudo_order <- order(pseudo_obj$proj_dt$psuper, decreasing = F)
sc_data_0to24_ordered <- sc_data_0to24[,pseudo_order]

### cell cycle prediction ####
## retrieve a lookup table for ensembl and mgi gene names
biomart_cc <- biomart[which(!is.na(biomart$entrezgene_id)),]
biomart_cc <- biomart_cc[complete.cases(biomart_cc),]

names_lookup <- data.frame(orig_name=rownames(sc_data_0to24),
                           mgi=NA,
                           ensembl=sapply(rownames(sc_data_0to24),function(x) strsplit(x,"_")[[1]][1]))
for(d in 1:dim(names_lookup)[1]){
  if(length(biomart_cc[which(biomart_cc$ensembl==names_lookup[d,3]),1])==1){
    names_lookup[d,2] <- biomart_cc[which(biomart_cc$ensembl==names_lookup[d,3]),3][1]
  }
}
## sort out esembl genes that were not unique
unique_ensembl <- names_lookup[which(!names_lookup$ensembl%in%
                                       unique(names_lookup[duplicated(names_lookup[,3]),3])),3]

## which genes are unique in the ensembl lookup and can be used in the cellcycle
## calculation
sc_data_0to24_cellcycle <- sc_data_0to24[which(row.names(sc_data_0to24)%in%unique_ensembl),]


mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
cell_cycle <- cyclone(as.matrix(sc_data_0to24_cellcycle),mm.pairs)
names(cell_cycle$phases) <- colnames(sc_data_0to24)

gg_cell_cycle <- data.frame(cellcycle=cell_cycle$phases,
                            pseudotime=pseudo_obj$proj_dt$psuper,
                            timepoint=time_points_0to24,
                            norm=NA)
gg_cell_cycle[which(gg_cell_cycle$timepoint=="0h"),4] <-
  gg_cell_cycle[which(gg_cell_cycle$timepoint=="0h"),2] -
  mean(gg_cell_cycle[which(gg_cell_cycle$timepoint=="0h"),2])
gg_cell_cycle[which(gg_cell_cycle$timepoint=="6h"),4] <-
  gg_cell_cycle[which(gg_cell_cycle$timepoint=="6h"),2] -
  mean(gg_cell_cycle[which(gg_cell_cycle$timepoint=="6h"),2])
gg_cell_cycle[which(gg_cell_cycle$timepoint=="12h"),4] <-
  gg_cell_cycle[which(gg_cell_cycle$timepoint=="12h"),2] -
  mean(gg_cell_cycle[which(gg_cell_cycle$timepoint=="12h"),2])
gg_cell_cycle[which(gg_cell_cycle$timepoint=="24h"),4] <-
  gg_cell_cycle[which(gg_cell_cycle$timepoint=="24h"),2] -
  mean(gg_cell_cycle[which(gg_cell_cycle$timepoint=="24h"),2])

gg_cell_cycle <- gg_cell_cycle[pseudo_order,]

ggplot(gg_cell_cycle, aes(x=pseudotime, col=cellcycle, fill=cellcycle))+
  geom_density(alpha=0.3)+
  facet_wrap(~factor(timepoint, levels=c("0h", "6h", "12h", "24h")), ncol=1)

ggplot(gg_cell_cycle, aes(x=norm, col=cellcycle, fill=cellcycle))+
  geom_density(alpha=0.3)+
  facet_wrap(~factor(timepoint, levels=c("0h", "6h", "12h", "24h")), ncol=1)

wilcox.test(gg_cell_cycle[which(gg_cell_cycle$cellcycle=="G1" &
                                  gg_cell_cycle$timepoint=="0h"),4],
            gg_cell_cycle[which(gg_cell_cycle$cellcycle=="G1" &
                                  gg_cell_cycle$timepoint=="24h"),4])

wilcox.test(gg_cell_cycle[which(gg_cell_cycle$cellcycle=="G1" &
                                  gg_cell_cycle$timepoint=="0h"),4],
            gg_cell_cycle[which(gg_cell_cycle$cellcycle=="G1" &
                                  gg_cell_cycle$timepoint=="6h"),4])

wilcox.test(gg_cell_cycle[which(gg_cell_cycle$cellcycle=="G1" &
                                  gg_cell_cycle$timepoint=="0h"),4],
            gg_cell_cycle[which(gg_cell_cycle$cellcycle=="G1" &
                                  gg_cell_cycle$timepoint=="12h"),4])

### identify genes that expressed at certain level in most cells ####
A_relation_b_hm_consist <- readRDS("RDS/A_relation_b_hm_consist_perc_genes.rds")
## set rownames and colnames to ensembl for easy handling
A_relation_b_hm_consist <- A_relation_b_hm_consist[row.names(A_relation_b_hm_consist)%in%SC_lookup$mgi,
                                                   colnames(A_relation_b_hm_consist)%in%SC_lookup$mgi]
for(d in 1:nrow(A_relation_b_hm_consist)){
  row.names(A_relation_b_hm_consist)[d] <- SC_lookup[which(SC_lookup$mgi==row.names(A_relation_b_hm_consist)[d]),1]
  colnames(A_relation_b_hm_consist)[d] <- row.names(A_relation_b_hm_consist)[d]
}

hist(unlist(as.matrix(log2(sc_data)),use.names = F))

## summarize how big the percentage of cells is that make a certain TPM
## cutoff either for all cells or for certain time points.
summary_thres <- data.frame(gene=row.names(sc_data),
## columns summarizing the percentage of cells making the cut off over all
## cells
                            bigger_1=apply(sc_data,1,function(x) (sum(x>1,na.rm=T)/
                                                                                                length(x))),
                            bigger_2=apply(sc_data,1,function(x) (sum(x>2,na.rm=T)/
                                                                                                length(x))),
                            bigger_3=apply(sc_data,1,function(x) (sum(x>3,na.rm=T)/
                                                                                                length(x))),
                            bigger_4=apply(sc_data,1,function(x) (sum(x>4,na.rm=T)/
                                                                                                length(x))),
                            bigger_5=apply(sc_data,1,function(x) (sum(x>5,na.rm=T)/
                                                                                                length(x))),
##  the percentage of cells having NAs over all cells
                            nas=apply(sc_data,1,function(x) (sum(x<0.1)/length(x))),#(sum(is.na(x))/length(x))),
## columns summarizing the percentage of cells making the cut off over all
## cells per time point
                            bigger_1_T0=apply(sc_data[,which(SC_RC9@meta.data$type=="0h")]
                                               ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_2_T0=apply(sc_data[,which(SC_RC9@meta.data$type=="0h")]
                                               ,1,function(x) (sum(x>2,na.rm=T)/length(x))),
                            bigger_1_T6=apply(sc_data[,which(SC_RC9@meta.data$type=="6h")]
                                               ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_2_T6=apply(sc_data[,which(SC_RC9@meta.data$type=="6h")]
                                               ,1,function(x) (sum(x>2,na.rm=T)/length(x))),
                            bigger_1_T12=apply(sc_data[,which(SC_RC9@meta.data$type=="12h")]
                                               ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_2_T12=apply(sc_data[,which(SC_RC9@meta.data$type=="12h")]
                                               ,1,function(x) (sum(x>2,na.rm=T)/length(x))),
                            bigger_1_T24=apply(sc_data[,which(SC_RC9@meta.data$type=="24h")]
                                              ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_2_T24=apply(sc_data[,which(SC_RC9@meta.data$type=="24h")]
                                              ,1,function(x) (sum(x>2,na.rm=T)/length(x))),
                            bigger_1_T48=apply(sc_data[,which(SC_RC9@meta.data$type=="48h")]
                                               ,1,function(x) (sum(x>1,na.rm=T)/length(x))),
                            bigger_2_T48=apply(sc_data[,which(SC_RC9@meta.data$type=="48h")]
                                               ,1,function(x) (sum(x>2,na.rm=T)/length(x))))
## prepare plot to compare cut offs
gg_thresh <- data.frame(percentage=c(summary_thres$bigger_1,
                                     summary_thres$bigger_2,
                                     summary_thres$bigger_3),
                        group=c(rep("count > 1",length(summary_thres$bigger_1)),
                                rep("count > 2",length(summary_thres$bigger_2)),
                                rep("count > 3",length(summary_thres$bigger_3))))

ggplot(gg_thresh, aes(x=percentage,color=group,fill=group)) +
  geom_density(alpha=0.2)+
  ggtitle("density plot of genes passing different expression thresholds") +
  theme(plot.title=element_text(hjust=0.5)) +
  xlim(0,0.05)

## read in time course data and define which genes change in the TC
gpr_list_shrunkenFCS <- readRDS("RDS/gpr_list_shrunkenFCS.rds")
log_norm_heatmap_gpr <- gpr_list_shrunkenFCS[[1]][,2:18]
log_norm_heatmap_gpr <- log_norm_heatmap_gpr-log_norm_heatmap_gpr[,1]
log_norm_heatmap_gpr_changed <- log_norm_heatmap_gpr[which(rowMax(log_norm_heatmap_gpr)-
                                                             rowMin(log_norm_heatmap_gpr)>=0.5),]
changed_genes_ensembl <- row.names(log_norm_heatmap_gpr_changed)
changed_genes_mgi <- biomart[which(biomart$ensembl_gene_id%in%changed_genes_ensembl),3]
## further summary what genes have a certain expression in at least
## 70 or 80 percent of the cells
genes_summary <- data.frame(row.names=row.names(sc_data),
                            gene=SC_lookup$mgi,
## columns summarizing whether a gene has at least 50 or 60
## percent of cells above a TPM cutoff
                            TPM_bigger_1_20_perc=summary_thres$bigger_1>=0.2,
                            TPM_bigger_1_30_perc=summary_thres$bigger_1>=0.3,
## columns summarizing whether a gene has less than 10
## percent of NAs per gene, mean, median and variance TPM per gene
                            Na_perc_lower_10=summary_thres$nas<0.1,
                            mean_TPM=apply(sc_data,1,function(x) mean(x,na.rm=T)),
                            median_TPM=apply(sc_data,1,function(x) median(x,na.rm=T)),
                            var_TPM=apply(sc_data,1,function(x) var(x,na.rm=T)),
## columns summarizing whether a gene has at least 70 or 80
## percent of cells above a TPM cutoff per time point
                            TPM_bigger_1_20_perc_T0=summary_thres$bigger_1_T0>=0.2,
                            TPM_bigger_1_30_perc_T0=summary_thres$bigger_1_T0>=0.3,
                            TPM_bigger_1_20_perc_T6=summary_thres$bigger_1_T6>=0.2,
                            TPM_bigger_1_30_perc_T6=summary_thres$bigger_1_T6>=0.3,
                            TPM_bigger_1_20_perc_T12=summary_thres$bigger_1_T12>=0.2,
                            TPM_bigger_1_30_perc_T12=summary_thres$bigger_1_T12>=0.3,
                            TPM_bigger_1_20_perc_T24=summary_thres$bigger_1_T24>=0.2,
                            TPM_bigger_1_30_perc_T24=summary_thres$bigger_1_T24>=0.3,
                            TPM_bigger_1_20_perc_T48=summary_thres$bigger_1_T48>=0.2,
                            TPM_bigger_1_30_perc_T48=summary_thres$bigger_1_T48>=0.3)


## for every gene check where highest and lowest expression is in TC (window)
## then give each gene a sign based on that + high expression late, - low expression late
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
## use windows (3 time points together) to be less vulnerable to a single outlier
for(d in 1:nrow(gpr_TPMs_short)){
  windows <- c(rowMeans(gpr_TPMs_short[d,1:3]),
               rowMeans(gpr_TPMs_short[d,2:4]),
               rowMeans(gpr_TPMs_short[d,3:5]),
               rowMeans(gpr_TPMs_short[d,4:6]),
               rowMeans(gpr_TPMs_short[d,5:7]),
               rowMeans(gpr_TPMs_short[d,6:8]),
               rowMeans(gpr_TPMs_short[d,7:9]),
               rowMeans(gpr_TPMs_short[d,8:10]),
               rowMeans(gpr_TPMs_short[d,9:11]),
               rowMeans(gpr_TPMs_short[d,10:12]),
               rowMeans(gpr_TPMs_short[d,11:13]),
               rowMeans(gpr_TPMs_short[d,12:14]),
               rowMeans(gpr_TPMs_short[d,13:15]),
               rowMeans(gpr_TPMs_short[d,14:16]),
               rowMeans(gpr_TPMs_short[d,15:17]))
## if all windows are min and max there is no change and we
## fill NA as min and max
  if(length(which(windows==min(windows)))==15 &&
     length(which(windows==max(windows)))==15){
    gpr_TPMs_short[d,19] <- NA
    gpr_TPMs_short[d,20] <- NA
## fill in min and max of gene window per gene
  } else{
    gpr_TPMs_short[d,19] <- which(windows==min(windows))[1]
    gpr_TPMs_short[d,20] <- which(windows==max(windows))[1]
## if window min is located before window max sign is positive
    if(gpr_TPMs_short[d,19] < gpr_TPMs_short[d,20]){
      gpr_TPMs_short[d,21] <- 1
## if window min is located after window max sign is negative
    } else if(gpr_TPMs_short[d,19] > gpr_TPMs_short[d,20]){
      gpr_TPMs_short[d,21] <- -1
## else no sign
    } else{
      gpr_TPMs_short[d,21] <- 0
    }

  }
}

## remove the directions that are inconsistent with SC data
## here we look at the changes between E5.5 and E3.5 as E4.5 has
## a very small number of cells
gain_loss_ggplot$change_24_vs_0 <- NA
for(d in 1:nrow(gain_loss_ggplot)){
## log2FCs E4.5 vs E3.5 and E5.5 vs E3.5 should not be na
  if(!is.na(gain_loss_ggplot[d,5])){
## if log2FC E5.5 vs E3.5 is positive sign is positive
    if(sign(gain_loss_ggplot[d,5])==1){
      gain_loss_ggplot[d,7] <- 1
    }gg_cell_cycle
## if log2FC E5.5 vs E3.5 is negative sign is negative
    if(sign(gain_loss_ggplot[d,5])==-1){
      gain_loss_ggplot[d,7] <- -1
    }
  }
}
genes_one_direct <- row.names(gain_loss_ggplot[which(gain_loss_ggplot$change_sign_all%in%c(-1,1)),])

## prepare a plot to summarize log2FCs in SC data
gg_SC_TC <- gain_loss_ggplot[,c(4,5)]
gg_SC_TC$TC_sign <- NA
for(d in 1:nrow(gg_SC_TC)){
  if(row.names(gg_SC_TC)[d]%in%row.names(gpr_TPMs_short)){
    if(length(gpr_TPMs_short[which(row.names(gpr_TPMs_short)==row.names(gg_SC_TC)[d]),21])==1){
      gg_SC_TC[d,3] <- gpr_TPMs_short[which(row.names(gpr_TPMs_short)==row.names(gg_SC_TC)[d]),21]
    }
  }
}

ggplot(gg_SC_TC, aes(x=log2FC_12_0,y=log2FC_24_0,
                     col=as.factor(TC_sign),fill=as.factor(TC_sign)))+
  geom_point()+
  scale_alpha_manual(values=c(0.2,0.6,1)) +
  theme_bw() +
  labs(x="log2FC 12 vs 0",
       y="log2FC 24 vs 0") +
  ggtitle("consisctency change in SC data") +
  theme(plot.title=element_text(hjust=0.5))


## are directions of change between time course and SC data the same
gpr_TPMs_short_consist <- gpr_TPMs_short
for(d in 1:nrow(gpr_TPMs_short_consist)){
  if(gpr_TPMs_short_consist[d,21]%in%c(1,-1)){
    if(row.names(gpr_TPMs_short_consist)[d] %in% row.names(gain_loss_ggplot)){
      if(sign(gain_loss_ggplot[row.names(gpr_TPMs_short_consist)[d],5])!=gpr_TPMs_short_consist[d,21]){
        gpr_TPMs_short_consist[d,21] <- 0
      }
    } else {
      gpr_TPMs_short_consist[d,21] <- 0
    }
  }
}

#changed_genes_signs <- changed_genes_ensembl[which(changed_genes_ensembl %in%
#                                                     row.names(gpr_TPMs_short_consist[which(gpr_TPMs_short_consist[,21]%in%c(-1,1)),]))]
genes_one_direct_consist <- genes_one_direct[which(genes_one_direct %in%
                             row.names(gpr_TPMs_short_consist[which(gpr_TPMs_short_consist[,21]%in%c(-1,1)),]))]

### get thresholds and min max from the rocs ####

gpr_TPMs_short_consist$thresh <- NA
gpr_TPMs_short_consist$roc_FPR <- NA
gpr_TPMs_short_consist$roc_TPR <- NA
gpr_TPMs_short_consist$max_roc_FPR <- NA
gpr_TPMs_short_consist$min_roc_TPR <- NA

for(d in 1:nrow(gpr_TPMs_short_consist)){
  if(d%%1000 == 0){
    print(d)
  }
  if(gpr_TPMs_short_consist[d,21]%in%c(1,-1)){ #&gpr_TPMs_short_consist[d,18]%in%genes_pass
    if(sum(is.na(sc_data_0to24[row.names(gpr_TPMs_short_consist)[d],
                                                     c(T0_cells,T24_cells)]))==0){
      holder <- ROC_gg(gene=as.character(row.names(gpr_TPMs_short_consist)[d]),
                       SC_frame=sc_data_0to24, #log2_CP10k_0to24
                       direction=gpr_TPMs_short_consist[d,21],
                       nbreaks = 30, mode="thresh",
                       early = T0_cells, late = T24_cells)
      gpr_TPMs_short_consist[d,22] <- holder[1]
      gpr_TPMs_short_consist[d,23] <- holder[2]
      gpr_TPMs_short_consist[d,24] <- holder[3]

      holder <- ROC_gg(gene=as.character(row.names(gpr_TPMs_short_consist)[d]),
                       SC_frame=sc_data_0to24,
                       direction=gpr_TPMs_short_consist[d,21],
                       nbreaks = 30, mode="thresh",
                       NA_cut = 0.1,
                       early = T0_cells, late = T24_cells, min_TPR=T)
      gpr_TPMs_short_consist[d,25] <- holder[2]
      gpr_TPMs_short_consist[d,26] <- holder[3]
    }

  }
}

gg_threshs <- gpr_TPMs_short_consist[,c(18,21:26)]
gg_threshs <- na.omit(gg_threshs)
dim(gg_threshs)
gg_threshs[which(gg_threshs$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                     "Klf4", "Prdm14","Fgf5", "Pou3f1", "Otx2",
                                     "Dnmt3a", "Dnmt3b")),]

dim(gg_threshs[which(gg_threshs$max_roc_FPR <= 0.3 &
                       gg_threshs$min_roc_TPR >= 0.7),])

dim(gg_threshs[which(gg_threshs$max_roc_FPR > 0.3 &
                       gg_threshs$min_roc_TPR < 0.7),])


#genes_thresh_short <- genes_thresh[which(genes_thresh$SC_thresh_TPs_seperate==1),]

SC_relations_left_thresh <- A_relation_b_hm_consist[which(row.names(A_relation_b_hm_consist)%in%genes_one_direct_consist),
                                                    which(colnames(A_relation_b_hm_consist)%in%genes_one_direct_consist)]


## preparing to plot possible dependencies left after sorting out
group <- rep(NA, nrow(SC_relations_left_thresh))
group[which(row.names(SC_relations_left_thresh)%in%
              SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                 "Klf4", "Prdm14")),1])]  <- "naive markers"
group[which(row.names(SC_relations_left_thresh)%in%
              SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                 "Dnmt3a", "Dnmt3b")),1])] <- "formative markers"

ha <- rowAnnotation(group=group,
                    col=list(group=c("naive markers"="#57BC77","formative markers"="#719AC4")))

clust_nmf <- readRDS("RDS/clusters_h.rds")
clust_nmf$mgi <- row.names(clust_nmf)
clust_nmf_short <- clust_nmf[which(clust_nmf$mgi%in%
                                     SC_lookup[which(SC_lookup$ensembl%in%row.names(SC_relations_left_thresh)),2]),]
clust_nmf_short$ensembl <- NA
for(d in 1:nrow(clust_nmf_short)){
  clust_nmf_short[d,3] <- SC_lookup[which(SC_lookup$mgi==clust_nmf_short[d,2]),1]
}

Heatmap(SC_relations_left_thresh,
        col=colorRamp2(c(-1.2,0,1.2), c("blue",'white',"red")),
        show_row_dend = F, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = F, # hclust(dist(SC_relations_left_thresh),method="ward.D2"),
        cluster_columns = F, #hclust(dist(SC_relations_left_thresh),method="ward.D2"),
        row_split = clust_nmf_short$cluster,
        column_split = clust_nmf_short$cluster,
        name = "mean activity diff A to B", left_annotation = ha, #right_annotation = ha,
        column_title = "dependencies between A (row) and B (column) \nfrom KO data, one direction in SC data\n")

SC_relations_left_thresh_marker <- SC_relations_left_thresh[which(row.names(SC_relations_left_thresh)%in%
                                 SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                    "Klf4", "Prdm14","Fgf5", "Pou3f1", "Otx2",
                                                                    "Dnmt3a", "Dnmt3b")),1]),
                         which(row.names(SC_relations_left_thresh)%in%
                                 SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                    "Klf4", "Prdm14","Fgf5", "Pou3f1", "Otx2",
                                                                    "Dnmt3a", "Dnmt3b")),1])]

row.names(SC_relations_left_thresh_marker) <- SC_lookup[which(SC_lookup$ensembl %in%
                  row.names(SC_relations_left_thresh_marker)),2][c(7,1:6)]
colnames(SC_relations_left_thresh_marker) <- row.names(SC_relations_left_thresh_marker)

Heatmap(SC_relations_left_thresh_marker,
        col=colorRamp2(c(-0.7,0,0.7), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = T, show_column_names = T,
        cluster_rows = hclust(dist(SC_relations_left_thresh_marker),method="ward.D2"),
        cluster_columns = hclust(dist(SC_relations_left_thresh_marker),method="ward.D2"),
        name = "mean activity diff A to B",
        column_title = "dependencies between A (row) and B (column) \nfrom KO data, one direction in SC data\n")


### use early late with threshold from SC ####
## create a holder for later results
## new contingency based on rank test ####
SC_dependencies <- get_relations_from_SC(SC_relations = SC_relations_left_thresh,
                                         SC_data = sc_data_0to24_ordered, 
                                         gpr_TPMs_short = gpr_TPMs_short_consist,
                                         early =  T0_cells,
                                         late = T24_cells,
                                         return_dependency = T,
                                         return_extreme_mean = F,
                                         extreme_cut = 0.2,
                                         contingency_out = F,
                                         early_late_cont = T,
                                         padj_cut = 0.05,
                                         counts = T,
                                         buffer_genes = 600)
### log2 of fraction on diagonal do positive and negative signs have effect on SC relations?2392

gg_contingency <- data.frame(Process_A=rep(NA,length(SC_dependencies[[2]])),
                             Process_B=rep(NA,length(SC_dependencies[[2]])),
                             dependency_SC=rep(NA,length(SC_dependencies[[2]])),
                             dependency_KO=rep(NA,length(SC_dependencies[[2]])),
                             padj_SC=rep(NA,length(SC_dependencies[[2]])),
                             #zero_relation=rep(NA,length(SC_dependencies[[2]])),
                             #zero_relation_pval=rep(NA,length(SC_dependencies[[2]])),
                             mgi_A=rep(NA,length(SC_dependencies[[2]])),
                             mgi_B=rep(NA,length(SC_dependencies[[2]])))

for(d in 1:nrow(gg_contingency)){
  gg_contingency[d,1] <- row.names(SC_dependencies[[1]])[which(SC_dependencies[[1]]==d, arr.ind = T)[1]]
  gg_contingency[d,2] <- row.names(SC_dependencies[[1]])[which(SC_dependencies[[1]]==d, arr.ind = T)[2]]

  gg_contingency[d,3] <- SC_dependencies[[2]][[d]]
  gg_contingency[d,4] <- SC_relations_left_thresh[which(SC_dependencies[[1]]==d, arr.ind = T)[1],
                                           which(SC_dependencies[[1]]==d, arr.ind = T)[2]]
  gg_contingency[d,5] <- SC_dependencies[[4]][[d]]
  #gg_contingency[d,6] <- log2((SC_dependencies[[6]][[d]][1,2]+1)/
  #                               (SC_dependencies[[6]][[d]][2,1]+1))
  #gg_contingency[d,7] <- fisher.test(SC_dependencies[[6]][[d]])$p.value
  gg_contingency[d,6] <- SC_lookup[which(SC_lookup$ensembl==gg_contingency[d,1]),2]
  gg_contingency[d,7] <- SC_lookup[which(SC_lookup$ensembl==gg_contingency[d,2]),2]



}

ggplot(gg_contingency)+
  geom_point(aes(x=dependency_KO,y=dependency_SC,
                 col=-log10(padj_SC),
                            fill=-log10(padj_SC)
                 ),alpha=0.3,
             size=2) +
  # scale_fill_gradient2(low = "bl|ue", mid = "white",
  #                       high = "red")+
  # scale_colour_gradient2(low = "blue", mid = "white",
  #                         high = "red")+
  theme_bw() +
  labs(x="KO based dependency A to B",
       y="sc based dependency A to B") +
  ggtitle("sc consistency matrices vs. process relations (using thresholds from sc data)") +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_vline(xintercept = 0, colour="black") +
  geom_hline(yintercept = 0, colour="black")


## Fig 5.27
ggplot(gg_contingency)+
  geom_point(aes(x=dependency_KO,y=dependency_SC),color="grey40",fill="grey40",alpha=0.3,
             size=2) +
  geom_point(data=subset(gg_contingency, Process_A %in% SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                           "Klf4", "Prdm14")),1]),
             aes(x=dependency_KO, y=dependency_SC),color="#57BC77",fill="#57BC77",alpha=0.5,
             size=2) +
  geom_point(data=subset(gg_contingency, Process_B %in% SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                           "Klf4", "Prdm14")),1]),
             aes(x=dependency_KO, y=dependency_SC),color="#57BC77",fill="#57BC77",alpha=0.5,
             size=2) +
  geom_point(data=subset(gg_contingency, Process_A %in% SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                           "Dnmt3a", "Dnmt3b")),1]),
             aes(x=dependency_KO, y=dependency_SC),color="#719AC4",fill="#719AC4",alpha=0.5,
             size=2) +
  geom_point(data=subset(gg_contingency, Process_B %in% SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                           "Dnmt3a", "Dnmt3b")),1]),
             aes(x=dependency_KO, y=dependency_SC),color="#719AC4",fill="#719AC4",alpha=0.5,
             size=2) +
  theme_bw() +
  labs(x="KO based dependency A to B",
       y="sc based deendency A to B") +
  ggtitle("sc consistency matrices vs. process relations (using thresholds from sc data)") +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_vline(xintercept = 0, colour="black") +
  geom_hline(yintercept = 0, colour="black") +

  annotate(geom="text",x=0.9, y=0.79,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC>0&
                                             gg_contingency$dependency_KO>0),]),
           colour="grey40") +
  annotate(geom="text",x=0.9, y=0.75,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC>0&
                                             gg_contingency$dependency_KO>0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                             "Klf4", "Prdm14")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                               "Klf4", "Prdm14")),1])),]),
           colour="#57BC77") +
  annotate(geom="text",x=0.9, y=0.71,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC>0&
                                             gg_contingency$dependency_KO>0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                             "Dnmt3a", "Dnmt3b")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                               "Dnmt3a", "Dnmt3b")),1])),]),
           colour="#719AC4") +

  annotate(geom="text",x=-0.9, y=0.79,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC>0&
                                             gg_contingency$dependency_KO<0),]),
           colour="grey40") +
  annotate(geom="text",x=-0.9, y=0.75,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC>0&
                                             gg_contingency$dependency_KO<0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                             "Klf4", "Prdm14")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                               "Klf4", "Prdm14")),1])),]),
           colour="#57BC77") +
  annotate(geom="text",x=-0.9, y=0.71,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC>0&
                                             gg_contingency$dependency_KO<0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                             "Dnmt3a", "Dnmt3b")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                               "Dnmt3a", "Dnmt3b")),1])),]),
           colour="#719AC4") +

  annotate(geom="text",x=-0.9, y=-0.71,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC<0&
                                             gg_contingency$dependency_KO<0),]),
           colour="grey40") +
  annotate(geom="text",x=-0.9, y=-0.75,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC<0&
                                             gg_contingency$dependency_KO<0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                             "Klf4", "Prdm14")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                               "Klf4", "Prdm14")),1])),]),
           colour="#57BC77") +
  annotate(geom="text",x=-0.9, y=-0.79,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC<0&
                                             gg_contingency$dependency_KO<0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                             "Dnmt3a", "Dnmt3b")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                               "Dnmt3a", "Dnmt3b")),1])),]),
           colour="#719AC4") +

  annotate(geom="text",x=0.9, y=-0.71,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC<0&
                                             gg_contingency$dependency_KO>0),]),
           colour="grey40") +
  annotate(geom="text",x=0.9, y=-0.75,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC<0&
                                             gg_contingency$dependency_KO>0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                             "Klf4", "Prdm14")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                                                               "Klf4", "Prdm14")),1])),]),
           colour="#57BC77") +
  annotate(geom="text",x=0.9, y=-0.79,
           label=nrow(gg_contingency[which(gg_contingency$dependency_SC<0&
                                             gg_contingency$dependency_KO>0&
                                             (gg_contingency$Process_A%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                             "Dnmt3a", "Dnmt3b")),1]|
                                                gg_contingency$Process_B%in%SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                                                                               "Dnmt3a", "Dnmt3b")),1])),]),
           colour="#719AC4")


gg_contingency$group <- NA
gg_contingency[which(gg_contingency$dependency_KO > 0),8] <- "positive"
gg_contingency[which(gg_contingency$dependency_KO < 0),8] <- "negative"
gg_contingency$group <- as.factor(gg_contingency$group)

gg_contingency$seperatability <- NA
gg_contingency$minCP10K <- NA
for(d in 1:nrow(gg_contingency)){
  gg_contingency$seperatability[d] <- min(gg_threshs[as.character(gg_contingency[d,1:2]),7])
  gg_contingency$minCP10K[d] <- min(genes_summary[as.character(gg_contingency[d,1:2]),5])
}


## safe all process process relations depending on which quadrants they are located in
Q1 <- gg_contingency[which(gg_contingency$dependency_SC > 0 & gg_contingency$dependency_KO > 0),]
Q2 <- gg_contingency[which(gg_contingency$dependency_SC > 0 & gg_contingency$dependency_KO < 0),]
Q3 <- gg_contingency[which(gg_contingency$dependency_SC < 0 & gg_contingency$dependency_KO < 0),]
Q4 <- gg_contingency[which(gg_contingency$dependency_SC < 0 & gg_contingency$dependency_KO > 0),]

(nrow(Q1)+nrow(Q3))/(nrow(Q2)+nrow(Q4))
nrow(Q1)+nrow(Q3)
### QC plots for process process relations ####

gg_contingency$quadrant <- NA

gg_contingency[which(gg_contingency$dependency_SC > 0 & gg_contingency$dependency_KO > 0),11] <- 1
gg_contingency[which(gg_contingency$dependency_SC > 0 & gg_contingency$dependency_KO < 0),11] <- 2
gg_contingency[which(gg_contingency$dependency_SC < 0 & gg_contingency$dependency_KO < 0),11] <- 3
gg_contingency[which(gg_contingency$dependency_SC < 0 & gg_contingency$dependency_KO > 0),11] <- 4


## sort out possible process process dependencies by quadrants
## only pairs in the first and third quadrant are consistent between KO data and SC
gg_contingency_Q1_Q3 <- gg_contingency[which(gg_contingency$quadrant%in%c(1,3)),]


SC_relations_consist_Q13 <- SC_relations_left_thresh
SC_relations_consist_Q13[1:383,1:383] <- 0

for(d in 1:nrow(gg_contingency_Q1_Q3)){
  SC_relations_consist_Q13[gg_contingency_Q1_Q3[d,1],gg_contingency_Q1_Q3[d,2]] <- 1
  SC_relations_consist_Q13[gg_contingency_Q1_Q3[d,2],gg_contingency_Q1_Q3[d,1]] <- 1
}

SC_relations_Q1_Q3 <- SC_relations_left_thresh * SC_relations_consist_Q13
#process_thresh[which(process_thresh$Processes%in%processes_left),2]
genes_left <- row.names(SC_relations_Q1_Q3)[which(rowSums(SC_relations_Q1_Q3)!=0)]

SC_relations_Q1_Q3 <- SC_relations_Q1_Q3[genes_left,genes_left]

group <- rep(NA, nrow(SC_relations_Q1_Q3))
group[which(row.names(SC_relations_Q1_Q3)%in%
              SC_lookup[which(SC_lookup$mgi%in%c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                 "Klf4", "Prdm14")),1])]  <- "naive markers"
group[which(row.names(SC_relations_Q1_Q3)%in%
              SC_lookup[which(SC_lookup$mgi%in%c("Fgf5", "Pou3f1", "Otx2",
                                                 "Dnmt3a", "Dnmt3b")),1])]  <- "formative markers"

ha <- rowAnnotation(group=group,
                    col=list(group=c("naive markers"="#57BC77","formative markers"="#719AC4")))

Heatmap(SC_relations_left_thresh[genes_left, genes_left],
        col=colorRamp2(c(-1.2,0,1.2), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = hclust(dist(SC_relations_Q1_Q3),method="ward.D2"),
        cluster_columns = hclust(dist(SC_relations_Q1_Q3),method="ward.D2"),
        name = "red: row independent of column\nblue: column independent of row",
        left_annotation = ha, #right_annotation = ha,
        column_title = "dependencies between A (row) and B (column) \nnot filtered\n")

Heatmap(SC_relations_Q1_Q3,
        col=colorRamp2(c(-1.2,0,1.2), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = hclust(dist(SC_relations_Q1_Q3),method="ward.D2"),
        cluster_columns = hclust(dist(SC_relations_Q1_Q3),method="ward.D2"),
        name = "red: row independent of column\nblue: column independent of row",
        left_annotation = ha, #right_annotation = ha,
        column_title = "dependencies between A (row) and B (column) \nconsistent between KO and SC data\n")
## define dependencies where contradiction is removed but rest is kept ####
A_relation_b_hm_nocontra <- readRDS("RDS/A_relation_b_hm_consist_perc_genes.rds")

for(d in 1:nrow(Q2)){
  A_relation_b_hm_nocontra[Q2[d,6],Q2[d,7]] <- 0
  A_relation_b_hm_nocontra[Q2[d,7],Q2[d,6]] <- 0
}
for(d in 1:nrow(Q4)){
  A_relation_b_hm_nocontra[Q4[d,6],Q2[d,7]] <- 0
  A_relation_b_hm_nocontra[Q4[d,7],Q2[d,6]] <- 0
}


saveRDS(SC_relations_Q1_Q3, "RDS/A_relation_b_hm_consist_KO_SC_genes_Christa_sct_counts_600.rds")
saveRDS(genes_one_direct_consist, "RDS/genes_one_direct_consist.rds")
saveRDS(A_relation_b_hm_consist, "RDS/A_relation_b_hm_consist_ensembl.rds")
saveRDS(gg_cell_cycle, "RDS/gg_cell_cycle.rds")
saveRDS(A_relation_b_hm_nocontra, "RDS/A_relation_b_hm_nocontra.rds")


### ####
### plot cells outside buffer results ####
library(scales)
gg_buff <- data.frame(buff=rep(c(100, 250,
                                 400, 450,
                                 500, 550,
                                 600, 650,
                                 700, 750,
                                 1000, 1250,
                                 1500, 2000),2),
                      res=c(1.51, 1.67,
                            1.81, 1.87,
                            1.94, 2.02,
                            2.11, 2.20,
                            2.26, 2.35,
                            2.65, 2.76,
                            2.41, 1.32,

                            rescale(c(11092, 9947, 8325, 7887,
                                      7469, 7155, 6893, 6642,
                                      6281, 6011, 4483, 2911,
                                      1199, 308), to=c(1.32,2.76))),
                      type=rep(c("ratio", "total"),each = 14))
## Fig 5.26
ggplot(gg_buff, aes(x=buff, y=res, color=type))+
  geom_vline(xintercept = 600, alpha = 0.3, lty = 2)+
  geom_line(linewidth = 1.2)+
  geom_point(size=2)+
  labs(x = "cut-off (cells outside buffer zone)",
       y = "ratio (agreement/disaggrement)",
       color ="read out")+
  scale_color_manual(values = c(
    "ratio" ="#FB8072",
    "total" = "#80B1D3")) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5),
        legend.position = "right")+
  scale_y_continuous(sec.axis = dup_axis(name ="# of dependencies in agreement",
                                         breaks=rescale(c(11092, 308, 500,
                                                          2500, 5000, 7500, 10000),
                                                        to=c(1.32,2.76))[3:7],
                                         labels=c("500", "2500", "5000", "7500", "10000")))

## test gene pairs ###

## Pou3f1 ENSMUSG00000090125
## Klf4 ENSMUSG00000003032
## Dnmt3a ENSMUSG00000020661
## Nanog ENSMUSG00000012396

## Fig 5.25 A & B
test <- gene_gene_dependencies(Gene_A = "ENSMUSG00000090125",
                               Gene_B = "ENSMUSG00000003032",
                               SC_data = sc_data_0to24_ordered,#noise_dat,#test_mat,#
                               gpr_TPMs_short=gpr_TPMs_short_consist,
                               #NA_cut = -0.5,
                               early = T0_cells, late = T24_cells,
                               return_dependency = F,
                               return_extreme_mean=T,
                               extreme_cut=0.2,
                               counts = T)

gg_test <- data.frame(Gene_A=test[[1]],#noise_dat["ENSMUSG00000038648",],#
                      Gene_B=test[[2]],#noise_dat["ENSMUSG00000020661",],#
                      time=test[[3]],
                      nonmt.libsize=SC_RC9[["nonmt.libsize"]][as.numeric(sapply(names(test[[3]]),
                                                                                function(x) strsplit(x,"_")[[1]][2])),1],
                      number.genes=SC_RC9@meta.data$nFeature_RNA[as.numeric(sapply(names(test[[3]]),
                                                                                   function(x) strsplit(x,"_")[[1]][2]))],
                      avg.CP10K=test[[4]],
                      counts_A=SC_RC9@assays[["RNA"]]@counts["ENSMUSG00000090125",
                                                             as.numeric(sapply(names(test[[3]]),
                                                                               function(x) strsplit(x,"_")[[1]][2]))],
                      counts_B=SC_RC9@assays[["RNA"]]@counts["ENSMUSG00000003032",
                                                             as.numeric(sapply(names(test[[3]]),
                                                                               function(x) strsplit(x,"_")[[1]][2]))])
gg_test$time <- factor(gg_test$time, levels=c("0h","6h","12h","24h","48h","Negative"))

ggplot(gg_test,
       aes(x=rank(Gene_A,ties.method = "first"), #
           y=rank(Gene_B,ties.method = "first"), alpha=0.7,#
           col=time, fill=time))+
  geom_point() +
  theme_bw() +
  labs(x=paste("rank of ","Pou3f1",sep = ""),
       y=paste("rank of ","Klf4",sep = "")) +
  ggtitle("") +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_abline(slope=1,intercept=0) +
  coord_fixed()

rank_diff <- rank(gg_test$Gene_A,ties.method = "first")-
                rank(gg_test$Gene_B,ties.method = "first")

buffer <- max(max(abs(rank_diff[which(gg_test$counts_A==0 &
                                        gg_test$counts_B==0 &
                                        gg_test$time == "24h")])),
              max(abs(rank_diff[which(gg_test$counts_A==0 &
                                        gg_test$counts_B==0 &
                                        gg_test$time == "24h")])))

buffer <- buffer * 1.1

ggplot(gg_test,
       aes(x=rank(Gene_A,ties.method = "first"), #
           y=rank(Gene_B,ties.method = "first"), alpha=0.7,#
           col=time, fill=time))+
  geom_point() +
  theme_bw() +scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                                          "darkorange3","coral3",
                                          "darkred")) +
  scale_fill_manual(values=c("darkgoldenrod1","darkorange1",
                             "darkorange3","coral3",
                             "darkred")) +
  labs(x=paste("rank of ","Dnmt3a",sep = ""),
       y=paste("rank of ","Nanog",sep = "")) +
  ggtitle("") +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_abline(slope=1,intercept=0) +
  geom_abline(slope=1,intercept=buffer) +
  geom_abline(slope=1,intercept=-buffer) +
  coord_fixed()


