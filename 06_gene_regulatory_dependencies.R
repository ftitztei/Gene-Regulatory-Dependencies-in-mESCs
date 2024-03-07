## Here we aim to find possible gene regulatory dependencies based on
## behaviour of connected genes across the 73 KOs.
## The main idea is to look at how much of the change that happens in the WT
## to a gene already happened in a KO. This results in a completion percentage
## per gene per KO and one can see if there is a consistent relation between
## two genes (i.e. in gene A genes show higher completion percentage across KOs
## than in gene B thus B might depend on A)

library(ggplot2)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)

### functions and data ####
## load functions for this script and other scripts in the analysis
source("R/functions.R")
## loading objects from previous scripts and analysis
CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
gpr_list <- readRDS("RDS/gpr_list_shrunkenFCS.rds")
enrichment_set <- readRDS("RDS/enrichment_set.rds")
ggplot_hist_repr_induct_complete <- readRDS("RDS/ggplot_hist_repr_induct_complete.rds")

### prepare objects needed ####
## biomart lookup between ensembl and mgi symbols
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol"),
                 mart=ensembl)
biomart <- biomart[which(!is.na(biomart$mgi_symbol)),]
#biomart <- biomart[complete.cases(biomart),]

## revtrieve log2FCs and adjusted pvalues for 24 hours vs 0 hours (2i)
## for all KOs and WT
percentage_table <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
KOvKO_padj <- CLIM2$singlecontrasts$fit_KOVKO$adj.P.Value

#saveRDS(FPKM, "RDS/FPKM.rds")
#FPKM <- read.csv("/data/public/ftitztei/100KO/2i_FPKM.csv")
#FPKM <- FPKM[which(FPKM$symbol%in%enrichment_set),
#             c(1,2,25,26,49,50,85,86,115,116,145,146,159,160)]
FPKM <- readRDS("RDS/FPKM.rds")

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

## determine which genes are regulated how often
summary_reg <- as.data.frame(table(ggplot_hist_repr_induct_complete$gene))

mono_reg_ensembl <- as.character(summary_reg[which(summary_reg$Freq==1),1])
dual_reg_ensembl <- as.character(summary_reg[which(summary_reg$Freq==2),1])
tri_reg_ensembl <- as.character(summary_reg[which(summary_reg$Freq==3),1])

mono_reg_mgi <- unique(biomart[which(biomart$ensembl_gene_id%in%mono_reg_ensembl),3])
dual_reg_mgi <- unique(biomart[which(biomart$ensembl_gene_id%in%dual_reg_ensembl),3])
tri_reg_mgi <- unique(biomart[which(biomart$ensembl_gene_id%in%tri_reg_ensembl),3])

mono_reg_diff_ensembl <- mono_reg_ensembl[which(abs(log_norm_heatmap_gpr[mono_reg_ensembl,14])>=0.5)]
dual_reg_diff_ensembl <- dual_reg_ensembl[which(abs(log_norm_heatmap_gpr[dual_reg_ensembl,14])>=0.5)]
tri_reg_diff_ensembl <- tri_reg_ensembl[which(abs(log_norm_heatmap_gpr[tri_reg_ensembl,14])>=0.5)]

mono_reg_diff_mgi <-  unique(biomart[which(biomart$ensembl_gene_id%in%mono_reg_diff_ensembl),3])
dual_reg_diff_mgi <-  unique(biomart[which(biomart$ensembl_gene_id%in%dual_reg_diff_ensembl),3])
tri_reg_diff_mgi <-  unique(biomart[which(biomart$ensembl_gene_id%in%tri_reg_diff_ensembl),3])
## calculate how much of the WT change (24 hours vs 2i) has already happened
## in the KOs between 24 hours and 2i
percentage_list <- list()

for(g in 1:length(enrichment_set)){
  percentage_list <- c(percentage_list, list(KOvKO[enrichment_set[g],1:73]/
                                               KOvKO[enrichment_set[g],74]))
}

names(percentage_list) <- enrichment_set
## average across genes per KO for each of the processes
percentage_table <- as.data.frame(matrix(NA,nrow=sum(!is.na(percentage_list)),ncol=73))
row.names(percentage_table) <- names(percentage_list[!is.na(percentage_list)])
colnames(percentage_table) <- sapply(colnames(KOvKO), function(x) strsplit(x,split = "_")[[1]][2])[1:73]

for(p in 1:nrow(percentage_table)){
  percentage_table[p,] <- percentage_list[[row.names(percentage_table)[p]]]
}

percentage_table <- as.matrix(percentage_table)

gg_bias <- data.frame(gene = character(),
                      median_wt_compl = numeric(),
                      mean_wt_compl = numeric(),
                      direction_wt = numeric(),
                      FPKM_2i = numeric())

for(gene in enrichment_set){
  if(gene %in% FPKM$symbol){
    gg_bias <- rbind(gg_bias,
                     data.frame(gene = gene,
                                median_wt_compl = median(as.numeric(percentage_table[gene,])),
                                mean_wt_compl = mean(as.numeric(percentage_table[gene,])),
                                direction_wt = KOvKO[gene,74],
                                FPKM_2i = mean(as.numeric(FPKM[which(FPKM$symbol==gene),
                                                    3:14]))))
  }
}



ggplot(gg_bias, aes(x=direction_wt, y=median_wt_compl)) +
  geom_point() +
  xlab("WT change N24 vs 2i")+
  ylab("median % of WT completion")+
  ggtitle("WT completion vs WT change") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ggplot(gg_bias) +
  geom_density(data=subset(gg_bias, direction_wt < 0),
               aes(x=median_wt_compl), fill="blue",
               colour="blue", alpha=0.3) +
  geom_density(data=subset(gg_bias, direction_wt > 0),
               aes(x=median_wt_compl), fill="red",
               colour="red", alpha=0.3) +
  xlab("median % of WT completion")+
  ggtitle("Density median % of WT completion \npositive (red) and negative (blue) WT change\n") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ks.test(gg_bias[which(gg_bias$direction_wt>0),2],
        gg_bias[which(gg_bias$direction_wt<0),2])


ggplot(gg_bias, aes(x=log10(FPKM_2i+1), y=median_wt_compl)) +
  geom_point() +
  xlab("log10(FPKM+1) in 2i")+
  ylab("median % of WT completion")+
  ggtitle("WT completion vs expression level in 2i") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))
##
ggplot(gg_bias, aes(x=direction_wt, y=mean_wt_compl)) +
  geom_point() +
  xlab("WT change N24 vs 2i")+
  ylab("mean % of WT completion")+
  ggtitle("WT completion vs WT change") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

ggplot(gg_bias) +
  geom_density(data=subset(gg_bias, direction_wt < 0),
               aes(x=mean_wt_compl), fill="blue",
               colour="blue", alpha=0.3) +
  geom_density(data=subset(gg_bias, direction_wt > 0),
               aes(x=mean_wt_compl), fill="red",
               colour="red", alpha=0.3) +
  xlab("mean % of WT completion")+
  ggtitle("Density mean % of WT completion \npositive (red) and negative (blue) WT change\n") +
  theme(plot.title=element_text(hjust=0.5),
        ggplot(gg_bias) +
          geom_density(data=subset(gg_bias, direction_wt < 0),
                       aes(x=median_wt_compl), fill="blue",
                       colour="blue", alpha=0.3) +
          geom_density(data=subset(gg_bias, direction_wt > 0),
                       aes(x=median_wt_compl), fill="red",
                       colour="red", alpha=0.3) +
          xlab("median % of WT completion")+
          ggtitle("Density median % of WT completion \npositive (red) and negative (blue) WT change\n") +
          theme(plot.title=element_text(hjust=0.5),
                panel.background = element_rect(fill=NA),
                panel.grid.major = element_line(colour = "grey90")),
        panel.grid.major = element_line(colour = "grey90"))

ks.test(gg_bias[which(gg_bias$direction_wt>0),3],
        gg_bias[which(gg_bias$direction_wt<0),3])


ggplot(gg_bias, aes(x=log10(FPKM_2i+1), y=mean_wt_compl)) +
  geom_point() +
  xlab("log10(FPKM+1) in 2i")+
  ylab("mean % of WT completion")+
  ggtitle("WT completion vs expression level in 2i") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

## calculate relations from completion percentage with pre defined functions
relation_list <- get_relations_from_percentage_genewise(percentage_table=percentage_table)
## each list element corresponds to other to one of the following
## relation A to B (sign). Here, relations not significantly different
## from zero were already excluded
A_relation_B_matrix <- relation_list[[1]]
## mean relation between processes
A_mean_B_matrix <- relation_list[[2]]
## adjusted pval is the distribution A-B signigicantly different from 0
A_padj_B_matrix <- relation_list[[3]]
## pval is the distribution A minus B signigicantly different from 0
A_pval_B_matrix <- relation_list[[4]]
## how consisten is the directionality between A and B with the mean
consistency_frame <- relation_list[[5]]

## combine directionality with absolute mean relation
A_relation_b_hm <- A_relation_B_matrix * abs(A_mean_B_matrix)

## which relations are non 0
included_relations <- rep(NA, nrow(A_relation_b_hm)*ncol(A_relation_b_hm))
for(r in 1:dim(A_relation_b_hm)[1]){
  for(c in 1:dim(A_relation_b_hm)[2]){
    if(A_relation_b_hm[r,c]!=0){
      included_relations[as.numeric(((r-1)*nrow(A_relation_b_hm)+c))] <-
        paste(row.names(A_relation_b_hm)[r],
              colnames(A_relation_b_hm)[c], sep=";")
    }
  }
}
included_relations <- na.omit(included_relations)
## extend consistency frame for relation between processes across KOs
consistency_frame$merged <- paste(consistency_frame$A,consistency_frame$B,sep=";")
consistency_frame$direction <- consistency_frame$relation
consistency_frame[which(consistency_frame$mean<0),6] <- -consistency_frame[which(consistency_frame$mean<0),6]

## only summarize for non 0 relationships
summary_consistency <- relation_list[[6]]
row.names(summary_consistency) <- summary_consistency$pair
summary_consistency <- as.matrix(summary_consistency[which(summary_consistency$pair%in%included_relations),1:2])
## plot consistency of relation between process across KOs
ggplot(as.data.frame(summary_consistency),aes(x=fraction_same_side)) +
  geom_histogram() +
  xlab("fraction of KOs showing same direction as mean relation")+
  ggtitle("Consistency of direction in relations") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

## only keep dependencies where more than 70 percent
## of the KOs show the same direction
consist_mat <- data.frame(matrix(0, nrow = dim(A_relation_b_hm)[1],
                                 ncol = dim(A_relation_b_hm)[2]))
colnames(consist_mat) <- colnames(A_relation_b_hm)
row.names(consist_mat) <- row.names(A_relation_b_hm)

for(r in 1:dim(summary_consistency)[1]){
  if(r%%50000==0){
    print(r)
  }
  if(summary_consistency[r,1] >= 0.7){ ## cut off can be adjusted
    consist_mat[strsplit(row.names(summary_consistency)[r],";")[[1]][1],
                strsplit(row.names(summary_consistency)[r],";")[[1]][2]] <- 1
    consist_mat[strsplit(row.names(summary_consistency)[r],";")[[1]][2],
                strsplit(row.names(summary_consistency)[r],";")[[1]][1]] <- 1
  }
}

## trim down relations to those consistent in more than 70 percent of the KOs
A_relation_b_hm_consist <- consist_mat * A_relation_b_hm
group <- rep(NA, nrow(A_relation_b_hm_consist))
group[which(row.names(A_relation_b_hm_consist)%in%
              c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                "Klf4", "Prdm14"))]  <- "naive markers"
group[which(row.names(A_relation_b_hm_consist)%in%
              c("Fgf5", "Pou3f1", "Otx2",
                "Dnmt3a", "Dnmt3b"))]  <- "formative markers"

ha <- rowAnnotation(group=group,
                     col=list(group=c("naive markers"="#57BC77","formative markers"="#719AC4")))

Heatmap(A_relation_b_hm_consist,
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = hclust(dist(A_relation_b_hm),method="ward.D2"),
        cluster_columns = hclust(dist(A_relation_b_hm),method="ward.D2"),
        name = "mean activity diff A to B", left_annotation = ha, #right_annotation = ha,
        column_title = "dependencies between A (row) and B (column) \nfrom KO data\n")

A_relation_b_hm_consist_safe <- A_relation_b_hm_consist[which(row.names(A_relation_b_hm_consist)%in%
                                                                c(mono_reg_mgi,dual_reg_diff_mgi)),
                                                        which(colnames(A_relation_b_hm_consist)%in%
                                                                c(mono_reg_mgi,dual_reg_diff_mgi))]

sum(apply(A_relation_b_hm_consist, 1, function(x) x != 0))
nrow(A_relation_b_hm_consist)*ncol(A_relation_b_hm_consist)

### test if certain KOs are most strong or contradicting ####
indexes_mat <- which(A_relation_b_hm_consist!=0, arr.ind = TRUE)
# remove lower triangle as it is the same just opposing sign
indexes_mat <- indexes_mat[which(indexes_mat[,1] < indexes_mat[,2]),]

gene_pairs_nonzero <- rep("", nrow(indexes_mat))
for(d in 1:length(gene_pairs_nonzero)){
  gene_pairs_nonzero[d] <- paste(row.names(A_relation_b_hm_consist)[indexes_mat[d,1]],
                                 row.names(A_relation_b_hm_consist)[indexes_mat[d,2]],
                                 sep = ";")
}

dependency_genepairs_KO <- matrix(NA,nrow= length(gene_pairs_nonzero),
                                 ncol = ncol(percentage_table))
row.names(dependency_genepairs_KO) <- gene_pairs_nonzero
colnames(dependency_genepairs_KO) <- colnames(percentage_table)

for(d in 1:nrow(dependency_genepairs_KO)){
  dependency_genepairs_KO[d,] <- (consistency_frame[which(consistency_frame$merged==
                                                           row.names(dependency_genepairs_KO)[d]),3] -
    consistency_frame[which(consistency_frame$merged==
                              row.names(dependency_genepairs_KO)[d]),4]) *
    sign(consistency_frame[which(consistency_frame$merged==
                                   row.names(dependency_genepairs_KO)[d]),4])
}

dim(dependency_genepairs_KO)
rank_mat <- dependency_genepairs_KO
for(d in 1:nrow(rank_mat)){
  rank_mat[d,] <- rank(rank_mat[d,])
}


gg_wt_compl_KOs <- data.frame(gene_pair = rep(row.names(dependency_genepairs_KO), each = 73),
                              KO = rep(colnames(dependency_genepairs_KO), nrow(dependency_genepairs_KO)),
                              rel_wt_compl = as.vector(t(dependency_genepairs_KO)))

col_weak_pheno <- ifelse(levels(reorder(gg_wt_compl_KOs$KO, gg_wt_compl_KOs$rel_wt_compl, mean)) %in%
                           c("Msi2", "Rps6ka1", "Arih2", "Cabin1",
                             "Hprt", "Tet1", "Igf2bp1", "Nes",
                             "Ssr2", "Irak3", "Etl4", "Pum1",
                             "Myc", "Dido1", "Hnrnph1"), "red", "black")

ggplot(gg_wt_compl_KOs, aes(x=reorder(KO, rel_wt_compl, mean), y=rel_wt_compl)) +
  geom_violin()+
  stat_summary(fun = "mean",
               geom = "crossbar",
               color= "red")+
  xlab("KO") +
  ylab("relative WT completion") +
  #coord_cartesian(ylim = c(-1.5,1.5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   color = col_weak_pheno),
        plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))


gg_wt_compl_KOs_rank <- data.frame(gene_pair = rep(row.names(rank_mat), each = 73),
                              KO = rep(colnames(rank_mat), nrow(rank_mat)),
                              rank = as.vector(t(rank_mat)))

col_weak_pheno_rank <- ifelse(levels(reorder(gg_wt_compl_KOs_rank$KO, gg_wt_compl_KOs_rank$rank, mean)) %in%
                           c("Msi2", "Rps6ka1", "Arih2", "Cabin1",
                             "Hprt", "Tet1", "Igf2bp1", "Nes",
                             "Ssr2", "Irak3", "Etl4", "Pum1",
                             "Myc", "Dido1", "Hnrnph1"), "red", "black")
ggplot(gg_wt_compl_KOs_rank, aes(x=reorder(KO, rank, mean), y=rank)) +
  geom_violin()+
  stat_summary(fun = "mean",
               geom = "crossbar",
               color= "red")+
  xlab("KO") +
  ylab("rank") +
  #coord_cartesian(ylim = c(-1.5,1.5))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   color = gg_wt_compl_KOs_rank),
        plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

### plot naive and formative dependencies #####
A_relation_b_hm_consist_naive_formative <-
  A_relation_b_hm_consist[which(row.names(A_relation_b_hm_consist) %in%
                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                    "Tbx3", "Klf4", "Prdm14", "Fgf5",
                                    "Pou3f1", "Otx2", "Dnmt3a", "Dnmt3b")),]

A_relation_b_hm_consist_naive_formative <-
  A_relation_b_hm_consist_naive_formative[,which(apply(
    A_relation_b_hm_consist_naive_formative,2, function(x) sum(x != 0)) > 0)]

Heatmap(A_relation_b_hm_consist_naive_formative,
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        cluster_rows = T,
        cluster_columns = T,
        name = "mean activity diff A to B",
        column_title = "dependencies between A (row) and B (column) \nfrom KO data\n")

A_relation_b_hm_consist_naive_formative_only <-
  A_relation_b_hm_consist_naive_formative[,which(colnames(A_relation_b_hm_consist_naive_formative) %in%
                                                   c("Nanog", "Tfcp2l1", "Esrrb",
                                                     "Tbx3", "Klf4", "Prdm14", "Fgf5",
                                                     "Pou3f1", "Otx2", "Dnmt3a", "Dnmt3b"))]

Heatmap(A_relation_b_hm_consist_naive_formative_only,
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = T,
        show_row_names = T, show_column_names = T,
        cluster_rows = T,
        cluster_columns = T,
        name = "mean activity diff A to B",
        column_title = "dependencies between A (row) and B (coconsistency_framelumn) \nfrom KO data\n")

### same for KOs with diff phenotype ####
naive_markers_ext <- c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3", "Klf4", "Prdm14", "Zfp42" ,
                       "Klf2", "Klf5", "Nr0b1")
primed_markers <- c("Dnmt3a", "Dnmt3b", "Fgf5", "Otx2", "Pou3f1")
naive_markers <- c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3", "Klf4", "Prdm14", "Zfp42")

naive_clustering <- hclust( dist(t(CLIM2$singlecontrasts$fit_N1$coefficients[naive_markers_ext,])),
                            method="ward.D2" )

col_fun_hm <- colorRamp2(c(-3,0,3),c("blue","white","red"))

hm_frame <- CLIM2$singlecontrasts$fit_N1$coefficients
colnames(hm_frame) <- sapply(X = colnames(hm_frame), FUN = function(x) strsplit(x,split="_")[[1]][2])

hm_int_naivem <- Heatmap(data.frame(hm_frame[naive_markers_ext,]),
                         col=col_fun_hm,#row_labels = "naive markers",
                         show_row_dend=F,cluster_columns=naive_clustering,
                         column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                         column_title="",name="log2FC (N24)",show_heatmap_legend=T)

hm_int_primed <- Heatmap(data.frame(hm_frame[primed_markers,]),
                         col=col_fun_hm,#row_labels = "naive markers",
                         show_row_dend=F,cluster_columns=naive_clustering,
                         column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                         column_title="",name=" ",show_heatmap_legend=F)

draw(hm_int_naivem %v% hm_int_primed, ht_gap = unit(5,"mm"))

weak_pheno <- c("Msi2", "Rps6ka1", "Arih2", "Cabin1",
                "Hprt", "Tet1", "Igf2bp1", "Nes",
                "Ssr2", "Irak3", "Etl4", "Pum1",
                "Myc", "Dido1", "Hnrnph1")

hm_form_ext <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
colnames(hm_form_ext) <- sapply(colnames(hm_form_ext), function(x) strsplit(x,split = "_")[[1]][2])

hm_form_ext_wt <- hm_form_ext/hm_form_ext[,"wt"]

form_ext_clustering <- hclust( dist(t(hm_form_ext[c(naive_markers,
                                                    primed_markers),])),
                            method="ward.D2" )

hm_naivem <- Heatmap(data.frame(hm_form_ext_wt[naive_markers,]),
                         col=col_fun_hm,#row_labels = "naive markers",
                         show_row_dend=F,cluster_columns=form_ext_clustering,
                         column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                         column_title="",name="WT completion",show_heatmap_legend=T)

hm_primed <- Heatmap(data.frame(hm_form_ext_wt[primed_markers,]),
                         col=col_fun_hm,#row_labels = "naive markers",
                         show_row_dend=F,cluster_columns=form_ext_clustering,
                         column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                         column_title="",name=" ",show_heatmap_leg=T)
                     naive_markers_ext <- c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3", "Klf4", "Prdm14", "Zfp42" ,
                                            "Klf2", "Klf5", "Nr0b1")
                     primed_markers <- c("Dnmt3a", "Dnmt3b", "Fgf5", "Otx2", "Pou3f1")
                     naive_markers <- c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3", "Klf4", "Prdm14", "Zfp42")

                     naive_clustering <- hclust( dist(t(CLIM2$singlecontrasts$fit_N1$coefficients[naive_markers_ext,])),
                                                 method="ward.D2" )

                     col_fun_hm <- colorRamp2(c(-3,0,3),c("blue","white","red"))

                     hm_frame <- CLIM2$singlecontrasts$fit_N1$coefficients
                     colnames(hm_frame) <- sapply(X = colnames(hm_frame), FUN = function(x) strsplit(x,split="_")[[1]][2])

                     hm_int_naivem <- Heatmap(data.frame(hm_frame[naive_markers_ext,]),
                                              col=col_fun_hm,#row_labels = "naive markers",
                                              show_row_dend=F,cluster_columns=naive_clustering,
                                              column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                                              column_title="",name="log2FC (N24)",show_heatmap_legend=T)

                     hm_int_primed <- Heatmap(data.frame(hm_frame[primed_markers,]),
                                              col=col_fun_hm,#row_labels = "naive markers",
                                              show_row_dend=F,cluster_columns=naive_clustering,
                                              column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                                              column_title="",name=" ",show_heatmap_legend=F)

                     draw(hm_int_naivem %v% hm_int_primed, ht_gap = unit(5,"mm"))

                     weak_pheno <- c("Msi2", "Rps6ka1", "Arih2", "Cabin1",
                                     "Hprt", "Tet1", "Igf2bp1", "Nes",
                                     "Ssr2", "Irak3", "Etl4", "Pum1",
                                     "Myc", "Dido1", "Hnrnph1")

                     hm_form_ext <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
                     colnames(hm_form_ext) <- sapply(colnames(hm_form_ext), function(x) strsplit(x,split = "_")[[1]][2])

                     hm_form_ext_wt <- hm_form_ext/hm_form_ext[,"wt"]

                     form_ext_clustering <- hclust( dist(t(hm_form_ext[c(naive_markers,
                                                                         primed_markers),])),
                                                    method="ward.D2" )

                     hm_naivem <- Heatmap(data.frame(hm_form_ext_wt[naive_markers,]),
                                          col=col_fun_hm,#row_labels = "naive markers",
                                          show_row_dend=F,cluster_columns=form_ext_clustering,
                                          column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                                          column_title="",name="WT completion",show_heatmap_legend=T)

                     hm_primed <- Heatmap(data.frame(hm_form_ext_wt[primed_markers,]),
                                          col=col_fun_hm,#row_labels = "naive markers",
                                          show_row_dend=F,cluster_columns=form_ext_clustering,
                                          column_dend_gp=gpar(col="black"),show_column_dend=T,cluster_rows=F,
                                          column_title="",name=" ",show_heatmap_legend=F)

draw(hm_naivem %v% hm_primed, ht_gap = unit(5,"mm"))

draw(hm_naivem %v% hm_primed, ht_gap = unit(5,"mm"))

## calculate relations from completion percentage with pre defined functions
relation_list_strong <- get_relations_from_percentage_genewise(percentage_table=percentage_table[,which(!colnames(percentage_table)%in%weak_pheno)])
## each list element corresponds to other to one of the following
## relation A to B (sign). Here, relations not significantly different
## from zero were already excluded
A_relation_B_matrix_strong <- relation_list_strong[[1]]
## mean relation between processes
A_mean_B_matrix_strong <- relation_list_strong[[2]]
## adjusted pval is the distribution A-B signigicantly different from 0
A_padj_B_matrix_strong <- relation_list_strong[[3]]
## pval is the distribution A minus B signigicantly different from 0
A_pval_B_matrix_strong <- relation_list_strong[[4]]
## how consisten is the directionality between A and B with the mean
consistency_frame_strong <- relation_list_strong[[5]]

## combine directionality with absolute mean relation
A_relation_b_hm_strong <- A_relation_B_matrix_strong * abs(A_mean_B_matrix_strong)

## which relations are non 0
included_relations_strong <- rep(NA, nrow(A_relation_b_hm_strong)*ncol(A_relation_b_hm_strong))
for(r in 1:dim(A_relation_b_hm_strong)[1]){
  for(c in 1:dim(A_relation_b_hm_strong)[2]){
    if(A_relation_b_hm_strong[r,c]!=0){
      included_relations_strong[as.numeric(((r-1)*nrow(A_relation_b_hm_strong)+c))] <-
        paste(row.names(A_relation_b_hm_strong)[r],
              colnames(A_relation_b_hm_strong)[c], sep=";")
    }
  }
}
included_relations_strong <- na.omit(included_relations_strong)
## extend consistency frame for relation between processes across KOs
consistency_frame_strong$merged <- paste(consistency_frame_strong$A,consistency_frame_strong$B,sep=";")
consistency_frame_strong$direction <- consistency_frame_strong$relation
consistency_frame_strong[which(consistency_frame_strong$mean<0),6] <- -consistency_frame_strong[which(consistency_frame_strong$mean<0),6]

## only summarize for non 0 relationships
summary_consistency_strong <- relation_list_strong[[6]]
row.names(summary_consistency_strong) <- summary_consistency_strong$pair
summary_consistency_strong <- as.matrix(summary_consistency_strong[which(summary_consistency_strong$pair%in%included_relations_strong),1:2])
## plot consistency of relation between process across KOs
ggplot(as.data.frame(summary_consistency_strong),aes(x=fraction_same_side)) +
  geom_histogram() +
  xlab("fraction of KOs showing same direction as mean relation")+
  ggtitle("Consistency of direction in relations") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

## only keep process process relations where more than 70 percent
## of the KOs show the same direction
consist_mat_strong <- data.frame(matrix(0, nrow = dim(A_relation_b_hm_strong)[1],
                                 ncol = dim(A_relation_b_hm_strong)[2]))
colnames(consist_mat_strong) <- colnames(A_relation_b_hm_strong)
row.names(consist_mat_strong) <- row.names(A_relation_b_hm_strong)

for(r in 1:dim(summary_consistency_strong)[1]){
  if(r%%50000==0){
    print(r)
  }
  if(summary_consistency_strong[r,1] >= 0.7){ ## cut off can be adjusted
    consist_mat_strong[strsplit(row.names(summary_consistency_strong)[r],";")[[1]][1],
                strsplit(row.names(summary_consistency_strong)[r],";")[[1]][2]] <- 1
    consist_mat_strong[strsplit(row.names(summary_consistency_strong)[r],";")[[1]][2],
                strsplit(row.names(summary_consistency_strong)[r],";")[[1]][1]] <- 1
  }
}

## trim down relations to those consistent in more than 70 percent of the KOs
A_relation_b_hm_consist_strong <- consist_mat_strong * A_relation_b_hm_strong
group <- rep(NA, nrow(A_relation_b_hm_consist_strong))
group[which(row.names(A_relation_b_hm_consist_strong)%in%
              c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                "Klf4", "Prdm14"))]  <- "naive markers"
group[which(row.names(A_relation_b_hm_consist_strong)%in%
              c("Fgf5", "Pou3f1", "Otx2",
                "Dnmt3a", "Dnmt3b"))]  <- "formative markers"

ha <- rowAnnotation(group=group,
                    col=list(group=c("naive markers"="#57BC77","formative markers"="#719AC4")))
## Fig 5.15
Heatmap(A_relation_b_hm_consist_strong,
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = hclust(dist(A_relation_b_hm_consist_strong),method="ward.D2"),
        cluster_columns = hclust(dist(A_relation_b_hm_consist_strong),method="ward.D2"),
        name = "mean activity diff A to B", left_annotation = ha)#, #right_annotation = ha,
        #column_title = "dependencies between A (row) and B (column) \nfrom KO data\n")


sum(apply(A_relation_b_hm_consist_strong, 1, function(x) x != 0))
nrow(A_relation_b_hm_consist_strong)*ncol(A_relation_b_hm_consist_strong)

### plot naive and formative dependencies #####
A_relation_b_hm_consist_naive_formative_strong <-
  A_relation_b_hm_consist_strong[which(row.names(A_relation_b_hm_consist_strong) %in%
                                  c("Nanog", "Tfcp2l1", "Esrrb",
                                    "Tbx3", "Klf4", "Prdm14", "Fgf5",
                                    "Pou3f1", "Otx2", "Dnmt3a", "Dnmt3b")),]

A_relation_b_hm_consist_naive_formative_strong <-
  A_relation_b_hm_consist_naive_formative_strong[,which(apply(
    A_relation_b_hm_consist_naive_formative_strong,2, function(x) sum(x != 0)) > 0)]
## Fig 5.17
Heatmap(A_relation_b_hm_consist_naive_formative_strong,
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        cluster_rows = T,
        cluster_columns = T,
        name = "mean activity diff A to B",
        column_title = "dependencies between A (row) and B (column) \nfrom KO data\n")

A_relation_b_hm_consist_naive_formative_only_strong <-
  A_relation_b_hm_consist_naive_formative_strong[,which(colnames(A_relation_b_hm_consist_naive_formative_strong) %in%
                                                   c("Nanog", "Tfcp2l1", "Esrrb",
                                                     "Tbx3", "Klf4", "Prdm14", "Fgf5",
                                                     "Pou3f1", "Otx2", "Dnmt3a", "Dnmt3b"))]
## Fig 5.18
Heatmap(A_relation_b_hm_consist_naive_formative_only_strong,
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = T,
        show_row_names = T, show_column_names = T,
        cluster_rows = T,
        cluster_columns = T,
        name = "mean activity diff A to B",
        column_title = "dependencies between A (row) and B (column) \nfrom KO data\n")

### compare strong phenos results to all ####
shared_n <- sum(row.names(summary_consistency_strong)%in%row.names(summary_consistency))
strong_only <- sum(!row.names(summary_consistency_strong)%in%row.names(summary_consistency))
all_only <- sum(!row.names(summary_consistency)%in%row.names(summary_consistency_strong))

summary_compare <- as.data.frame(matrix(0,nrow=shared_n+strong_only+all_only,
                          ncol=4))
colnames(summary_compare) <- c("consistency_all","mean_all",
                               "consistency_strong","mean_strong")

summary_compare[1:shared_n,1:2] <-
  summary_consistency[which(
  row.names(summary_consistency)%in%
    row.names(summary_consistency_strong)),1:2]

summary_compare[1:shared_n,3:4] <-
  summary_consistency_strong[which(
  row.names(summary_consistency_strong)%in%
    row.names(summary_consistency)),1:2]

row.names(summary_compare[1:shared_n,]) <-
  row.names(summary_consistency)[which(
    row.names(summary_consistency)%in%
      row.names(summary_consistency_strong))]

summary_compare[(shared_n+1):(shared_n+strong_only),3:4] <-
  summary_consistency_strong[which(
  !row.names(summary_consistency_strong)%in%
    row.names(summary_consistency)),1:2]

row.names(summary_compare[(shared_n+1):(shared_n+strong_only),]) <-
  row.names(summary_consistency_strong)[which(
    !row.names(summary_consistency_strong)%in%
      row.names(summary_consistency))]

summary_compare[(shared_n+strong_only+1):(shared_n+strong_only+all_only),1:2] <-
  summary_consistency[which(
  !row.names(summary_consistency)%in%
    row.names(summary_consistency_strong)),1:2]

row.names(summary_compare[(shared_n+strong_only+1):(shared_n+strong_only+all_only),]) <-
  row.names(summary_consistency)[which(
    !row.names(summary_consistency)%in%
      row.names(summary_consistency_strong))]

library(ggExtra)
## Fig 5.14 A
ggplot(summary_compare, aes(x=mean_all,y=mean_strong))+
  geom_bin_2d(bins=85)+
  xlab("effect size based on all KOs")+
  ylab("effect size based on KOs with strong phenotype")+
  ggtitle("comparing effect size based on strong and all phenotypes") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

## Fig 5.14 B
ggplot(summary_compare, aes(x=consistency_all,y=consistency_strong))+
  geom_bin_2d(bins=58)+
  xlab("consistency based on all KOs")+
  ylab("consistency based on KOs with strong phenotype")+
  ggtitle("comparing consistency based on strong and all phenotypes") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))

p <- ggplot(summary_compare, aes(x=consistency_all,y=consistency_strong))+
  geom_point(alpha=0) +
  geom_bin_2d(bins=55) +
  scale_fill_gradient(trans= "log10") +
  xlab("consistency based on all KOs")+
  ylab("consistency based on KOs with strong phenotype")+
  ggtitle("comparing consistency based on strong and all phenotypes") +
  coord_fixed(ratio=1)+
  geom_abline(slope = 1, intercept = 0, color="red")+
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"),
        legend.position = "left")

ggMarginal(p, type="density", bw=0.05)

summary_dens <- data.frame(consistency = c(summary_compare[,1],summary_compare[,3]),
                           mean = c(summary_compare[,2],summary_compare[,4]),
                           type = c(rep("all",nrow(summary_compare)),
                                    rep("strong phenotypes",nrow(summary_compare))))

ggplot(summary_dens)+
  geom_density(aes(x=consistency, col = type, fill = type),
               alpha=0.5,bw=0.05) +
  xlab("consistency of KOs")+
  ggtitle("comparing consistency based on strong and all phenotypes") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))
## Fig 5.13A
percentage_plot(percentage_table = percentage_table,
                Gene_A = "Dnmt3a",
                Gene_B = "Nanog")
## Fig 5.13
percentage_plot(percentage_table = percentage_table,
                Gene_A = "Pou3f1",
                Gene_B = "Klf4")


hist_pheno_test <- pheno_hist_test(percentage_table = percentage_table,
                                   dependency_mat = A_relation_b_hm_consist_strong)

summary(hist_pheno_test[which(hist_pheno_test$type==1),1])
summary(hist_pheno_test[which(hist_pheno_test$type==0),1])

ks.test(hist_pheno_test[which(hist_pheno_test$type==1),1],
        hist_pheno_test[which(hist_pheno_test$type==0),1],
        alternative = "two.sided")

ggplot(hist_pheno_test) +
  geom_density(data=subset(hist_pheno_test, type == 1),
               aes(x=position), fill="blue",
               colour="blue", alpha=0.3) +
  geom_density(data=subset(hist_pheno_test, type == 0),
               aes(x=position), fill="red",
               colour="red", alpha=0.3) +
  xlab("position in comparison to mean of WT completion differences")+
  ggtitle("\n") +
  theme(plot.title=element_text(hjust=0.5),
        panel.background = element_rect(fill=NA),
        panel.grid.major = element_line(colour = "grey90"))



#### NMF ####
library(foreach)
hclust_relations_unscaled <- hclust(dist(A_relation_b_hm_consist_strong),method="ward.D2")
## Roberts function to retrieve withtin cluster variance
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
WCV_cluster_size = foreach(i = 3:50, .combine=c) %do% {
  rs.withinClusterVariance( A_relation_b_hm_consist_strong,
                            cutree(hclust_relations_unscaled,i) )
}
## plot change of within cluster variance to find good k
plot(diff(WCV_cluster_size) ~ c(4:50), type="b", xlim=c(0,50),
     xlab="number of clusters (k)", main="clusters of processes",
     ylab=expression( W[k] ~ "-" ~ W[k-1]), pch=17, cex=1.5, col="black" )
abline(v=c(10,12), lwd=2, lty=2)
number_clusters <- 11

library(NMF)

pos_mat <- abs(A_relation_b_hm_consist_strong)
pos_mat_2 <- A_relation_b_hm_consist_strong
pos_mat_2[pos_mat_2 < 0 ] <- 0

connectivity <- apply(pos_mat, 1, function(x) sum(x != 0))
more_than_one <- which(connectivity>1)

pos_mat_rand <- randomize(pos_mat)

test <- nmf(x = pos_mat, rank = 5:12,
            nrun = 50,seed = 123,
            method = "lee", .opt = "vp4")
test_rand <- nmf(x = pos_mat_rand, rank = 5:12,
                 nrun = 50,seed = 123,
                 method = "lee", .opt = "vp4")
plot(test,test_rand)

nmf_rank <- 6
nmf_res <- nmf(x = pos_mat, rank = nmf_rank,
               nrun = 200,seed = 1234,
               method = "lee", .opt = "vp4")


hierarchichal_clust_nmf <- consensushc(nmf_res, dendrogram = F)
clusters_hclust_nmf <- data.frame(cutree(hierarchichal_clust_nmf,k=nmf_rank))

W_nmf <- nmf_res@fit@W
H_nmf <- nmf_res@fit@H
centroids_w <- data.frame(row.names(W_nmf)[apply(W_nmf, 2, function(x) which(x==max(x)))])
clusters_h <- data.frame(apply(H_nmf, 2, function(x) which(x==max(x))))
clusters_h$gene <- row.names(clusters_h)
colnames(clusters_h) <- c("cluster", "gene")

H_norm <- apply(H_nmf, 2, function(x) x/max(x))

#W_norm <- apply(W_nmf, 2, function(x) x/sum(x))
#S <- diag(1, nrow=nrow(W_norm), ncol=nrow(W_norm))

clusters_nmf <- data.frame(gene = row.names(clusters_h),
                           cluster_h = clusters_h[,1],
                           centroids_w = NA,
                           clusters_hclust=clusters_hclust_nmf[,1])
for(d in 1:nrow(centroids_w)){
  clusters_nmf[which(clusters_nmf[,2]==d),3] <- centroids_w[d,1]
}

NMF_list <- list(Clusters=clusters_nmf,
                 H_nmf=H_nmf,
                 W_nmf=W_nmf)

group <- rep(NA, nrow(A_relation_b_hm_consist_strong))
names(group) <- names(clusters_h)
group[which(row.names(A_relation_b_hm_consist_strong)%in%
              c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                "Klf4", "Prdm14"))] <- "naive markers"
group[which(row.names(A_relation_b_hm_consist_strong)%in%
              c("Fgf5", "Pou3f1", "Otx2",
                "Dnmt3a", "Dnmt3b"))]  <- "formative markers"


ha_row <- rowAnnotation(group=group,
                        col=list(group=c("naive markers"="#57BC77",
                                         "formative markers"="#719AC4")))

ha_col <- columnAnnotation(group=group,
                           col=list(group=c("naive markers"="#57BC77",
                                            "formative markers"="#719AC4")),
                           show_legend=F)


# Heatmap(A_relation_b_hm_consist_strong,
#         col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
#         show_row_dend = T, show_column_dend = F,
#         show_row_names = F, show_column_names = F,
#         cluster_rows = hclust_relations_unscaled,
#         cluster_columns = hclust_relations_unscaled,
#         row_split = nmf_rank,
#         row_title_rot = 0,NMF_list
#         row_title_gp = gpar(fontsize = 8),
#         row_gap = unit(1,"mm"),
#         column_split = nmf_rank,
#         column_title_rot = 0,
#         column_title_gp = gpar(fontsize = 8),
#         column_gap = unit(1,"mm"),
#         name = "mean activity diff A to B", left_annotation = ha_row, top_annotation = ha_col,
#         column_title = paste("consistent dependencies between A (row) and B (column) \n ",
#                              nmf_rank, " hierarchichal clusters\n",sep=""))

ht = Heatmap(as.matrix(A_relation_b_hm_consist_strong),
        col=colorRamp2(c(-1.5,0,1.5), c("blue",'white',"red")),
        show_row_dend = F, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = F,
        cluster_columns = F,
        row_split = factor(clusters_h$cluster),
        cluster_row_slices = T,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 8),
        row_gap = unit(1,"mm"),
        column_split = factor(clusters_h$cluster),
        cluster_column_slices = T,
        column_title_rot = 0,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(1,"mm"),
        name = "red: row independent of column\nblue: column independent of row",
        left_annotation = ha_row,
        top_annotation = ha_col)

draw(ht, column_title = paste("consistent dependencies between A (row) and B (column) \n ",
                                 nmf_rank, " NMF clusters\n",sep=""))

### export objects ####

saveRDS(A_relation_b_hm_consist_strong, "RDS/A_relation_b_hm_consist_perc_genes.rds")
saveRDS(summary_consistency_strong, "RDS/summary_consistency_perc_genes.rds")
saveRDS(percentage_table,"RDS/percentage_table_genes.rds")
saveRDS(enrTblListFull, "RDS/enrTblListFull.rds")
saveRDS(clusters_h, "RDS/clusters_h.rds")
saveRDS(NMF_list, "RDS/NMF_list.rds")



