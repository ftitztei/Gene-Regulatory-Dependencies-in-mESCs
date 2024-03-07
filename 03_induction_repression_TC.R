library(ggplot2)
library(GO.db)
library(reactome.db)
library(Matrix)
library(cluster)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

library(Hmisc)
library(biomaRt)
library(topGO)
library(org.Mm.eg.db)
library(Mus.musculus)

source("R/functions.R")

### load objects from previous scripts ####
tpms_TC_gpr_shrunken <- readRDS("RDS/tpms_TC_gpr.rds")
gpr_list_shrunkenFCS <- readRDS("RDS/gpr_list_shrunkenFCS.rds")
shrunken_foldchanges_complete <- readRDS("RDS/shrunken_foldchanges_complete.rds")

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol"),
                 mart=ensembl)
#biomart <- biomart[which(!is.na(biomart$entrezgene)),]
biomart <- biomart[complete.cases(biomart[,c(1,3)]),]

### get slopes and 2nd derivative ####
TPMs_dev_max_shrunken <- as.matrix(tpms_TC_gpr_shrunken[,c(1,2,5:21)])
TPMs_dev_max_shrunken <- TPMs_dev_max_shrunken/rowMax(TPMs_dev_max_shrunken)

scaled_tpms_TC_gpr_shrunken <- as.matrix(tpms_TC_gpr_shrunken[,c(1,2,5:21)])


slopes_frame_TC <- data.frame(matrix(NA,nrow=dim(scaled_tpms_TC_gpr_shrunken)[1],ncol=16))
colnames(slopes_frame_TC) <- c("slope_0_2", "slope_2_4", "slope_4_6", "slope_6_8",
                               "slope_8_10", "slope_10_12", "slope_12_14", "slope_14_16",
                               "slope_16_18", "slope_18_20", "slope_20_22", "slope_22_24",
                               "slope_24_26", "slope_26_28", "slope_28_30", "slope_30_32")
row.names(slopes_frame_TC) <- row.names(scaled_tpms_TC_gpr_shrunken)
slopes_frame_TC <- as.matrix(slopes_frame_TC)
for(r in 1:dim(slopes_frame_TC)[1]){
  for(c in 1:dim(slopes_frame_TC)[2]){
    if(c == 1){
      slopes_frame_TC[r,c] <- log2(((scaled_tpms_TC_gpr_shrunken[r,3]+0.01)/(mean(scaled_tpms_TC_gpr_shrunken[r,1:2])+0.01)))
    }
    if(c == 16){
      slopes_frame_TC[r,c] <- log2(((mean(scaled_tpms_TC_gpr_shrunken[r,18:19])+0.01)/(scaled_tpms_TC_gpr_shrunken[r,17]+0.01)))
    }
    else{
      slopes_frame_TC[r,c] <- log2(((scaled_tpms_TC_gpr_shrunken[r,c+2]+0.01)/(scaled_tpms_TC_gpr_shrunken[r,c+1]+0.01)))
    }
  }
}

second_derivate <- data.frame(matrix(NA,nrow=dim(scaled_tpms_TC_gpr_shrunken)[1],
                                     ncol=16))
colnames(second_derivate) <- c("0h", "2h", "4h", "6h", "8h", "10h",
                               "12h", "14h", "16h", "18h", "20h",
                               "22h", "24h", "26h", "28h", "30h")
row.names(second_derivate) <- row.names(scaled_tpms_TC_gpr_shrunken)
second_derivate <- as.matrix(second_derivate)

for(r in 1:dim(second_derivate)[1]){
  if(0>slopes_frame_TC[r,1]){
    second_derivate[r,1] <- -dist(c(0,slopes_frame_TC[r,1]))
  } else {
    second_derivate[r,1] <- dist(c(0,slopes_frame_TC[r,1]))
  }

  for(c in 2:dim(second_derivate)[2]){
    if(slopes_frame_TC[r,c-1]>slopes_frame_TC[r,c]){
      second_derivate[r,c] <- -dist(c(slopes_frame_TC[r,c-1],
                                      slopes_frame_TC[r,c]))
    } else {
      second_derivate[r,c] <- dist(c(slopes_frame_TC[r,c-1],
                                     slopes_frame_TC[r,c]))
    }

  }
}


### calculate induction and repression points ####
gpr_list_shrunken_mean <- gpr_list_shrunkenFCS[[1]][,2:18]
for(d in 1:dim(gpr_list_shrunken_mean)[1]){
  gpr_list_shrunken_mean[d,] <- gpr_list_shrunken_mean[d,]-gpr_list_shrunken_mean[d,1]
}

changed_genes <- row.names(gpr_list_shrunken_mean[which(rowMax(gpr_list_shrunken_mean)-
                                                          rowMin(gpr_list_shrunken_mean)>=0.5),])
changed_genes <- changed_genes[which(changed_genes%in%
                                       row.names(tpms_TC_gpr_shrunken[which(tpms_TC_gpr_shrunken$max>=5),]))]

induct_repr_simple <- get_regulation_from_runs(gene_names_cutoff = changed_genes,
                                               slopes_frame_TC=slopes_frame_TC,
                                               log2FC_frame_TC=gpr_list_shrunken_mean,
                                               thresh_0=0.1)

ggplot_hist_simple_repr_induct <- data.frame(gene=character(),
                                             time_point=numeric(),
                                             type=character())

for(l in 1:length(induct_repr_simple)){
  if(!is.na(induct_repr_simple[[l]])[1]){
    for(d in 1:length(induct_repr_simple[[l]])){
      if(induct_repr_simple[[l]][d]>0){
        ggplot_hist_simple_repr_induct <- rbind(ggplot_hist_simple_repr_induct,
                                             data.frame(gene=names(induct_repr_simple)[l],
                                                        time_point=(as.numeric(names(induct_repr_simple[[l]])[d])),
                                                        type="induction"))
      } else {
        ggplot_hist_simple_repr_induct <- rbind(ggplot_hist_simple_repr_induct,
                                             data.frame(gene=names(induct_repr_simple)[l],
                                                        time_point=(as.numeric(names(induct_repr_simple[[l]])[d])),
                                                        type="repression"))
      }

    }
  }
}

ggplot(ggplot_hist_simple_repr_induct,aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on slope runs") +
  ggtitle("Histogram of induction or repression TPs\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme(plot.title=element_text(hjust=0.5))

length(unique(ggplot_hist_simple_repr_induct$gene))

hist(lengths(induct_repr_simple))

more_than_2_simple <- names(induct_repr_simple)[which(lengths(induct_repr_simple)>2)]
length(more_than_2_simple)

ggplot_hist_simple_repr_induct_short <- ggplot_hist_simple_repr_induct[which(!ggplot_hist_simple_repr_induct$gene%in%more_than_2_simple),]

ggplot(ggplot_hist_simple_repr_induct_short,aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on slope runs (1 or 2 runs)") +
  ggtitle("Histogram of induction or repression TPs\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme(plot.title=element_text(hjust=0.5))



changed_genes_not_simple <- changed_genes[which(!changed_genes%in%ggplot_hist_simple_repr_induct_short$gene)]

induct_repr_results <- get_induction_repression(second_derivate = second_derivate,
                                                slopes_frame_TC = slopes_frame_TC,
                                                gene_names_cutoff = changed_genes_not_simple,
                                                MinMaxFrame = scaled_tpms_TC_gpr_shrunken,
                                                MinMaxType = "TPM",
                                                dynamic_thresh_perc = 0.9411765) #(16/17 to get from 17 to 16 as I used -1 before)

ggplot_hist_all_repr_induct <- data.frame(gene=character(),
                                          time_point=numeric(),
                                          type=character())


for(l in 1:length(induct_repr_results)){
  if(!is.na(induct_repr_results[[l]])[1]){
    for(d in 1:length(induct_repr_results[[l]])){
      if(induct_repr_results[[l]][d]>0){
        ggplot_hist_all_repr_induct <- rbind(ggplot_hist_all_repr_induct,
                                             data.frame(gene=names(induct_repr_results)[l],
                                                        time_point=(as.numeric(names(induct_repr_results[[l]])[d])),
                                                        type="induction"))
      } else {
        ggplot_hist_all_repr_induct <- rbind(ggplot_hist_all_repr_induct,
                                             data.frame(gene=names(induct_repr_results)[l],
                                                        time_point=(as.numeric(names(induct_repr_results[[l]])[d])),
                                                        type="repression"))
      }

    }
  }
}

ggplot(ggplot_hist_all_repr_induct,aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of induction or repression TPs\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme(plot.title=element_text(hjust=0.5))

length(unique(ggplot_hist_all_repr_induct$gene))
hist(lengths(induct_repr_results))

more_than_3 <- names(induct_repr_results)[which(lengths(induct_repr_results)>3)]
length(more_than_3)

if(0){
  for(d in 1:length(more_than_3)){
  plot(plotGprKinetics(gene_name = more_than_3[d],
                  gpr_list = gpr_list_shrunkenFCS,
                  mgi=F, unfitted = shrunken_foldchanges_complete,
                  type="shrunken",
                  induction_repression = ggplot_hist_all_repr_induct[which(ggplot_hist_all_repr_induct$gene==
                                                                             more_than_3[d]),]))

    }
}


ggplot_hist_all_repr_induct_short <- ggplot_hist_all_repr_induct[which(!ggplot_hist_all_repr_induct$gene%in%more_than_3),]

ggplot(ggplot_hist_all_repr_induct_short,aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of induction or repression TPs\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme(plot.title=element_text(hjust=0.5))

ggplot_hist_repr_induct_complete <- rbind(ggplot_hist_simple_repr_induct_short,
                                          data.frame(gene=ggplot_hist_all_repr_induct_short$gene,
                                                     time_point=ggplot_hist_all_repr_induct_short$time_point,
                                                     type=ggplot_hist_all_repr_induct_short$type))

ggplot(ggplot_hist_repr_induct_complete,aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on runs and 2nd derivative") +
  ggtitle("Histogram of induction or repression TPs\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))


### TPs of regulation 0 hrs, 6 hrs late repression####
induction_0 <- as.character(ggplot_hist_repr_induct_complete[which(ggplot_hist_repr_induct_complete$time_point==0&
                                                                     ggplot_hist_repr_induct_complete$type=="induction"),1])

induction_6 <- as.character(ggplot_hist_repr_induct_complete[which(ggplot_hist_repr_induct_complete$time_point==6&
                                                                ggplot_hist_repr_induct_complete$type=="induction"),1])

induction_late <- as.character(ggplot_hist_repr_induct_complete[which(ggplot_hist_repr_induct_complete$time_point>=12&
                                                                        ggplot_hist_repr_induct_complete$type=="induction"),1])


repression_0 <- as.character(ggplot_hist_repr_induct_complete[which(ggplot_hist_repr_induct_complete$time_point==0&
                                                                 ggplot_hist_repr_induct_complete$type=="repression"),1])

repression_6 <- as.character(ggplot_hist_repr_induct_complete[which(ggplot_hist_repr_induct_complete$time_point==6&
                                                                      ggplot_hist_repr_induct_complete$type=="repression"),1])

repression_late <- as.character(ggplot_hist_repr_induct_complete[which(ggplot_hist_repr_induct_complete$time_point>=12&
                                                                    ggplot_hist_repr_induct_complete$type=="repression"),1])


sum(repression_0%in%induction_6)/length(induction_6)

sum(induction_0%in%repression_6)/length(repression_6)

groups_table <- data.frame(gene=c(induction_0,
                                  induction_6,
                                  induction_late,
                                  repression_0,
                                  repression_6,
                                  repression_late),
                           group=(c(rep("induction_0",length(induction_0)),
                                   rep("induction_6",length(induction_6)),
                                   rep("induction_late",length(induction_late)),
                                   rep("repression_0",length(repression_0)),
                                   rep("repression_6",length(repression_6)),
                                   rep("repression_late",length(repression_late)))))

## plot individual group histograms
## Fig 5.5A
ggplot(data= subset(ggplot_hist_repr_induct_complete, gene %in% induction_0),aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of genes induced at 0 hours\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
## Fig 5.5C
ggplot(data= subset(ggplot_hist_repr_induct_complete, gene %in% induction_6),aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of genes induced at 6 hours\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))

ggplot(data= subset(ggplot_hist_repr_induct_complete, gene %in% induction_late),aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of genes induced from 12 to 32 hours\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))


## Fig 5.5 B
ggplot(data= subset(ggplot_hist_repr_induct_complete, gene %in% repression_0),aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of genes repressed at 0 hours\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
## Fig 5.5 D
ggplot(data= subset(ggplot_hist_repr_induct_complete, gene %in% repression_6),aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of genes repressed at 6 hours\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))

ggplot(data= subset(ggplot_hist_repr_induct_complete, gene %in% repression_late),aes(x = time_point)) +
  geom_histogram(binwidth = 2,alpha=0.7) +
  xlim(c(-1,32)) +
  xlab("timing based on 2nd derivative") +
  ggtitle("Histogram of genes repressed from 12 to 32 hours\nshrunken FCs") +
  facet_grid(. ~ type) +
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))

### add ribbon plots ####
groups <- list(induction_0hrs = induction_0,
               induction_6hrs = induction_6,
               induction_late = induction_late,
               repression_0hrs = repression_0,
               repression_6hrs = repression_6,
               repression_late = repression_late)

gg_TC_groups <- data.frame(name=character(),
                           group=character(),
                           time=numeric(),
                           log2FC=numeric(),
                           type=character())

for(g in 1:length(groups)){
  holder <- data.frame(name=character(),
                       group=character(),
                       time=numeric(),
                       log2FC=numeric(),
                       type=character())

  for(d in 1:length(groups[[g]])){
    if(names(groups)[g]%in%c("induction_0hrs","induction_6hrs",
                             "induction_late")){
      holder <- rbind(holder,
                      data.frame(name=rep(groups[[g]][d],17),
                                 group=rep(names(groups)[g],17),
                                 time=seq(0,32,2),
                                 log2FC=gpr_list_shrunkenFCS[[1]][groups[[g]][d],2:18]-
                                   gpr_list_shrunkenFCS[[1]][groups[[g]][d],2],
                                 type="induction"))
    }else{
      holder <- rbind(holder,
                      data.frame(name=rep(groups[[g]][d],17),
                                 group=rep(names(groups)[g],17),
                                 time=seq(0,32,2),
                                 log2FC=gpr_list_shrunkenFCS[[1]][groups[[g]][d],2:18]-
                                   gpr_list_shrunkenFCS[[1]][groups[[g]][d],2],
                                 type="repression"))
    }
  }

  gg_TC_groups <- rbind(gg_TC_groups, holder)
}

#delete object mean from workspace to make it work
##Fig 5.6
ggplot(gg_TC_groups, aes(x=time, y=log2FC,
                         fill=group, col=group)) +
  stat_summary(fun.data=mean_cl_normal,geom="ribbon",
               fun.args=list(conf.int=0.95),alpha=0.2) +
  stat_summary(fun=mean,size=1,geom="line") +
  scale_color_manual(values = brewer.pal(11, "Spectral")[c(1,8,10,11,5,2)]) +
  scale_fill_manual(values = brewer.pal(11, "Spectral")[c(1,8,10,11,5,2)]) +
  theme_bw()+
  facet_grid(type~.)

### cluster gene groups to get more homogeneous groups and insight to stress reaction ####
clust_instant_ind <- hclust(dist(gpr_list_shrunkenFCS[[1]][induction_0,2:18]-
                                   gpr_list_shrunkenFCS[[1]][induction_0,2]))
clust_instant_ind_cut <- cutree(clust_instant_ind,k=2)

gg_TC_clust_instant_ind <- data.frame(name=character(),
                                      time=numeric(),
                                      log2FC=numeric(),
                                      cluster=character())

for(d in 1:length(clust_instant_ind_cut)){
  gg_TC_clust_instant_ind <- rbind(gg_TC_clust_instant_ind,
                                   data.frame(name=rep(names(clust_instant_ind_cut)[d],17),
                                              time=seq(0,32,2),
                                              log2FC=gpr_list_shrunkenFCS[[1]][names(clust_instant_ind_cut)[d],2:18]-
                                                gpr_list_shrunkenFCS[[1]][names(clust_instant_ind_cut)[d],2],
                                              cluster=rep(as.character(clust_instant_ind_cut)[d],17)))
}

#delete object mean from workspace to make it work
ggplot(gg_TC_clust_instant_ind, aes(x=time, y=log2FC,
                                    fill=cluster, col=cluster)) +
  geom_line(data=gg_TC_clust_instant_ind, aes(x=time, y=log2FC,col=cluster,
                                              group=name), alpha=0.1) +
  stat_summary(fun.data=mean_cl_normal,geom="ribbon",
               fun.args=list(conf.int=0.95),alpha=0.5) +
  stat_summary(fun=mean,size=1,geom="line") +
  theme_bw()+
  geom_hline(yintercept = c(-0.5,0.5), lty=2)

summary(as.factor(clust_instant_ind_cut))

###
clust_6hrs_repr <- hclust(dist(gpr_list_shrunkenFCS[[1]][repression_6,2:18]-
                                 gpr_list_shrunkenFCS[[1]][repression_6,2]))
clust_6hrs_repr_cut <- cutree(clust_6hrs_repr,k=2)

gg_TC_clust_6hrs_repr <- data.frame(name=character(),
                                    time=numeric(),
                                    log2FC=numeric(),
                                    cluster=character())

for(d in 1:length(clust_6hrs_repr_cut)){
  gg_TC_clust_6hrs_repr <- rbind(gg_TC_clust_6hrs_repr,
                                 data.frame(name=rep(names(clust_6hrs_repr_cut)[d],17),
                                            time=seq(0,32,2),
                                            log2FC=gpr_list_shrunkenFCS[[1]][names(clust_6hrs_repr_cut)[d],2:18]-
                                              gpr_list_shrunkenFCS[[1]][names(clust_6hrs_repr_cut)[d],2],
                                            cluster=rep(as.character(clust_6hrs_repr_cut)[d],17)))
}


ggplot(gg_TC_clust_6hrs_repr, aes(x=time, y=log2FC,
                                  fill=cluster, col=cluster)) +
  geom_line(data=gg_TC_clust_6hrs_repr, aes(x=time, y=log2FC,col=cluster,
                                            group=name), alpha=0.1) +
  stat_summary(fun.data=mean_cl_normal,geom="ribbon",
               fun.args=list(conf.int=0.95),alpha=0.5) +
  stat_summary(fun=mean,size=1,geom="line") +
  geom_hline(yintercept = c(-0.5,0.5), lty=2)

summary(as.factor(clust_6hrs_repr_cut))

#### heatmap and clustering of all genes that change over time ####
gpr_list_shrunken_mean_changed <- gpr_list_shrunken_mean[changed_genes,]

ztrans_log2FCs <- gpr_list_shrunken_mean_changed
for(d in 1:nrow(gpr_list_shrunken_mean_changed)){
  ztrans_log2FCs[d,] <- ztrans(gpr_list_shrunken_mean_changed[d,])
}

ztrans_log2FCs_all <- gpr_list_shrunken_mean
for(d in 1:nrow(gpr_list_shrunken_mean)){
  ztrans_log2FCs_all[d,] <- ztrans(gpr_list_shrunken_mean[d,])
}

Heatmap(ztrans_log2FCs,
        col=colorRamp2(c(-5,0,5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = T,
        cluster_columns = F,
        name = "log2FC vs 2i",
        column_title = "z transformed log2FCs vs 2i\nchanged genes\n")

withinClusterVariance = function( m, clus ) {
  withinvar <- numeric()
  for(c in unique(clus)){
    centroid <- apply(m[which(clus==c),],2,mean)
    clus_var <- sum( apply(m[which(clus==c),],1,function(x) x - centroid)^2 )
    withinvar <- c(withinvar,clus_var)
  }

  return( sum(withinvar) )
}

hclust_ztrans_unscaled <- hclust(dist(ztrans_log2FCs),
                                    method="ward.D2")

WCV_cluster_size = foreach(i = 3:30, .combine=c) %do% {
  withinClusterVariance( ztrans_log2FCs,
                            cutree(hclust_ztrans_unscaled,i) )
}

## plot change of within cluster variance to find good k
plot(diff(WCV_cluster_size) ~ c(4:30), type="b", xlim=c(0,30),
     xlab="number of clusters (k)", main="clusters of processes",
     ylab=expression( W[k] ~ "-" ~ W[k-1]), pch=17, cex=1.5, col="black" )
abline(v=c(9), lwd=2, lty=2)

clusters_hclust_unscaled <- data.frame(cutree(hclust_ztrans_unscaled,k=9))

Heatmap(ztrans_log2FCs,
        col=colorRamp2(c(-5,0,5), c("blue",'white',"red")),
        show_row_dend = T, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        row_split = clusters_hclust_unscaled,
        cluster_rows = T,
        cluster_columns = F,
        use_raster = F,
        name = "log2FC vs 2i",
        column_title = "z transformed log2FCs vs 2i\nchanged genes\n")

gg_clusters <- data.frame(log2FC=c(ztrans_log2FCs),
                          time=c(sapply(seq(0,32,2), function(x) rep(x,nrow(ztrans_log2FCs)))),
                          gene=row.names(ztrans_log2FCs),
                          cluster=clusters_hclust_unscaled$cutree.hclust_ztrans_unscaled..k...9.)
## Fig 5.3
ggplot(gg_clusters, aes(x=time, y=log2FC,
                    group=gene)) +
  geom_line(alpha=0.1, color="grey40") +
  facet_wrap(~cluster) +
  ylab("ztransformed log2FCs") +
  theme_bw() +
  stat_summary(data=gg_clusters, aes(x=time, y=log2FC,
                                group=cluster,
                                color=as.factor(cluster),
                                fill=as.factor(cluster)),
               fun.data=mean_cl_normal,geom="ribbon",
               fun.args=list(conf.int=0.95),alpha=0.5) +
  scale_color_manual(values = brewer.pal(11, "Spectral")[c(2,7,10,3,5,1,9,8,11)]) +
  scale_fill_manual(values = brewer.pal(11, "Spectral")[c(2,7,10,3,5,1,9,8,11)])

out_frame <- clusters_hclust_unscaled
out_frame$mgi <- NA
out_frame$ensemble <- row.names(out_frame)
for(d in 1:nrow(out_frame)){
  out_frame[d,2] <- biomart[which(biomart$ensembl_gene_id==out_frame[d,3]),3][1]
}
colnames(out_frame) <- c("cluster","mgi","ensembl")

write.csv(out_frame, "output_tables/kinetic_clusters_201223.csv",
          quote = F,row.names = F)
### redo ribbon plot on ztransformed ####

groups <- list(induction_0hrs = induction_0[which(induction_0%in%
                                                    row.names(ztrans_log2FCs))],
               induction_6hrs = induction_6[which(induction_6%in%
                                                    row.names(ztrans_log2FCs))],
               induction_late = induction_late[which(induction_late%in%
                                                       row.names(ztrans_log2FCs))],
               repression_0hrs = repression_0[which(repression_0%in%
                                                      row.names(ztrans_log2FCs))],
               repression_6hrs = repression_6[which(repression_6%in%
                                                      row.names(ztrans_log2FCs))],
               repression_late = repression_late[which(repression_late%in%
                                                         row.names(ztrans_log2FCs))])

gg_TC_groups <- data.frame(name=character(),
                           group=character(),
                           time=numeric(),
                           log2FC=numeric(),
                           type=character())

for(g in 1:length(groups)){
  holder <- data.frame(name=character(),
                       group=character(),
                       time=numeric(),
                       log2FC=numeric(),
                       type=character())

  for(d in 1:length(groups[[g]])){
    if(names(groups)[g]%in%c("induction_0hrs","induction_6hrs",
                             "induction_late")){
      holder <- rbind(holder,
                      data.frame(name=rep(groups[[g]][d],17),
                                 group=rep(names(groups)[g],17),
                                 time=seq(0,32,2),
                                 log2FC=ztrans_log2FCs[groups[[g]][d],],
                                 type="induction"))
    }else{
      holder <- rbind(holder,
                      data.frame(name=rep(groups[[g]][d],17),
                                 group=rep(names(groups)[g],17),
                                 time=seq(0,32,2),
                                 log2FC=ztrans_log2FCs[groups[[g]][d],],
                                 type="repression"))
    }
  }

  gg_TC_groups <- rbind(gg_TC_groups, holder)
}

#delete object mean from workspace to make it work
ggplot(gg_TC_groups, aes(x=time, y=log2FC,
                         fill=group, col=group)) +
  stat_summary(fun.data=mean_cl_normal,geom="ribbon",
               fun.args=list(conf.int=0.95),alpha=0.2) +
  stat_summary(fun=mean,size=1,geom="line") +
  scale_color_manual(values = brewer.pal(11, "Spectral")[c(1,8,10,11,5,2)]) +
  scale_fill_manual(values = brewer.pal(11, "Spectral")[c(1,8,10,11,5,2)]) +
  theme_bw()+
  facet_grid(type~.)

### export objects ####
saveRDS(ggplot_hist_all_repr_induct,"RDS/ggplot_hist_all_repr_induct.rds")
saveRDS(ggplot_hist_all_repr_induct_short,"RDS/ggplot_hist_all_repr_induct_short.rds")

saveRDS(ggplot_hist_simple_repr_induct,"RDS/ggplot_hist_simple_repr_induct.rds")
saveRDS(ggplot_hist_simple_repr_induct_short,"RDS/ggplot_hist_simple_repr_induct_short.rds")

saveRDS(process_TP_table,"RDS/process_TP_table.rds")
saveRDS(process_TP_table_REACT,"RDS/process_TP_table_REACT.rds")

saveRDS(groups_table,"RDS/groups_table.rds")

saveRDS(gpr_list_shrunken_mean,"RDS/gpr_list_shrunken_mean.rds")

saveRDS(slopes_frame_TC,"RDS/slopes_frame_TC.rds")
saveRDS(second_derivate,"RDS/second_derivate.rds")

saveRDS(ggplot_hist_repr_induct_complete,"RDS/ggplot_hist_repr_induct_complete.rds")

saveRDS(ztrans_log2FCs, "RDS/ztrans_log2FCs_TC.rds")
saveRDS(ztrans_log2FCs_all, "RDS/ztrans_log2FCs_all_TC.rds")

