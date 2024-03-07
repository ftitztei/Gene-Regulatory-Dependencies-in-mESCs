## This script aims to build a network based on possible gene regulatory
## dependencies that were consistent between the SC data and the KO data.

library(biomaRt)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(foreach)
library(GeneOverlap)
library(topGO)
library(GO.db)
library(Mus.musculus)


## load functions for this script and other scripts in the analysis
source("R/functions.R")
## read in objects need for this script
A_relation_b_hm_consist_KO_SC <- readRDS("RDS/A_relation_b_hm_consist_KO_SC_genes_Christa_sct_counts_600.rds")
A_relation_b_hm_consist_KO_SC <- A_relation_b_hm_consist_KO_SC[which(rowSums(abs(A_relation_b_hm_consist_KO_SC))!=0),
                                                               which(colSums(abs(A_relation_b_hm_consist_KO_SC))!=0)]

A_relation_b_hm_consist_KO_SC_mgi <- A_relation_b_hm_consist_KO_SC
row.names(A_relation_b_hm_consist_KO_SC_mgi) <- SC_lookup[match(row.names(A_relation_b_hm_consist_KO_SC),SC_lookup$ensembl),2]
colnames(A_relation_b_hm_consist_KO_SC_mgi) <- SC_lookup[match(row.names(A_relation_b_hm_consist_KO_SC),SC_lookup$ensembl),2]

CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients
SC_lookup <- readRDS("RDS/count_sct_Christa_lookup.rds")
## build lookup table for mgi symbols and ensembl gene names
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id", "mgi_symbol",
                                "description"),
                 mart=ensembl)
biomart <- biomart[which(!is.na(biomart$entrezgene)),]
biomart <- biomart[complete.cases(biomart[,c(1,3)]),]

gpr_list <- readRDS("RDS/gpr_list_shrunkenFCS.rds")
log_norm_heatmap_gpr <- gpr_list[[1]]
log_norm_heatmap_gpr <- log_norm_heatmap_gpr-log_norm_heatmap_gpr[,2]
log_norm_heatmap_gpr_changed <- log_norm_heatmap_gpr[which(rowMax(log_norm_heatmap_gpr[,2:18])-
                                                             rowMin(log_norm_heatmap_gpr[,2:18])>=0.5),]
changed_genes_ensembl <- row.names(log_norm_heatmap_gpr_changed)
changed_genes_mgi <- biomart[which(biomart$ensembl_gene_id%in%changed_genes_ensembl),3]

## test connectivity of clusters which clusters contain genes that could be in other clusters
NMF_list <- readRDS("RDS/NMF_list.rds")

H_norm <- NMF_list[["H_nmf"]]
H_norm <- apply(H_norm, 2, function(x) x/max(x))

clusters_h <- readRDS("RDS/clusters_h.rds")

nmf_rank <- 6

connect_mat <- matrix(NA, nmf_rank, nmf_rank)
row.names(connect_mat) <- paste("cluster",1:nmf_rank,sep = "_")
colnames(connect_mat) <- paste("cluster",1:nmf_rank,sep = "_")
for(r in 1:nmf_rank){
  for(c in 1:nmf_rank){
    connect_mat[r,c] <- round(mean(H_norm[c,which(clusters_h[,1]==r)]),digits = 3)
  }
}
## Fig 5.30 C
Heatmap(connect_mat,
        col=colorRamp2(c(1,0), c("red",'white')),
        show_row_dend = T, show_column_dend = T,
        show_row_names = T, show_column_names = T,
        cluster_rows = hclust(dist(connect_mat),method="ward.D2"),
        cluster_columns = hclust(dist(connect_mat),method="ward.D2"),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 8),
        row_gap = unit(1,"mm"),
        column_title_rot = 0,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(1,"mm"),
        name = "mean cluster identity",
        column_title = paste("genes from which clusters \n have high identity with other clusters \n",
                             nmf_rank, " NMF clusters\n",sep=""))



### create heatmap for bigger clusters individually ####

clusters_mgi <- data.frame(mgi=NMF_list[[1]][which(NMF_list[[1]]$gene%in%
                                                     row.names(A_relation_b_hm_consist_KO_SC_mgi)),1],
                           ensembl=SC_lookup[match(NMF_list[[1]][which(NMF_list[[1]]$gene%in%
                                                                         row.names(A_relation_b_hm_consist_KO_SC_mgi)),1],SC_lookup$mgi),1],
                           cluster=NMF_list[[1]][which(NMF_list[[1]]$gene%in%
                                                         row.names(A_relation_b_hm_consist_KO_SC_mgi)),2])

### how do formative and naive markers behave ####

A_relation_b_naive_formative <- A_relation_b_hm_consist_KO_SC_mgi[which(row.names(A_relation_b_hm_consist_KO_SC_mgi)%in%
                                                                          c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                                                                            "Klf4", "Prdm14","Fgf5", "Pou3f1", "Otx2",
                                                                            "Dnmt3a", "Dnmt3b")),]
A_relation_b_naive_formative <- A_relation_b_naive_formative[,which(apply(A_relation_b_naive_formative,2, function(x) sum(x==0))!=8)]

group <- rep(NA, nrow(A_relation_b_naive_formative))
group[which(row.names(A_relation_b_naive_formative)%in%
              c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                "Klf4", "Prdm14"))]  <- "naive markers"
group[which(row.names(A_relation_b_naive_formative)%in%
              c("Fgf5", "Pou3f1", "Otx2",
                "Dnmt3a", "Dnmt3b"))]  <- "formative markers"

ha <- rowAnnotation(group=group,
                    col=list(group=c("naive markers"="#57BC77","formative markers"="#719AC4")))

Heatmap(A_relation_b_naive_formative,
        col=colorRamp2(c(-1.2,0,1.2), c("blue",'white',"red")),
        show_row_dend = F, show_column_dend = F,
        show_row_names = F, show_column_names = F,
        cluster_rows = F,#hclust_relations,
        cluster_columns = F,#hclust_relations, row_title_rot = 0,
        row_title_gp = gpar(fontsize = 8),
        row_gap = unit(1,"mm"), column_title_rot = 0,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(1,"mm"),
        name = "mean activity diff A to B", left_annotation = ha, #right_annotation = ha,
        column_title = "consistent dependencies between A (row) and B (column) \n naive and formative markers in rows\n")

### have a look at average relation between clusters ####
clusters_sorted <- sort(unique(clusters_mgi$cluster))#[c(1:3,6:25,4:5)]
clusters_hclust_relation <- data.frame(cluster = clusters_mgi$cluster,
                                       row.names = clusters_mgi$ensembl)
clusters_hclust_relation[,2:(length(clusters_sorted)+1)] <- NA
colnames(clusters_hclust_relation) <- c("cluster",paste("mean_relation_to_cluster_",clusters_sorted,sep=""))

## for each cluster get mean of realtion of all process process pairs forthe two clusters
for(d in 2:(length(clusters_sorted)+1)){
  print(d-1)
  for(t in 1:dim(clusters_hclust_relation)[1]){
    clusters_hclust_relation[t,d] <- mean(as.numeric(A_relation_b_hm_consist_KO_SC[t,which(clusters_hclust_relation$cluster==clusters_sorted[d-1])]))
  }
}


cluster_hclust_relation_summary <- data.frame(number_of_processes=rep(NA,length(clusters_sorted)),
                                              row.names=paste("cluster_",clusters_sorted,sep=""))
cluster_hclust_relation_summary[,2:(length(clusters_sorted)+1)] <- NA
colnames(cluster_hclust_relation_summary) <- c("number_of_genes",
                                               paste("mean_relation_to_cluster_",clusters_sorted,sep=""))
cluster_hclust_relation_summary$genes <- NA

for(d in 1:nrow(cluster_hclust_relation_summary)){
  cluster_hclust_relation_summary[d,1] <- length(row.names(clusters_hclust_relation)[which(
    clusters_hclust_relation[,1]==clusters_sorted[d])])
  cluster_hclust_relation_summary$genes[d] <- paste(row.names(clusters_hclust_relation)[which(
    clusters_hclust_relation[,1]==clusters_sorted[d])],collapse=", ")
}


for(d in 1:dim(cluster_hclust_relation_summary)[1]){
  for(c in 1:dim(cluster_hclust_relation_summary)[1]){
    if(length(which(clusters_hclust_relation$cluster==clusters_sorted[d]))==1|
       length(which(clusters_hclust_relation$cluster==clusters_sorted[c]))==1){
      cluster_hclust_relation_summary[d,c+1] <- mean(as.numeric(A_relation_b_hm_consist_KO_SC[which(clusters_hclust_relation$cluster==clusters_sorted[d]),
                                                                     which(clusters_hclust_relation$cluster==clusters_sorted[c])]))
    }
    else{
      cluster_hclust_relation_summary[d,c+1] <- mean(rowMeans(A_relation_b_hm_consist_KO_SC[which(clusters_hclust_relation$cluster==clusters_sorted[d]),
                                                                              which(clusters_hclust_relation$cluster==clusters_sorted[c])]))
    }
  }
}


ha5 <-rowAnnotation(group=cluster_hclust_relation_summary$number_of_genes)

## Fig 5.30 A
Heatmap(cluster_hclust_relation_summary[,2:(length(clusters_sorted)+1)],
        col=colorRamp2(c(-0.30,0,0.30), c("blue",'white',"red")),
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = T,
        cluster_rows = T, cluster_columns = T,
        name = "mean activity diff A to B",
        row_labels = paste("cluster",clusters_sorted),
        column_labels = paste("cluster",clusters_sorted),
        left_annotation = ha5,
        column_title = paste("relation between clusters \n",
                             "\n",sep=""))


## remove number of processes from table to only have matrix left
A_relation_b_hm_consist_clust <- cluster_hclust_relation_summary[,2:(length(clusters_sorted)+1)]


### get cytoscape output ####
cluster_hclust_relation_short <- cluster_hclust_relation_summary[,2:(length(clusters_sorted)+1)]

for(r in 1:nrow(cluster_hclust_relation_short)){
  for(c in 1:ncol(cluster_hclust_relation_short)){
    if(abs(cluster_hclust_relation_short[r,c])<0.01){
      cluster_hclust_relation_short[r,c] <- 0
    }
  }
}


cytoframe <- data.frame(out_node=character(),
                        in_node=character(),
                        directed=logical(),
                        mean_relation=double())



for(r in 1:dim(cluster_hclust_relation_short)[1]){
  for(c in 1:dim(cluster_hclust_relation_short)[2]){
## only look at non redundant relations
    if(c >= r){
## only add non zero relations
      if(cluster_hclust_relation_short[r,c]!=0){
## get direction if positive row is stronger delayed than column
## therefore arrow from column to row as row seems to need column
        if(cluster_hclust_relation_short[r,c]<0){
          cytoframe <- rbind(cytoframe,
                             data.frame(out_node=paste("cluster",clusters_sorted[c],sep="_"),
                                        in_node=paste("cluster",clusters_sorted[r],sep="_"),
                                        directed=TRUE,
                                        mean_relation=abs(cluster_hclust_relation_short[r,c])))
        }
        if(cluster_hclust_relation_short[r,c]>0){
          cytoframe <- rbind(cytoframe,
                             data.frame(out_node=paste("cluster",clusters_sorted[r],sep="_"),
                                        in_node=paste("cluster",clusters_sorted[c],sep="_"),
                                        directed=TRUE,
                                        mean_relation=abs(cluster_hclust_relation_short[r,c])))
        }
      }
    }
  }
}
## which edges can be removed??
## egdes that connect two nodes that are connected through
## intermediated edges and nodes can be removed from the graph

## get essential edges for all nodes that do not have an incoming edge
cyto_cluster_2 <- get_essential_edges(cytoframe,"cluster_2")
# cyto_cluster_4 <- get_essential_edges(cytoframe,"cluster_4")
# cyto_cluster_5 <- get_essential_edges(cytoframe,"cluster_5")
## merge all graphs based on no in edge node graphs
cyto_union <- cyto_cluster_2 #unique(rbind(cyto_cluster_2,cyto_cluster_4,cyto_cluster_5))

graph_frame <- cyto_union[,c(1,2,4)]
colnames(graph_frame) <- c("out_node","in_node","weight")


net <- graph_from_data_frame(d=graph_frame,
                             vertices = paste("cluster",clusters_sorted,sep = "_"))
## not do it manually anymore but has do be done once to get coordinates ##
tkid <- tkplot(net)
layout <- tkplot.getcoords(tkid)
#layout <- matrix(data=c(200, 111,
#                        322, 216,
#                        265, 38,
#                        366, 151,
#                        261, 154,
#                        197, 0),nrow=6,ncol=2,byrow = T)

plot(net, layout = layout,
     edge.width = E(net)$weight*100,
     edge.arrow.size=1.5)

### add biological information to genes ####
clusters_mgi$is_TF <- 0
clusters_mgi$is_NAG <- 0
clusters_mgi$is_upNAG <- 0
clusters_mgi$is_downNAG <- 0
#clusters_mgi$gene_function <- NA

CLIM2 <- readRDS("/data/public/ftitztei/Rscripts/TCpackage/data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients

NAGs_frame <- read.table("/data/public/ftitztei/100KO/updated_results_Kdm6a/basic_tables/results_marker_multiple_regression_updated.tsv",
                         sep = "\t",colClasses = c("character","numeric","numeric",
                                                   "numeric","numeric","numeric",
                                                   "numeric","numeric","numeric",
                                                   "numeric","numeric"),
                         skip = 1)
colnames(NAGs_frame) <- c("gene", "Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
                          "Klf4", "Prdm14", "Zfp42", "rsq", "nondet", "nondet_z")
NAGs <- NAGs_frame[which(NAGs_frame$rsq >= 0.65),1]

upNAGs <- NAGs[which(NAGs%in%row.names(KOvKO[which(KOvKO[,74]>0),]))]

downNAGs <- NAGs[which(NAGs%in%row.names(KOvKO[which(KOvKO[,74]<0),]))]


clusters_mgi[which(clusters_mgi$mgi%in%TFs),4] <- 1
clusters_mgi[which(clusters_mgi$mgi%in%NAGsmmu04630),5] <- 1
clusters_mgi[which(clusters_mgi$mgi%in%upNAGs),6] <- 1
clusters_mgi[which(clusters_mgi$mgi%in%downNAGs),7] <- 1
write.csv(clusters_mgi, "output_tables/cluster_genes_NAGs_TFs_280623.csv",
          quote = F,row.names = F)


### plots ####
A_relation_b_hm_consist <- readRDS("RDS/A_relation_b_hm_consist_perc_genes.rds")

GO_mapping_entrez <- as.list(org.Mm.egGO2ALLEGS)
Reactome_mapping_entrez <- loadGeneSet(organism="mmusculus",enrichDatabase = "pathway_Reactome")
KEGG_mapping_entrez <- loadGeneSet(organism="mmusculus",enrichDatabase = "pathway_KEGG")
WikiPW_mapping_entrez <- loadGeneSet(organism="mmusculus",enrichDatabase = "pathway_Wikipathway")

naive <- c("Nanog", "Tfcp2l1", "Esrrb", "Tbx3",
           "Klf4", "Prdm14")

formative <- c("Fgf5", "Pou3f1", "Otx2",
               "Dnmt3a", "Dnmt3b")



## all no subset
## Fig 5.29 B
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist_KO_SC_mgi,
                       group_names=c("naive markers","formative markers"),
                       group_genes=list(A=naive,
                                        B=formative),
                       group_colors=c("#57BC77","#719AC4"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist_KO_SC_mgi)),],
                       title_top=paste("all consistent dependencies from SC and KO analysis \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = F,
                       show_names=F)
## FIg 5.31 B
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist_KO_SC_mgi,
                       group_names=c("naive markers","formative markers"),
                       group_genes=list(A=naive,
                                        B=formative),
                       group_colors=c("#57BC77","#719AC4"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist_KO_SC_mgi)),],
                       title_top=paste("all consistent dependencies from SC and KO analysis \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = T,
                       show_names=T)
## Fig 5.16
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist,
                       group_names=c("naive markers","formative markers"),
                       group_genes=list(A=naive,
                                        B=formative),
                       group_colors=c("#57BC77","#719AC4"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist)),],
                       title_top=paste("consistent dependencies from KO analysis only \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = F,
                       show_names=F)

A_relation_b_hm_nocontra <- readRDS("RDS/A_relation_b_hm_nocontra.rds")
## Fig 5.29 A
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_nocontra,
                       group_names=c("naive markers","formative markers"),
                       group_genes=list(A=naive,
                                        B=formative),
                       group_colors=c("#57BC77","#719AC4"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_nocontra)),],
                       title_top=paste("consistent dependencies from KO analysis filtered by SC \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = F,
                       show_names=F)
## Fig 5.31 A
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_nocontra,
                       group_names=c("naive markers","formative markers"),
                       group_genes=list(A=naive,
                                        B=formative),
                       group_colors=c("#57BC77","#719AC4"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_nocontra)),],
                       title_top=paste("consistent dependencies from KO analysis filtered by SC \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = T,
                       show_names=T)


### Boroviak
genes_boroviak_ICM_preimpl <- c("Nanog", "Bmp4","Zfp42",
                                "Gdf3","Prdm14", "Egr1",
                                "Tdgf1", "Fgf3","Mybl2",
                                "Il6st","Klf6","Lefty1",
                                "Zfp57")

genes_boroviak_preimpl <- c("Dnmt3l", "Lifr","Etv4",
                            "Nodal","Nr0b1","Pdgfa",
                            "Zfp36l1","Apobec2")

genes_boroviak_preimpl_postimpl <- c("Sox2","Tcf7l1","Utf1",
                                     "Foxd3","Zic3","Zscan10",
                                     "Etv5","Hif1a","Fgf4",
                                     "Fgf15","Map2k5","Lgr4",
                                     "Jub","Pcgf2")

genes_boroviak_postimpl <- c("Pou3f1", "Otx2", "Sox3",
                             "Sall2","Tead2","Pim2",
                             "Bex1","Hes6","Fgf5",
                             "Mapk12","Grb10","Fzd2",
                             "Fzd7","Rspo1","Smo",
                             "Notch3","Dnmt3b","Dnmt3",
                             "Hmga1")


## Fig 5.20 B
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist,
                       group_names=c("ICM preimplantation EPI",
                                     "preimplantation EPI",
                                     "preimplantation postimplantation EPI",
                                     "postimplantation EPI"),
                       group_genes=list(A=genes_boroviak_ICM_preimpl,
                                        B=genes_boroviak_preimpl,
                                        C=genes_boroviak_preimpl_postimpl,
                                        D=genes_boroviak_postimpl),
                       group_colors=c("#57BC77","#336941","#719AC4","#4A7BAF"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist)),],
                       title_top=paste("consistent dependencies from KO analysis only \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = T,
                       show_names=T)

plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist_KO_SC_mgi,
                       group_names=c("ICM preimplantation EPI",
                                     "preimplantation EPI",
                                     "preimplantation postimplantation EPI",
                                     "postimplantation EPI"),
                       group_genes=list(A=genes_boroviak_ICM_preimpl,
                                        B=genes_boroviak_preimpl,
                                        C=genes_boroviak_preimpl_postimpl,
                                        D=genes_boroviak_postimpl),
                       group_colors=c("#57BC77","#336941","#719AC4","#4A7BAF"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist_KO_SC_mgi)),],
                       title_top=paste("consistent dependencies from KO analysis only \n ",
                                       nmf_rank, " NMF clusters\n",sep=""),
                       plot_out=TRUE,subset_genes = T,
                       show_names=T)


### Signaling pathways regulating pluripotency
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist_KO_SC_mgi,
                       group_names=c("naive markers","formative markers",
                                     "Signaling pathways regulating pluripotency"),
                       group_genes=list(A=naive,
                                        B=formative,
                                        C=biomart[which(biomart$entrezgene_id%in%
                                                          KEGG_mapping_entrez[[1]][which(KEGG_mapping_entrez[[1]][,1]=="mmu04550"),3]),3]),
                       group_colors=c("#57BC77","#719AC4","orange"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist_KO_SC_mgi)),],
                       title_top=paste("subset consistent dependencies from SC and KO analysis \n ",
                                       "Signaling pathways regulating pluripotency\n",sep=""),
                       plot_out=TRUE,subset_genes = T,
                       show_names=T)

## Fig 5.20 A
plotDependenciesGroups(heatmap_matrix=A_relation_b_hm_consist,
                       group_names=c("naive markers","formative markers",
                                     "Signaling pathways regulating pluripotency"),
                       group_genes=list(A=naive,
                                        B=formative,
                                        C=biomart[which(biomart$entrezgene_id%in%
                                                          KEGG_mapping_entrez[[1]][which(KEGG_mapping_entrez[[1]][,1]=="mmu04550"),3]),3]),
                       group_colors=c("#57BC77","#719AC4","orange"),
                       clusters=clusters_h[which(clusters_h$gene%in%
                                                   row.names(A_relation_b_hm_consist)),],
                       title_top=paste("subset consistent dependencies from KO analysis only \n ",
                                       "Signaling pathways regulating pluripotency\n",sep=""),
                       plot_out=TRUE,subset_genes = T,
                       show_names=T)
