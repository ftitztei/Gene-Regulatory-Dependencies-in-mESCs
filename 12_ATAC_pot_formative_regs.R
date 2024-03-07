library(GenomicFeatures)
library(Biostrings)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ChIPseeker)
library(ggplot2)
library(TFTargetCaller)
library(biomaRt)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(universalmotif)
library(foreach)
library(RColorBrewer)
library(bezier)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggpmisc)
source("R/sankey.R")



txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
Mmusculus <- BSgenome.Mmusculus.UCSC.mm10

### read in fimo files and annotate ####
fimo2GRanges <- function(fimo_frame){
  GRanges_frame <- data.frame(seqnames=sapply(fimo_frame$sequence_name, FUN= function(x) strsplit(x,"_")[[1]][1]),
                              start=as.numeric(sapply(fimo_frame$sequence_name, FUN= function(x) strsplit(x,"_")[[1]][2]))+fimo_frame$start,
                              end=as.numeric(sapply(fimo_frame$sequence_name, FUN= function(x) strsplit(x,"_")[[1]][2]))+fimo_frame$stop,
                              width=as.numeric(sapply(fimo_frame$sequence_name, FUN= function(x) strsplit(x,"_")[[1]][2]))+fimo_frame$stop-
                                as.numeric(sapply(fimo_frame$sequence_name, FUN= function(x) strsplit(x,"_")[[1]][2]))+fimo_frame$start+1,
                              strand="*",
                              TF=fimo_frame$TF,
                              Motif=fimo_frame$Motif,
                              q.value=fimo_frame$q.value)
  GRanges_frame <- makeGRangesFromDataFrame(GRanges_frame,keep.extra.columns = T)
}

TF_table <- read.table("/cellfile/datapublic/ftitztei/Luis_ATAC/CisBP/TF_Information.txt",
                       sep = "\t", header = T)

fimo_files <- c("/cellfile/datapublic/ftitztei/Luis_ATAC/FIMO/fimo_naive_formative_all_3/fimo.tsv")

fimo_list <- list()
GRanges_fimo_list <- list ()
Anno_fimo_list <- list()

for(f in fimo_files){
  fimo_holder <- read.table(f,
                            header =T, colClasses = c("character",
                                                      "character",
                                                      "character",
                                                      "numeric",
                                                      "numeric",
                                                      "character",
                                                      "numeric",
                                                      "numeric",
                                                      "numeric",
                                                      "character"))

  #fimo_holder <- fimo_holder[which(fimo_holder$q.value<=0.01),]


  fimo_holder$TF <- sapply(fimo_holder$motif_id, FUN= function(x) strsplit(x,"_")[[1]][1])
  fimo_holder$Motif <- sapply(fimo_holder$motif_id, FUN= function(x) strsplit(x,"_")[[1]][2])

  fimo_list <- c(fimo_list, list(fimo_holder))

  GRanges_holder <- fimo2GRanges(fimo_holder)
  GRanges_fimo_list <- c(GRanges_fimo_list, list(GRanges_holder))
  Anno_fimo_list <- c(Anno_fimo_list, list(annotatePeak(GRanges_holder, TxDb=txdb)))
}

names(fimo_list) <- c("naive_formative_all_3")

names(GRanges_fimo_list) <- c("naive_formative_all_3")

names(Anno_fimo_list) <- c("naive_formative_all_3")
Anno_df_fimo_list <- lapply(Anno_fimo_list,as.data.frame)

#length(unique(fimo_naive$TF))
#length(unique(fimo_naive$Motif))
#length(unique(fimo_naive$sequence_name))


gg_TSS_fimo <- data.frame(distanceToTSS=c(Anno_df_fimo_list$naive_formative_all_3$distanceToTSS),
                          q.value=c(Anno_df_fimo_list$naive_formative_all_3$q.value),
                          group=c(rep("FIMO binding naive formative all (3)",
                                      length(Anno_df_fimo_list$naive_formative_all_3$q.value))))

ggplot(gg_TSS_fimo, aes(x=log10(abs(distanceToTSS)+1), color = group)) +
  geom_density()+
  ggtitle("distance to next TSS")+
  theme(plot.title = element_text(hjust = 0.5))

plotAnnoPie(Anno_fimo_list[["naive_formative_all_3"]],
            main="\n \n \n naive formative peaks\noverlap 3")

### TFTargetCaller ####
TFTcall_split_TFs <- function(fimo_frame, positionFrame,
                              exclude=NA, method="ClosestGene",
                              size=10000){
  TFs <- unique(fimo_frame$TF)
  holder <- list()
  for(d in 1:length(TFs)){
    if(!is.na(exclude)){
      if(d %in% exclude){

      }
      else{
        print(paste("processing", as.character(d), "of", length(TFs), "TFs", collapse = " "))
        fimo_frame_TF <- fimo_frame[which(fimo_frame$TF==TFs[d]),]
        peaks_position <- data.frame("chromosome"=sapply(fimo_frame_TF$sequence_name,
                                                         FUN= function(x) strsplit(x,"_")[[1]][1]),
                                     "center"=(apply(fimo_frame_TF, 1, FUN= function(x) {
                                       as.numeric(strsplit(as.character(x[3]),"_")[[1]][2])+as.numeric(x[4])})+
                                         apply(fimo_frame_TF, 1, FUN= function(x) {
                                           as.numeric(strsplit(as.character(x[3]),"_")[[1]][2])+as.numeric(x[5])}))/2)

        if(method == "ClosestGene"){
          holder <- c(holder, list(TFTargetCaller(peaks_position, positionFrame,
                                                  method = method, ClosestGeneScore = "qvalue")))
        }else{
          holder <- c(holder, list(TFTargetCaller(peaks_position, positionFrame,
                                                  method = method, n = size)))
        }

      }
    }
    else{
      print(paste("processing", as.character(d), "of", length(TFs), "TFs", collapse = " "))
      fimo_frame_TF <- fimo_frame[which(fimo_frame$TF==TFs[d]),]
      peaks_position <- data.frame("chromosome"=sapply(fimo_frame_TF$sequence_name,
                                                       FUN= function(x) strsplit(x,"_")[[1]][1]),
                                   "center"=(apply(fimo_frame_TF, 1, FUN= function(x) {
                                     as.numeric(strsplit(as.character(x[3]),"_")[[1]][2])+as.numeric(x[4])})+
                                       apply(fimo_frame_TF, 1, FUN= function(x) {
                                         as.numeric(strsplit(as.character(x[3]),"_")[[1]][2])+as.numeric(x[5])}))/2)

      if(method == "ClosestGene"){
        holder <- c(holder, list(TFTargetCaller(peaks_position, positionFrame,
                                                method = method, ClosestGeneScore = "qvalue")))
      }else{
        holder <- c(holder, list(TFTargetCaller(peaks_position, positionFrame,
                                                method = method, n = size)))
      }
    }



  }
  if(!is.na(exclude)){
    names(holder) <- TFs[which(!1:length(TFs)%in%exclude)]
  }
  else{
    names(holder) <- TFs
  }

  return(holder)
}
target_table_from_TFTres <- function(TFT_res, mart, cutoff = 0.05, binary=FALSE){
  holder <- list()
  for(d in 1:length(TFT_res)){
    if(binary==FALSE){
      if(length(names(TFT_res[[d]])[which(TFT_res[[d]]<cutoff)]>0)){

        target_table <- data.frame(ensembl=names(TFT_res[[d]])[which(TFT_res[[d]]<cutoff)],
                                 mgi=NA,
                                 qval=TFT_res[[d]][which(TFT_res[[d]]<cutoff)])

        for(r in 1:dim(target_table)[1]){

          if(length(mart[which(mart$ensembl_gene_id==target_table[r,1]),3])==1){

            target_table[r,2] <- mart[which(mart$ensembl_gene_id==target_table[r,1]),3]
            }
          }
        holder <- c(holder, list(target_table))
      }
      else{
        holder <- c(holder,list(NA))
      }
      }else{
      if(length(names(TFT_res[[d]])[which(TFT_res[[d]]==1)]>0)){

        target_table <- data.frame(ensembl=names(TFT_res[[d]])[which(TFT_res[[d]]==1)],
                                   mgi=NA,
                                   qval=TFT_res[[d]][which(TFT_res[[d]]==1)])

        for(r in 1:dim(target_table)[1]){

          if(length(mart[which(mart$ensembl_gene_id==target_table[r,1]),3])==1){

            target_table[r,2] <- mart[which(mart$ensembl_gene_id==target_table[r,1]),3]
          }
        }

        holder <- c(holder, list(target_table))
      }
        else{
          holder <- c(holder,list(NA))
        }
    }




  }
  names(holder) <- names(TFT_res)


  return(holder)
}

ensembl102 <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",
                         version = "102")
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol",
                                "chromosome_name","start_position","end_position"),
                 mart=ensembl)
biomart102 <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol",
                                "chromosome_name","start_position","end_position"),
                 mart=ensembl102)


files <- list.files("/cellfile/datapublic/ftitztei/Luis_ATAC/ATACseq_peaks/final_bed/")
files <- files[13]
TFTfiles <- files

data(genePosition)
gene_position_102 <- getBM(attributes = c("ensembl_gene_id","chromosome_name",
                                          "strand","transcription_start_site","gene_biotype"),
                           mart=ensembl102)
gene_position_102 <- gene_position_102[which(gene_position_102$gene_biotype %in%
                                               c("protein_coding", "lincRNA", "miRNA")),1:4]
gene_position_102$keep <- 0
for(gene in unique(gene_position_102$ensembl_gene_id)){
  rows_gene <- which(gene_position_102$ensembl_gene_id==gene)
  min_row <- which(gene_position_102[rows_gene, 4] == min(gene_position_102[rows_gene, 4]))
  gene_position_102[rows_gene[min_row[1]],5] <- 1
}
gene_position_102 <- gene_position_102[which(gene_position_102$keep==1),1:4]

gene_position_102[,2] <- paste0("chr",gene_position_102[,2])
colnames(gene_position_102) <- colnames(genePosition)
row.names(gene_position_102) <- gene_position_102$geneID

# saveRDS(ensembl102, "RDS/ensembl102.rds")
# saveRDS(biomart102, "RDS/biomart102.rds")
# saveRDS(gene_position_102, "RDS/gene_position_102.rds")

peaks_position_list <- list()

for(f in TFTfiles){
  peak_holder <- getPeakPosition(paste("/cellfile/datapublic/ftitztei/Luis_ATAC/ATACseq_peaks/final_bed/", f, sep =""))
  peaks_position_list <- c(peaks_position_list, list(peak_holder))
}

names(peaks_position_list) <- c("naive_formative_all_3")
TFT_list <- list()
TFT_target_tables_list <- list()

TFT_names <- c("naive_formative_all_3")

for(n in TFT_names){
  TFT_call_holder <- TFTcall_split_TFs(fimo_list[[n]], gene_position_102)
  TFT_list <- c(TFT_list, list(TFT_call_holder))

  # TFT_call_holder <- TFT_list[[TFT_names]]

  TFT_table_holder <- target_table_from_TFTres(TFT_call_holder, biomart102, cutoff = 0.05)
  TFT_table_holder <- TFT_table_holder[which(!is.na(TFT_table_holder))]
  TFT_target_tables_list <- c(TFT_target_tables_list, list(TFT_table_holder))

}

names(TFT_list) <- c("naive_formative_all_3")

names(TFT_target_tables_list) <- c("naive_formative_all_3")
TFT_target_tables_list <- TFT_target_tables_list[[1]]



extended_form_KO <- readRDS("RDS/extended_form_KO.rds")
rev_form_KO <- readRDS("RDS/rev_form_KO.rds")



TFT_target_tables_list_ext_form_holder <- TFT_target_tables_list
TFT_target_tables_list_rev_form_holder <- TFT_target_tables_list

for(TF in 1:length(TFT_target_tables_list_ext_form_holder)){
  TFT_target_tables_list_ext_form_holder[[TF]] <- TFT_target_tables_list_ext_form_holder[[TF]][which(TFT_target_tables_list_ext_form_holder[[TF]][,2]
                                                                       %in%extended_form_KO),]

  TFT_target_tables_list_rev_form_holder[[TF]] <- TFT_target_tables_list_rev_form_holder[[TF]][which(TFT_target_tables_list_rev_form_holder[[TF]][,2]
                                                                                                     %in%rev_form_KO),]
}

TFT_target_tables_list_ext_form <- TFT_target_tables_list_ext_form_holder[which(lapply(TFT_target_tables_list_ext_form_holder, nrow)>0)]
TFT_target_tables_list_ext_form_5 <- TFT_target_tables_list_ext_form_holder[which(lapply(TFT_target_tables_list_ext_form_holder, nrow)>4)]


TFT_target_tables_list_rev_form <- TFT_target_tables_list_rev_form_holder[which(lapply(TFT_target_tables_list_rev_form_holder, nrow)>0)]
TFT_target_tables_list_rev_form_5 <- TFT_target_tables_list_rev_form_holder[which(lapply(TFT_target_tables_list_rev_form_holder, nrow)>4)]


gg_reg_ext_form <- data.frame(motif=rep(NA,sum(sapply(TFT_target_tables_list_ext_form_5,nrow))),
                              target=rep(NA,sum(sapply(TFT_target_tables_list_ext_form_5,nrow))),
                              TFT_qval=rep(NA,sum(sapply(TFT_target_tables_list_ext_form_5,nrow))))
index <- 0
for(l in 1:length(TFT_target_tables_list_ext_form_5)){
  gg_reg_ext_form[(index+1):(index+nrow(TFT_target_tables_list_ext_form_5[[l]])),1] <- names(TFT_target_tables_list_ext_form_5)[l]
  gg_reg_ext_form[(index+1):(index+nrow(TFT_target_tables_list_ext_form_5[[l]])),2:3] <- TFT_target_tables_list_ext_form_5[[l]][,2:3]

  index <- index+nrow(TFT_target_tables_list_ext_form_5[[l]])
}

gg_reg_rev_form <- data.frame(motif=rep(NA,sum(sapply(TFT_target_tables_list_rev_form_5,nrow))),
                              target=rep(NA,sum(sapply(TFT_target_tables_list_rev_form_5,nrow))),
                              TFT_qval=rep(NA,sum(sapply(TFT_target_tables_list_rev_form_5,nrow))))
index <- 0
for(l in 1:length(TFT_target_tables_list_rev_form_5)){
  gg_reg_rev_form[(index+1):(index+nrow(TFT_target_tables_list_rev_form_5[[l]])),1] <- names(TFT_target_tables_list_rev_form_5)[l]
  gg_reg_rev_form[(index+1):(index+nrow(TFT_target_tables_list_rev_form_5[[l]])),2:3] <- TFT_target_tables_list_rev_form_5[[l]][,2:3]

  index <- index+nrow(TFT_target_tables_list_rev_form_5[[l]])
}

TFT_target_tables_list_rev_form_5

gg_reg_form <- rbind(gg_reg_ext_form, gg_reg_rev_form)

binary_TF_target_form <- matrix(data = 0, nrow = length(unique(gg_reg_form$motif)),
                                    ncol = length(unique(gg_reg_form$target)))
colnames(binary_TF_target_form) <- unique(gg_reg_form$target)
row.names(binary_TF_target_form) <- unique(gg_reg_form$motif)

for(r in 1:nrow(gg_reg_form)){
  binary_TF_target_form[gg_reg_form[r,1], gg_reg_form[r,2]] <- 1
}

clust_binary_form_targets <- hclust(dist(t(binary_TF_target_form)), method="ward.D2" )
clust_binary_form_TFs <- hclust(dist(binary_TF_target_form), method="ward.D2" )
TF_split <- cutree(clust_binary_form_TFs, 7)

Heatmap(binary_TF_target_form,
        show_row_names = F,
        #cluster_rows = clust_binary_rev_form_TFs,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        cluster_columns = clust_binary_form_targets,
        row_split = TF_split,
        col=c("white","red"),
        name="target",
        border = TRUE)


### restrict TFs to those that are expressed in 2i ####
CRISP2 <- readRDS("data/RDS_data/CRISP2_Kdm6a_updated_091020.rds")

FPKM_mat_short <- data.frame(RC9_N2=rowMeans(CRISP2$fpkm[,which(colnames(CRISP2$fpkm)%in%
                                                                  row.names(CRISP2$design_rall[which(CRISP2$design_rall$gene_x2=="wtN1"),]))]),
                             RC9_2i=rowMeans(CRISP2$fpkm[,which(colnames(CRISP2$fpkm)%in%
                                                                  row.names(CRISP2$design_rall[which(CRISP2$design_rall$gene_x2=="wt2i"),]))]))

FPKM_mat_short_TFs <- FPKM_mat_short[which(row.names(FPKM_mat_short)%in%
                                             row.names(binary_TF_target_form)),]

expressed_TFs <- row.names(FPKM_mat_short_TFs)[which(FPKM_mat_short_TFs$RC9_2i>=1)]


TF_split_expr <- as.data.frame(TF_split)[expressed_TFs,]

expr_binary_TF_target <- binary_TF_target_form[expressed_TFs,]

clust_binary_expr_targets <- hclust(dist(t(expr_binary_TF_target)), method="ward.D2" )
clust_binary_expr_TFs <- hclust(dist(expr_binary_TF_target), method="ward.D2" )
TF_split_expr <- as.data.frame(cutree(clust_binary_expr_TFs, 7))

Heatmap(expr_binary_TF_target,
        show_row_names = F,
        #cluster_rows = clust_binary_rev_form_TFs,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        cluster_columns = clust_binary_expr_targets,
        row_split = TF_split_expr,
        col=c("white","red"),
        name="target",
        border = TRUE)

write.csv(expr_binary_TF_target, "output_tables/binary_TF_target_171123.csv",
          quote = F)

targets_TFs <- colnames(expr_binary_TF_target)[which(colnames(expr_binary_TF_target)%in%
                                                       unique(TF_table$TF_Name))]
target_TFs_binary_mat <- expr_binary_TF_target[,targets_TFs]
target_TFs_binary_mat <- target_TFs_binary_mat[which(rowSums(target_TFs_binary_mat)!=0),]

TF_clust_FPKM <- TF_split_expr
TF_clust_FPKM$FPKM <- FPKM_mat_short_TFs[which(FPKM_mat_short_TFs$RC9_2i>=1),2]

motifs_TFs_interest <- fimo_list[["naive_formative_all_3"]][which(fimo_list[["naive_formative_all_3"]]$TF %in%
                                                                    expressed_TFs),11:12]
motifs_TFs_interest <- motifs_TFs_interest[!duplicated(motifs_TFs_interest),]


motif_files <- sapply(unique(motifs_TFs_interest$Motif), function(x)
  paste0("/cellfile/datapublic/ftitztei/Luis_ATAC/CisBP/pwms_all_motifs/",x,"_2.00.txt"))

motifs_list <- lapply(motif_files, function(x) read_cisbp(x))
for(m in 1:length(motifs_list)){
  motifs_list[[m]]@name <- names(motifs_list)[m]
}

motif_compare <- compare_motifs(motifs_list, method = "EUCL", min.mean.ic = 0)

motif_hclust <- hclust(dist(motif_compare), method = "ward.D2")

gplots::heatmap.2(motif_compare,
                  trace = "none",
                  hclustfun = function(x) hclust(x, method = "ward.D2"))

gplots::heatmap.2(motif_compare,
                  trace = "none",
                  Rowv = as.dendrogram(motif_hclust),
                  Colv = as.dendrogram(motif_hclust))


## Roberts function to retrieve withtin cluster variance
withinClusterdist = function( m, clus) {
  require(foreach)
  withindist = foreach( c = unique(clus), .combine = c ) %do% {
    centroid = apply(m[,which(clus==c),drop=F],1,mean)
    return( centroid )
  }

  return( mean(withindist) )
}
## test cluster sizes in reagrads to within cluster variance
## on unscales data
WCV_cluster_size = foreach(i = 2:50, .combine=c) %do% {
  withinClusterdist( motif_compare,
                            cutree(motif_hclust,i) )
}
## plot change of within cluster variance to find good k
plot(WCV_cluster_size, type="b", xlim=c(0,50),
     xlab="number of clusters (k)", main="clusters of motifs",
     ylab="mean within cluster distance", pch=17, cex=1.5, col="black" )
abline(v=c(4,8), lwd=2, lty=2)


motif_cut_tree <- as.data.frame(cutree(motif_hclust, k = 8))
colnames(motif_cut_tree) <- "motif_cluster"
motif_cut_tree$TF <- NA
motif_cut_tree$TF_cluster <- NA

for(r in 1:nrow(motif_cut_tree)){
  TF_holder <- motifs_TFs_interest[which(motifs_TFs_interest$Motif==
                                           row.names(motif_cut_tree)[r]),1]
  motif_cut_tree[r,2] <- paste(TF_holder,
                               collapse = ";")
  motif_cut_tree[r,3] <- paste(TF_split_expr[TF_holder,1],
                               collapse = ";")
}


gplots::heatmap.2(motif_compare,
                  trace = "none",
                  hclustfun = function(x) hclust(x, method = "ward.D2"),
                  RowSideColors=brewer.pal(8, "Set3")[as.numeric(motif_cut_tree$motif_cluster)],
                  ColSideColors=brewer.pal(8, "Set3")[as.numeric(motif_cut_tree$motif_cluster)])
legend(x="topleft", legend=unique(as.numeric(motif_cut_tree$motif_cluster)),
       col=brewer.pal(8, "Set3")[unique(as.numeric(motif_cut_tree$motif_cluster))],
       pch=15)


ha = rowAnnotation(m1=as.numeric(row.names(expr_binary_TF_target)%in%
                              motif_cut_tree[which(motif_cut_tree$motif_cluster=="1"),2]),
                   m2=as.numeric(row.names(expr_binary_TF_target)%in%
                                   c(motif_cut_tree[which(motif_cut_tree$motif_cluster=="2"),2],
                                     "Tfdp1", "Tfdp2")),
                   m3=as.numeric(row.names(expr_binary_TF_target)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="3"),2]),
                   m4=as.numeric(row.names(expr_binary_TF_target)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="4"),2]),
                   m5=as.numeric(row.names(expr_binary_TF_target)%in%
                                   c(motif_cut_tree[which(motif_cut_tree$motif_cluster=="5"),2],
                                     "Mxi1", "Mxd1")),
                   m6=as.numeric(row.names(expr_binary_TF_target)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="6"),2]),
                   m7=as.numeric(row.names(expr_binary_TF_target)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="7"),2]),
                   m8=as.numeric(row.names(expr_binary_TF_target)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="8"),2]),
                   col=list(m1=c("1"="#8DD3C7", "0"="white"),
                             m2=c("1"="#FFFFB3", "0"="white"),
                             m3=c("1"="#BEBADA", "0"="white"),
                             m4=c("1"="#FB8072", "0"="white"),
                             m5=c("1"="#80B1D3", "0"="white"),
                             m6=c("1"="#FDB462", "0"="white"),
                             m7=c("1"="#B3DE69", "0"="white"),
                             m8=c("1"="#FCCDE5", "0"="white")))
## Fig 5.35 B
Heatmap(expr_binary_TF_target,
        show_row_names = F,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        cluster_columns = clust_binary_expr_targets,
        row_split = TF_split_expr,
        col=c("white","red"),
        name="target",
        border = TRUE,
        right_annotation = ha)

if(0){
### which TFs are important for differentiation ####
summary(as.factor(motif_cut_tree[which(motif_cut_tree$TF_cluster==3),1]))

summary(as.factor(motif_cut_tree[which(motif_cut_tree$TF_cluster==5),1]))

summary(as.factor(motif_cut_tree[which(motif_cut_tree$TF_cluster==6),1]))

## mostly 2 & 3 maybe cluster them a bit more
motif_compare_interest <- motif_compare[which(row.names(motif_compare)%in%
                                                c(row.names(motif_cut_tree[which(motif_cut_tree$TF_cluster%in%
                                                                                 c(3,5,6)),]),
                                                  row.names(target_TFs_binary_mat))),
                                        which(row.names(motif_compare)%in%
                                                c(row.names(motif_cut_tree[which(motif_cut_tree$TF_cluster%in%
                                                                                 c(3,5,6)),]),
                                                  row.names(target_TFs_binary_mat)))]

motif_hclust_interest <- hclust(dist(motif_compare_interest), method = "ward.D2")

motif_cut_tree_interest <- as.data.frame(cutree(motif_hclust_interest, k = 7))
colnames(motif_cut_tree_interest) <- "motif_cluster"
motif_cut_tree_interest$TF <- NA
motif_cut_tree_interest$TF_cluster <- NA

for(r in 1:nrow(motif_cut_tree_interest)){
  TF_holder <- motifs_TFs_interest[which(motifs_TFs_interest$Motif==
                                           row.names(motif_cut_tree_interest)[r]),1]
  motif_cut_tree_interest[r,2] <- paste(TF_holder,
                               collapse = ";")
  motif_cut_tree_interest[r,3] <- paste(TF_split_expr[TF_holder,1],
                               collapse = ";")
}
## Fig 5.35 A
gplots::heatmap.2(motif_compare_interest,
                  trace = "none",
                  hclustfun = function(x) hclust(x, method = "ward.D2"),
                  RowSideColors=brewer.pal(7, "Set3")[as.numeric(motif_cut_tree_interest$motif_cluster)],
                  ColSideColors=brewer.pal(7, "Set3")[as.numeric(motif_cut_tree_interest$motif_cluster)])
legend(x="topleft", legend=unique(as.numeric(motif_cut_tree_interest$motif_cluster)),
       col=brewer.pal(7, "Set3")[unique(as.numeric(motif_cut_tree_interest$motif_cluster))],
       pch=15)

summary(as.factor(motif_cut_tree_interest[which(motif_cut_tree_interest$TF_cluster==3),1]))

summary(as.factor(motif_cut_tree_interest[which(motif_cut_tree_interest$TF_cluster==5),1]))

summary(as.factor(motif_cut_tree_interest[which(motif_cut_tree_interest$TF_cluster==6),1]))

motif_cut_tree_interest$combi_TF_cluster_motif_cluster <- paste0(motif_cut_tree_interest$TF_cluster,
                                                                 ";",
                                                                 motif_cut_tree_interest$motif_cluster)

motif_cut_tree_interest_nonred <- motif_cut_tree_interest[!duplicated(motif_cut_tree_interest),]
motif_cut_tree_interest_nonred$FPKM <- NA
motif_cut_tree_interest_nonred$old_motif_clust <- NA

for(r in 1:nrow(motif_cut_tree_interest_nonred)){
  motif_cut_tree_interest_nonred$FPKM[r] <- TF_clust_FPKM[motif_cut_tree_interest_nonred[r,2],2]
  motif_cut_tree_interest_nonred$old_motif_clust[r] <- motif_cut_tree[row.names(motif_cut_tree_interest_nonred)[r],1]
}

motif_cut_tree_interest_nonred$scaled_clust_FPKM <- NA
for(comb in unique(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster)){
  rows <- which(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster==
                  comb)
  motif_cut_tree_interest_nonred[rows,7] <- motif_cut_tree_interest_nonred[rows,5]/
    max(motif_cut_tree_interest_nonred[rows,5])
}

sankey_frame <- motif_cut_tree_interest_nonred[,c(6,1)]
sankey_frame$flow <- NA
sankey_frame$from_val <- NA
sankey_frame$to_val <- NA
sankey_frame$col1 <- NA
sankey_frame$col2 <- NA


from_vals <- round(as.numeric(summary(as.factor(sankey_frame$old_motif_clust)))/nrow(sankey_frame)*100,2)
to_vals <- round(as.numeric(summary(as.factor(sankey_frame$motif_cluster)))/nrow(sankey_frame)*100,2)

for(r in 1:nrow(sankey_frame)){
  sankey_frame$flow[r] <- nrow(sankey_frame[which(sankey_frame$old_motif_clust==sankey_frame[r,1] &
                                                    sankey_frame$motif_cluster==sankey_frame[r,2]),])/
    nrow(sankey_frame[which(sankey_frame$old_motif_clust==sankey_frame[r,1]),])
  sankey_frame$from_val[r] <- from_vals[sankey_frame[r,1]]
  sankey_frame$to_val[r] <- to_vals[sankey_frame[r,2]]
  sankey_frame$col1[r] <- brewer.pal(9, "Set3")[sankey_frame[r,1]]
  sankey_frame$col2[r] <- brewer.pal(7, "Set3")[sankey_frame[r,2]]
}



p <- with(sankey_frame, ggSankeyGrad(c1 = old_motif_clust,
                               c2 = motif_cluster,
                               col1 = col1,
                               col2 = col2,
                               values = flow,
                               order1 = old_motif_clust,
                               order2 = motif_cluster,
                               #label2_col = from_val,
                               #label1_col = to_val,
                               alpha = 1,
                               label_color = F,
                               label_adjust = 2
                               ))
p
}
## sort by expression and motifs
## collect data
## RNA-seq TC log2FCs & TPMs
gpr_list <- readRDS("RDS/gpr_list_shrunkenFCS.rds")
tpms_TC_gpr <- readRDS("RDS/tpms_TC_gpr.rds")
tpm_table_TC <- readRDS("RDS/tpm_table.rds")

log_norm_heatmap_gpr <- gpr_list[[1]]
## set everything in relation to 2i
log_norm_heatmap_gpr <- log_norm_heatmap_gpr-log_norm_heatmap_gpr[,2]

## RNA-seq MC log2FCs
tpm_table_MC <- readRDS("RDS/tpms_MC.rds")
medium_change_shrunken <- readRDS("RDS/log2FCs_shrunken_MC.rds")
medium_change_results <- readRDS("RDS/MC_DESeq2_results.rds")

## prot TC log2FCs
prot_gpr <- read.table("/cellfile/datapublic/pkuzneco/master_thesis/data/processed/1_re-analyzed_data/proteomics_time_course/5_gaussian_process_regression/output_gp_FCs.tssv",
                      row.names = 1)
colnames(prot_gpr) <- prot_gpr[1,]
prot_gpr <- prot_gpr[2:nrow(prot_gpr),]
for(c in 1:ncol(prot_gpr)){
  prot_gpr[,c] <- as.numeric(prot_gpr[,c])
}

## prot MC log2FCs
prot_MC <- read.table("/cellfile/datapublic/pkuzneco/master_thesis/data/processed/1_re-analyzed_data/medium change/2_FC_shrinkage/1_on_lm_FCs/output_all_FCs.tsv",
                       row.names = 1)
colnames(prot_MC) <- prot_MC[1,]
prot_MC <- prot_MC[2:nrow(prot_MC),]
for(c in 1:ncol(prot_MC)){
  prot_MC[,c] <- as.numeric(prot_MC[,c])
}

## z-transformed TC log2FCs
ztrans_log2FCs <- readRDS("RDS/ztrans_log2FCs_TC.rds")
ztrans_log2FCs_all <- readRDS("RDS/ztrans_log2FCs_all_TC.rds")

## how do they behave in TC MC and Proteomics TC and MC
TF_overview_plot = function(tpms_gpr_rna,
                            tpms_unsmoothed_rna,
                            log2FCs_rna,
                            log2FCs_MC_rna,
                            log2FCs_TC_prot,
                            log2FCs_MC_prot,
                            ztrans_log2FCs,
                            Gene_info=FALSE,
                            Gene_info_frame,
                            Gene_info_name="scaled TPM: ",
                            plot_org_dots = TRUE,
                            targets = character(),
                            Genes=character(),
                            title="",
                            subtitle="",
                            numcol=3,
                            axis_ext=0.1) {

  ## make sure at least one gene provided was found in the list
  if(length(which(Genes%in%tpms_gpr_rna$mgi)) == 0){
    warning("Genes not found in data")
    return(NULL)
  }

  ## limit provided genes to those present in at least the RNA time course
  Genes <- Genes[which(Genes%in%tpms_gpr_rna$mgi)]

  ## create empty objects for the plot to then fill up
  gg_labels <- character()
  gg_frame <- data.frame(value=numeric(),
                         time=numeric(),
                         type=character(),
                         TPM=numeric(),
                         Gene=character(),
                         row.names=numeric())


  ## for each gene fill the ggplot table
  for(g in Genes){

    ## look up ensembl name for the corresponding gene
    Gene_ensembl <- row.names(tpms_gpr_rna)[which(tpms_gpr_rna$mgi==g)]

    ## add label in case info is provided
    if(Gene_info == TRUE){
      gg_labels <- c(gg_labels, paste0(g, "\n ", Gene_info_name,
                                     round(Gene_info_frame[which(Gene_info_frame[,1]==g),2][1],2)))
    }


    ## add original TPM measurements if wanted but convert them to log2FCs first
    if(plot_org_dots == TRUE){
      org_FCs <- log2(tpms_unsmoothed_rna[Gene_ensembl,c(1,2,5:21)]/tpms_gpr_rna[Gene_ensembl,1])
      }

  ## add data of TC after GPR on RNA level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=as.numeric(log2FCs_rna[Gene_ensembl,2:18]),
                         time=seq(0,32,2),
                         type="TC RNA log2FC GPR",
                         TPM=NA,
                         Gene=g,
                         row.names=1:length(seq(0,32,2))))

  ## add data of original TC change measurements in log2FCS on RNA level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=as.numeric(org_FCs),
                               time=c(0,seq(0,32,2),32),
                               type="TC RNA log2FC unsmoothed",
                               TPM=as.numeric(tpms_unsmoothed_rna[Gene_ensembl,
                                                                  c(1,2,5:21)]),
                               Gene=g,
                               row.names=(nrow(gg_frame)+1):
                                 (nrow(gg_frame)+length(c(0,seq(0,32,2),32)))))

  ## add data for MC 2i to 2i on RNA level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=c(0,as.numeric(log2FCs_MC_rna[Gene_ensembl,1:2])),
                               time=c(0,4,8),
                               type="MC RNA log2FC 2i",
                               TPM=NA,
                               Gene=g,
                               row.names=(nrow(gg_frame)+1):
                                 (nrow(gg_frame)+length(c(0,4,8)))))

  ## add data for MC 2i to N2B27 on RNA level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=c(0,as.numeric(log2FCs_MC_rna[Gene_ensembl,3:4])),
                               time=c(0,4,8),
                               type="MC RNA log2FC N2B27",
                               TPM=NA,
                               Gene=g,
                               row.names=(nrow(gg_frame)+1):
                                 (nrow(gg_frame)+length(c(0,4,8)))))

  ## add data for TC on protein level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=c(0,as.numeric(log2FCs_TC_prot[g,2:6])),
                               time=c(0,4,8,12,24,32),
                               type="TC Prot log2FC",
                               TPM=NA,
                               Gene=g,
                               row.names=(nrow(gg_frame)+1):
                                 (nrow(gg_frame)+length(c(0,4,8,12,24,32)))))

  ## add data for MC 2i to 2i on protein level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=c(0,as.numeric(log2FCs_MC_prot[g,10:11])),
                               time=c(0,4,8),
                               type="MC Prot log2FC 2i",
                               TPM=NA,
                               Gene=g,
                               row.names=(nrow(gg_frame)+1):
                                 (nrow(gg_frame)+length(c(0,4,8)))))

  ## add data for MC 2i to N2B27 on protein level
  gg_frame <- rbind(gg_frame,
                    data.frame(value=c(0,as.numeric(log2FCs_MC_prot[g,12:13])),
                               time=c(0,4,8),
                               type="MC Prot log2FC N2B27",
                               TPM=NA,
                               Gene=g,
                               row.names=(nrow(gg_frame)+1):
                                 (nrow(gg_frame)+length(c(0,4,8)))))
  }

  ## add Gene info in case of facetting
  if(Gene_info == TRUE){
    gg_frame$Gene <- factor(gg_frame$Gene, levels = Genes,
                     labels = gg_labels)
  }



 ## create output plot
  plot_out <- ggplot(gg_frame, aes(x = time, y = value, color = type)) +
    geom_point(data = gg_frame[gg_frame$type=="TC RNA log2FC unsmoothed",], color = "black") +
    geom_line(data = gg_frame[gg_frame$type=="TC RNA log2FC GPR",], size = 1.5, alpha = 0.8) +
    geom_line(data = gg_frame[gg_frame$type=="MC RNA log2FC 2i",], size = 1.5, lty = 2, alpha = 0.8) +
    geom_line(data = gg_frame[gg_frame$type=="MC RNA log2FC N2B27",], size = 1.5, alpha = 0.8) +
    geom_line(data = gg_frame[gg_frame$type=="TC Prot log2FC",], size = 1.5, alpha = 0.8) +
    geom_line(data = gg_frame[gg_frame$type=="MC Prot log2FC 2i",], size = 1.5, lty = 2, alpha = 0.8) +
    geom_line(data = gg_frame[gg_frame$type=="MC Prot log2FC N2B27",], size = 1.5, alpha = 0.8) +
    stat_peaks(data = gg_frame[gg_frame$type=="TC RNA log2FC unsmoothed",],
               span = NULL, geom = "text", color = "black",
               vjust = "bottom", position = position_nudge(y = 0.07),
               aes(label = round(TPM,0))) +
    stat_valleys(data = gg_frame[gg_frame$type=="TC RNA log2FC unsmoothed",],
                 span = NULL, geom = "text", color = "black",
                 vjust = "top", position = position_nudge(y = -0.07),
                 aes(label = round(TPM,0))) +
    labs(x = "Time",
         y = "log2FC",
         color ="Legend")+
    scale_color_manual(values = c(
      "TC RNA log2FC GPR" ="#FB8072",
      "MC RNA log2FC 2i" = "#FDB462",
      "MC RNA log2FC N2B27" = "#FDB462",
      "TC Prot log2FC" = "#80B1D3",
      "MC Prot log2FC 2i" =  "#8DD3C7",
      "MC Prot log2FC N2B27" = "#8DD3C7")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "right")

  if(length(Genes)==1){
      plot_out <- plot_out +
        ggtitle(paste0(title," \n",subtitle," \n",gg_labels))
    }
  if(length(Genes)>1){
    ##define labels for the grid and assign them
    plot_out <- plot_out +
      facet_wrap(~Gene, ncol = numcol, scales = "free") +
      scale_y_continuous(expand = expansion(mult = c(axis_ext, axis_ext)))+
      scale_x_continuous(expand = expansion(mult = c(axis_ext, axis_ext)))+
      ggtitle(paste0(title," \n",subtitle))
  }


 ## if targets are provided add grey lines indicating change of targets over TC
  ## this will lead to different scales as targets are z transformed while others aren't
    if(length(targets)>0){

      targets <- targets[which(targets%in%row.names(ztrans_log2FCs))]
      gg_frame <- rbind(gg_frame,
                        data.frame(value=c(t(ztrans_log2FCs[targets,])),
                                   time=rep(seq(0,32,2),length(targets)),
                                   type="TF targets",
                                   TPM=rep(targets, each = 17),
                                   Gene=g,
                                   row.names=(nrow(gg_frame)+1):
                                     (nrow(gg_frame)+length(rep(seq(0,32,2),length(targets))))))

      plot_out <- plot_out +
        geom_line(data = gg_frame[gg_frame$type=="TF targets",],
                  aes(x = time, y = value, group = TPM,
                      color = "black", alpha = 0.05))
      plot_out$layers <- plot_out$layers[c(length(plot_out$layers),
                                           1:(length(plot_out$layers)-1))]
    }

  return(plot_out)
}

## only if no TFs were removed by clusters
motif_cut_tree_interest_nonred <- motif_cut_tree[!duplicated(motif_cut_tree),]
motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster <- paste0(motif_cut_tree_interest_nonred$TF_cluster,
                                                                 ";",
                                                                 motif_cut_tree_interest_nonred$motif_cluster)

motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Mxi1;Mxd1"),
                               3] <- 7
motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Mxi1;Mxd1"),
                               4] <- "7;5"


motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Tfdp1;Tfdp2"),
                               3] <- 4
motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Tfdp1;Tfdp2"),
                               4] <- "4;2"

motif_cut_tree_interest_nonred$FPKM <- NA

for(r in 1:nrow(motif_cut_tree_interest_nonred)){
  motif_cut_tree_interest_nonred$FPKM[r] <- TF_clust_FPKM[motif_cut_tree_interest_nonred[r,2],2]
}
motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Mxi1;Mxd1"),
                               5] <- mean(TF_clust_FPKM[c("Mxi1", "Mxd1"),2])
motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Tfdp1;Tfdp2"),
                               5] <- mean(TF_clust_FPKM[c("Tfdp1", "Tfdp2"),2])

motif_cut_tree_interest_nonred$scaled_clust_FPKM <- NA
for(comb in unique(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster)){
  rows <- which(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster==
                  comb)
  motif_cut_tree_interest_nonred[rows,6] <- motif_cut_tree_interest_nonred[rows,5]/
    max(motif_cut_tree_interest_nonred[rows,5])
}

motif_cut_tree_interest_nonred <-
  motif_cut_tree_interest_nonred[order(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster),]

motif_cut_tree_interest_nonred$TPM_TC <- NA
for(r in 1:nrow(motif_cut_tree_interest_nonred)){
  motif_cut_tree_interest_nonred[r,7] <- mean(as.numeric(tpm_table_TC[which(tpm_table_TC$mgi==
                                                                   motif_cut_tree_interest_nonred[r,2]),1:2]))
}
motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Mxi1;Mxd1"),
                               7] <- mean(as.numeric(tpm_table_TC[which(tpm_table_TC$mgi%in%
                                                                          c("Mxi1", "Mxd1")),2]))
motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Tfdp1;Tfdp2"),
                               7] <- mean(as.numeric(tpm_table_TC[which(tpm_table_TC$mgi%in%
                                                                          c("Tfdp1", "Tfdp2")),2]))

motif_cut_tree_interest_nonred$scaled_TPM_TC <- NA
for(comb in unique(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster)){
  rows <- which(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster==
                  comb)
  motif_cut_tree_interest_nonred[rows,8] <- motif_cut_tree_interest_nonred[rows,7]/
    max(motif_cut_tree_interest_nonred[rows,7])
}

motif_cut_tree_interest_nonred$perc_changed_targets <- NA
for(r in 1:nrow(motif_cut_tree_interest_nonred)){
  targets_all <- unique(TFT_target_tables_list[[motif_cut_tree_interest_nonred[r,2]]][,1])
  targets_all <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs_all))]
  targets_changed <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs))]

  motif_cut_tree_interest_nonred$perc_changed_targets[r] <- length(targets_changed)/length(targets_all)
}
targets_all <- unique(c(TFT_target_tables_list[["Mxi1"]][,1],
                        TFT_target_tables_list[["Mxd1"]][,1]))
targets_all <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs_all))]
targets_changed <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs))]

motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Mxi1;Mxd1"),
                               9] <- length(targets_changed)/length(targets_all)

targets_all <- unique(c(TFT_target_tables_list[["Tfdp1"]][,1],
                        TFT_target_tables_list[["Tfdp2"]][,1]))
targets_all <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs_all))]
targets_changed <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs))]

motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF=="Tfdp1;Tfdp2"),
                               9] <- length(targets_changed)/length(targets_all)


### genes with high expression that represent motif target cluster combinations selected by hand ####
genes_selected <- c("Klf2", "E2f4", "Elf3", "Ascl2", "Sp1",
                    "Maz", "Zfp281", "Zfp341", "Zfp383",
                    "Rfx2", "Otx2", "E2f2", "Zscan10",
                    "Tfdp1", "Tfdp2", "Srebf2", "Stat6",
                    "Esrrb", "Rarg", "Sox2", "Sox4", "Sox3",
                    "Tcf3", "Tcf12", "Tcf4", "Neurod1", "Zic3",
                    "Nrf1", "Plagl1", "E2f3", "Dnmt1", "Zic2",
                    "Sall4", "Zfp740", "Max", "Mycn", "Myc")
expr_binary_TF_target_sel <- expr_binary_TF_target[which(row.names(expr_binary_TF_target)%in%
                                                           genes_selected),]
expr_binary_TF_target_sel <- expr_binary_TF_target_sel[,which(apply(expr_binary_TF_target_sel,
                                                                    2, max)>0)]

ha_sel = rowAnnotation(m1=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="1"),2]),
                   m2=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   c(motif_cut_tree[which(motif_cut_tree$motif_cluster=="2"),2],
                                     "Tfdp1", "Tfdp2")),
                   m3=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="3"),2]),
                   m4=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="4"),2]),
                   m5=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   c(motif_cut_tree[which(motif_cut_tree$motif_cluster=="5"),2],
                                     "Mxi1", "Mxd1")),
                   m6=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="6"),2]),
                   m7=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="7"),2]),
                   m8=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                   motif_cut_tree[which(motif_cut_tree$motif_cluster=="8"),2]),
                   col=list(m1=c("1"="#8DD3C7", "0"="white"),
                            m2=c("1"="#FFFFB3", "0"="white"),
                            m3=c("1"="#BEBADA", "0"="white"),
                            m4=c("1"="#FB8072", "0"="white"),
                            m5=c("1"="#80B1D3", "0"="white"),
                            m6=c("1"="#FDB462", "0"="white"),
                            m7=c("1"="#B3DE69", "0"="white"),
                            m8=c("1"="#FCCDE5", "0"="white")))
## Fig 5.36 A
Heatmap(expr_binary_TF_target_sel,
        show_row_names = T,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        row_split = TF_split_expr[which(row.names(TF_split_expr)%in%row.names(expr_binary_TF_target_sel)),],
        col=c("white","red"),
        name="target",
        border = TRUE,
        right_annotation = ha_sel)

motif_mat <- data.frame(m1=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                             motif_cut_tree[which(motif_cut_tree$motif_cluster=="1"),2]),
                        m2=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        c(motif_cut_tree[which(motif_cut_tree$motif_cluster=="2"),2],
                                          "Tfdp1", "Tfdp2")),
                        m3=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        motif_cut_tree[which(motif_cut_tree$motif_cluster=="3"),2]),
                        m4=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        motif_cut_tree[which(motif_cut_tree$motif_cluster=="4"),2]),
                        m5=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        c(motif_cut_tree[which(motif_cut_tree$motif_cluster=="5"),2],
                                          "Mxi1", "Mxd1")),
                        m6=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        motif_cut_tree[which(motif_cut_tree$motif_cluster=="6"),2]),
                        m7=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        motif_cut_tree[which(motif_cut_tree$motif_cluster=="7"),2]),
                        m8=as.numeric(row.names(expr_binary_TF_target_sel)%in%
                                        motif_cut_tree[which(motif_cut_tree$motif_cluster=="8"),2]))
row.names(motif_mat) <- row.names(expr_binary_TF_target_sel)

motif_mat_hclust <- hclust(dist(motif_mat), method = "ward.D2")
motif_mat_split <- as.data.frame(cutree(motif_mat_hclust, 14))

Heatmap(motif_mat,
        show_row_names = T,
        cluster_rows = F,
        cluster_columns = F,
        row_split = motif_mat_split,
        col=c("white","red"),
        name="cluster",
        border = TRUE)

expr_binary_TF_target_sel_hclust <- hclust(dist(expr_binary_TF_target_sel), method = "ward.D2")
## Fig 5.36 B
Heatmap(expr_binary_TF_target_sel,
        show_row_names = T,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        row_split = motif_mat_split[which(row.names(motif_mat_split)%in%row.names(expr_binary_TF_target_sel)),],
        col=c("white","red"),
        name="target",
        border = TRUE,
        right_annotation = ha_sel)

rowSums(expr_binary_TF_target_sel[,c("Pou3f1","Otx2",
                                      "Dnmt3a","Dnmt3b","Fgf5")])


expr_binary_TF_target_sel_short <- expr_binary_TF_target[c("E2f4","Zic3",#"Dnmt1",
                                                               "Sp1","Sall4",#"Zfp281",
                                                               "Esrrb","Klf4",#"Sox2",#"Ascl2","Myc",
                                                               "Tcf3","Sp5"),]
expr_binary_TF_target_sel_short <- expr_binary_TF_target_sel_short[,which(
  apply(expr_binary_TF_target_sel_short,2, max)>0)]
Heatmap(expr_binary_TF_target_sel_short,
        show_row_names = T,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        col=c("white","red"),
        name="target")

write.csv(expr_binary_TF_target_sel_short, "output_tables/binary_TF_target_final8_080124.csv",
          quote = F)

# for(t in 1:nrow(motif_cut_tree_interest_nonred)){
#   print(TF_overview_plot(tpms_gpr_rna = tpms_TC_gpr,
#                    tpms_unsmoothed_rna = tpm_table_TC,
#                    log2FCs_rna = log_norm_heatmap_gpr,
#                    log2FCs_MC_rna = medium_change_shrunken,
#                    log2FCs_TC_prot = prot_gpr,
#                    log2FCs_MC_prot = prot_MC,
#                    Genes=motif_cut_tree_interest_nonred[t,2],
#                    subtitle = paste0("taget; motif cluster:",
#                                      motif_cut_tree_interest_nonred[t,4],
#                                      "; scaled TPM: ",
#                                      round(motif_cut_tree_interest_nonred[t,9],2))))
# }



# TFs_vis <- c("Rbak", "Klf10", "Zfp740", "Zfp212", "Tcf3",
#              "Zfp661", "Tcfl5", "Dnmt1", "Zfp213",
#              "Zic3", "Zfp64", "Zfp566", "Ctcf",
#              "Rfx1", "Sall4", "E2f6", "E2f4", "E2f7",
#              "Zfp383", "Zfp341", "Klf5", "Zfp281", "Sp3", "Sp1")

# TFs_vis <- c("Trp53", "Zfp281", "Zfp423", "Myc", "Tcf7l1", "Rbpj", "Tet1")

TFs_vis <- c("Klf2", "E2f4", "Elf3", "Ascl2", "Sp1",
             "Maz", "Zfp281", "Zfp341", "Zfp383",
             "Rfx2", "Otx2", "E2f2", "Zscan10",
             "Tfdp1", "Tfdp2", "Srebf2", "Stat6",
             "Esrrb", "Rarg", "Sox2", "Sox4", "Sox3",
             "Tcf3", "Tcf12", "Tcf4", "Neurod1", "Zic3",
             "Nrf1", "Plagl1", "E2f3", "Dnmt1", "Zic2",
             "Sall4", "Zfp740", "Max", "Mycn", "Myc")
  #c("Zic3", "Tcf3", "E2f4", "Sp1", "Sall4", "Klf4", "Sp5", "Zfp281", "Zic2")

TFs_vis <- c("Zfp281", "Zfp383", "Maz")
TFs_vis <- c("Zfp740", "Sall4")
TFs_vis <- c("Zfp341", "Sp1")
TFs_vis <- c("Dnmt1", "Nrf1")
TFs_vis <- c("Tfdp1", "Tfdp2") ## not found?
TFs_vis <- c("Klf2", "E2f4")
TFs_vis <- c("Sox2", "Sox3", "Sox4")
TFs_vis <- c("Max", "Mycn", "Myc")
TFs_vis <- c("Tcf12", "Tcf4")
TFs_vis <- c("Zic3", "Zic2")
TFs_vis <- c("Tcf3", "Rarg")

for(t in TFs_vis){
  print(TF_overview_plot(tpms_gpr_rna = tpms_TC_gpr,
                         tpms_unsmoothed_rna = tpm_table_TC,
                         log2FCs_rna = log_norm_heatmap_gpr,
                         log2FCs_MC_rna = medium_change_shrunken,
                         log2FCs_TC_prot = prot_gpr,
                         log2FCs_MC_prot = prot_MC,
                         ztrans_log2FCs = ztrans_log2FCs,
                         #targets=unique(TFT_target_tables_list[[t]][,1]),
                         Genes=motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF==t)[1],2],
                         subtitle = paste0(motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF==t)[1],2],
                                           ": target; motif cluster:",
                                           motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF==t)[1],4],
                                           # "; scaled TPM: ",
                                           # round(motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF==t)[1],8],2),
                                           "; % targets changed: ",
                                           round(motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$TF==t)[1],9],2))
                         ))
}



unrelated_TFs <- names(TFT_target_tables_list)[which(!names(TFT_target_tables_list) %in%
                                                       unique(c(names(TFT_target_tables_list_rev_form_5),
                                                              names(TFT_target_tables_list_ext_form_5))))]

related_TFs <- unique(c(names(TFT_target_tables_list_rev_form_5),
                      names(TFT_target_tables_list_ext_form_5)))

gg_target_changed <- data.frame(perc_changed_targets=rep(NA, length(c(unrelated_TFs,related_TFs))),
                             type=rep(NA, length(c(unrelated_TFs,related_TFs))),
                             TF=rep(NA, length(c(unrelated_TFs,related_TFs))))

row.names(gg_target_changed) <- c(unrelated_TFs,related_TFs)
for(t in c(unrelated_TFs,related_TFs)){
  targets_all <- unique(TFT_target_tables_list[[t]][,1])
  targets_all <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs_all))]
  targets_changed <- targets_all[which(targets_all%in%row.names(ztrans_log2FCs))]

  gg_target_changed[t,1] <- length(targets_changed)/length(targets_all)
  if(t %in% related_TFs){
    gg_target_changed[t,2] <- "has formative extended targets"
  }else{
    gg_target_changed[t,2] <- "has no formative extended targets"
  }
  gg_target_changed[t,3] <- t
}

ggplot(gg_target_changed, aes(x = perc_changed_targets,
                              group = type, color = type, fill = type)) +
  geom_density(alpha = 0.3) +
  theme_bw()

for(comb in unique(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster)){
  print(comb)
  Genes <- motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster==comb),2]
  if(length(Genes) < 6){
    col_num <- 2
  }else if(length(Genes) < 12){
    col_num <- 3
  }else{
    col_num <- 4
  }

  print(TF_overview_plot(tpms_gpr_rna = tpms_TC_gpr,
                   tpms_unsmoothed_rna = tpm_table_TC,
                   log2FCs_rna = log_norm_heatmap_gpr,
                   log2FCs_MC_rna = medium_change_shrunken,
                   log2FCs_TC_prot = prot_gpr,
                   log2FCs_MC_prot = prot_MC,
                   ztrans_log2FCs = ztrans_log2FCs,
                   Gene_info=TRUE,
                   Gene_info_frame = motif_cut_tree_interest_nonred[which(motif_cut_tree_interest_nonred$combi_TF_cluster_motif_cluster==comb),c(2,8)],
                   Genes=Genes,
                   title="target cluster; motif cluster",subtitle = comb,
                   numcol = col_num))
}

CLIM2 <- readRDS("data/RDS_data/CLIM2_Kdm6a_updated_091020.rds")
KOvKO <- CLIM2$singlecontrasts$fit_KOVKO$coefficients

## how do all targets (not only extended formative) behave in TC

#write.xlsx(fimo_list[["naive_formative_all_3"]],
#           file = "R_outputs/fimo_naive_formative_all_3_peaks.xlsx")

### distance of extended to next gene ####
target_split <- cutree(clust_binary_form_targets, 2)

low_con_targets <- rep(0,length(biomart102[which(biomart102$mgi_symbol%in%
                                                names(target_split[which(target_split==1)])),1]))
names(low_con_targets) <- biomart102[which(biomart102$mgi_symbol%in%
                                          names(target_split[which(target_split==1)])),1]


high_con_targets <- rep(0,length(biomart102[which(biomart102$mgi_symbol%in%
                                                names(target_split[which(target_split==2)])),1]))
names(high_con_targets) <- biomart102[which(biomart102$mgi_symbol%in%
                                          names(target_split[which(target_split==2)])),1]

for(gene in names(low_con_targets)){
  chromosome <- gene_position_102[gene,2]
  tss <- gene_position_102[gene,4]

  genes_chrom <- gene_position_102[which(gene_position_102$chromosome==chromosome),]

  low_con_targets[gene] <- mean(sort(abs(tss-genes_chrom[,4]))[2])
}

for(gene in names(high_con_targets)){
  chromosome <- gene_position_102[gene,2]
  tss <- gene_position_102[gene,4]

  genes_chrom <- gene_position_102[which(gene_position_102$chromosome==chromosome),]

  high_con_targets[gene] <- mean(sort(abs(tss-genes_chrom[,4]))[2])
}

summary(low_con_targets)
summary(high_con_targets)
ks.test(low_con_targets,high_con_targets)

