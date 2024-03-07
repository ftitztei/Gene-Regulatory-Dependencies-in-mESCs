
library(Seurat)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(biomaRt)


matrix_dir = "data/998929_RC9_filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
#remove the -1 on cell IDs
colnames(mat) <-  sapply(strsplit(colnames(mat), split = "-", fixed=T), function(x) (x[1]))
head(colnames(mat))

tags <- read.csv("data/998929_RC9_filtered_feature_bc_matrix/RC9_seurat_hashID.csv")
table(tags$hash.ID)

mitochondrial <- c("ENSMUSG00000064348", "ENSMUSG00000064344",
                   "ENSMUSG00000064365", "ENSMUSG00000064366",
                   "ENSMUSG00000064371", "ENSMUSG00000064339",
                   "ENSMUSG00000064350", "ENSMUSG00000064340",
                   "ENSMUSG00000064349", "ENSMUSG00000064338",
                   "ENSMUSG00000064363", "ENSMUSG00000064367",
                   "ENSMUSG00000064361", "ENSMUSG00000064354",
                   "ENSMUSG00000064358", "ENSMUSG00000064352",
                   "ENSMUSG00000064337", "ENSMUSG00000064359",
                   "ENSMUSG00000064346", "ENSMUSG00000064353",
                   "ENSMUSG00000064357", "ENSMUSG00000064370",
                   "ENSMUSG00000064345", "ENSMUSG00000064364",
                   "ENSMUSG00000064351", "ENSMUSG00000064368",
                   "ENSMUSG00000064369", "ENSMUSG00000064336",
                   "ENSMUSG00000064341", "ENSMUSG00000064355",
                   "ENSMUSG00000064347", "ENSMUSG00000064343",
                   "ENSMUSG00000064356", "ENSMUSG00000064342",
                   "ENSMUSG00000064360", "ENSMUSG00000065947",
                   "ENSMUSG00000064372")

SC_RC9 <- CreateSeuratObject(mat[,which(colnames(mat)%in%tags$X)],
                             project = "SC_RC9", assay = "RNA")


SC_RC9[["percent.mt"]] <- PercentageFeatureSet(SC_RC9, features = mitochondrial[which(mitochondrial%in%
                                                                                       row.names(SC_RC9))])
SC_RC9[["nonmt.libsize"]] <- colSums(SC_RC9@assays$RNA@counts[which(!row.names(SC_RC9@assays$RNA@counts)%in%
                                                              mitochondrial),])

sum(tags[which(tags$X%in%colnames(SC_RC9)),1]==colnames(SC_RC9))
SC_RC9[["type"]] <- tags[which(tags$X%in%colnames(SC_RC9)),2]

VlnPlot(SC_RC9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,)

gg_qc <- cbind(SC_RC9[["nFeature_RNA"]],SC_RC9[["nCount_RNA"]],
               SC_RC9[["percent.mt"]],SC_RC9[["type"]])
gg_qc$type_fact <- as.factor(gg_qc$type)
levels(gg_qc$type_fact) <- levels(gg_qc$type_fact)[c(1,5,2,3,4,6,7)]

ggplot(gg_qc,aes(x=nCount_RNA,y=percent.mt,
                 alpha = 0.2, color=type_fact, fill=type_fact)) +
  geom_point() +
  theme_bw()

ggplot(gg_qc,aes(x=nFeature_RNA,y=percent.mt,
                 alpha = 0.2, color=type_fact, fill=type_fact)) +
  geom_point()+
  theme_bw()

ggplot(gg_qc,aes(x=nCount_RNA,y=nFeature_RNA,
                 alpha = 0.2, color=type_fact, fill=type_fact)) +
  geom_point()+
  theme_bw()
## Fig 5.21 B
ggplot(gg_qc,aes(x=type_fact,y=nFeature_RNA,
                 alpha = 0.2, color=type_fact, fill=type_fact)) +
  geom_violin()+
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "grey","black"))+
  scale_fill_manual(values=c("darkgoldenrod1","darkorange1",
                             "darkorange3","coral3",
                             "darkred","grey","black"))+
  theme_bw()

## Fig 5.21 A
ggplot(gg_qc,aes(x=type_fact,y=nCount_RNA,
                 alpha = 0.2, color=type_fact, fill=type_fact)) +
  geom_violin()+
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                               "darkorange3","coral3",
                              "darkred", "grey","black"))+
  scale_fill_manual(values=c("darkgoldenrod1","darkorange1",
                               "darkorange3","coral3",
                              "darkred","grey","black"))+
  theme_bw()

ggplot(gg_qc,aes(x=nCount_RNA,
                 alpha = 0.2, color=type_fact, fill=type_fact)) +
  geom_density()+
  theme_bw()

summary(gg_qc$nCount_RNA[which(gg_qc$type=="Doublet")])
summary(gg_qc$nCount_RNA[which(gg_qc$type!="Doublet")])

length(gg_qc$nCount_RNA[which(gg_qc$type!="Doublet")])

dim(subset(SC_RC9,
           subset = type != "Doublet" &
             nFeature_RNA > 200 &
             nFeature_RNA < 3500 &
             percent.mt < 20))

dim(subset(SC_RC9,
           subset = type != "Doublet" &
             nFeature_RNA > 200 &
             nFeature_RNA < 3500 &
             nCount_RNA > 1000 &
             nCount_RNA < 8000 &
             percent.mt < 20))


SC_RC9 <- subset(SC_RC9,
                 subset = type != "Doublet" &
                   nFeature_RNA > 200 &
                   nFeature_RNA < 3500 &
                   #nCount_RNA > 1000 &
                   #nCount_RNA < 9000 &
                   percent.mt < 20)

CP10K <- NormalizeData(SC_RC9, normalization.method = "RC",
                       scale.factor = 10000)@assays$RNA@data
CP10K <- as.matrix(CP10K)
CP10K <- CP10K+0.1

SC_RC9 <- NormalizeData(SC_RC9, normalization.method = "LogNormalize",
                        scale.factor = 10000)


#SC_RC9 <- FindVariableFeatures(SC_RC9, selection.method = "vst", nfeatures = 2000)

all.genes <- row.names(SC_RC9)
SC_RC9 <- ScaleData(SC_RC9, features = all.genes)

SC_RC9 <- SCTransform(SC_RC9, method = "glmGamPoi",
                      vars.to.regress = c("nFeature_RNA","percent.mt","nonmt.libsize"),
                      return.only.var.genes = F)#"nFeature_RNA","percent.mt",

CP10K_sct <- SC_RC9@assays[["SCT"]]@counts
CP10K_sct <- apply(CP10K_sct, 2, function(x) x/(sum(x)/10000))
CP10K_sct <- CP10K_sct + 0.1

count_sct <- SC_RC9@assays[["SCT"]]@counts
count_sct <- count_sct + 0.1

SC_RC9 <- RunPCA(SC_RC9,
                 features = VariableFeatures(object = SC_RC9))
eigValues <- (SC_RC9@reductions$pca@stdev)**2
varExpl <- eigValues/sum(eigValues)

gg_pca <- data.frame(sample=row.names(SC_RC9@reductions$pca@cell.embeddings),
                     PC1=SC_RC9@reductions$pca@cell.embeddings[,1],
                     PC2=SC_RC9@reductions$pca@cell.embeddings[,2],
                     PC3=SC_RC9@reductions$pca@cell.embeddings[,3],
                     PC4=SC_RC9@reductions$pca@cell.embeddings[,4],
                     PC5=SC_RC9@reductions$pca@cell.embeddings[,5],
                     time=NA,
                     perc.mt=NA,
                     nonmt.libsize=SC_RC9[["nonmt.libsize"]],
                     number.genes=SC_RC9@meta.data$nFeature_RNA)
for(r in 1:dim(gg_pca)[1]){
  gg_pca[r,7] <- as.character(tags[which(tags$X==gg_pca[r,1]),2])
  gg_pca[r,8] <- SC_RC9[["percent.mt"]][which(row.names(SC_RC9[["percent.mt"]])==gg_pca[r,1]),1]
}
gg_pca$time <- factor(gg_pca$time, levels = c("0h","6h","12h",
                                              "24h","48h","Negative"))

## Fig 5.23 A
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=time, fill=time),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))

ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=perc.mt, fill=perc.mt),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))
## Fig 5.22 A
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=log10(nonmt.libsize), fill=log10(nonmt.libsize)),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))
## Fig 5.22B
ggplot(gg_pca)+
  geom_point(aes(x=PC1,y=PC2,size=4, color=log10(number.genes), fill=log10(number.genes)),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC1: ",round(varExpl[1]*100,1),"%"),
       y=paste0("PC2: ",round(varExpl[2]*100,1),"%"))


## Fig 5.23B
ggplot(gg_pca)+
  geom_point(aes(x=PC2,y=PC3,size=4, color=time, fill=time),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC2: ",round(varExpl[2]*100,1),"%"),
       y=paste0("PC3: ",round(varExpl[3]*100,1),"%"))

ggplot(gg_pca)+
  geom_point(aes(x=PC3,y=PC4,size=4, color=time, fill=time),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x=paste0("PC3: ",round(varExpl[3]*100,1),"%"),
       y=paste0("PC4: ",round(varExpl[4]*100,1),"%"))

ggplot(gg_pca)+
  geom_point(aes(x=PC3,y=PC4,size=4, color=log10(nonmt.libsize), fill=log10(nonmt.libsize)),
             size=2, alpha=0.3) +
  theme_bw() +
  labs(x=paste0("PC3: ",round(varExpl[3]*100,1),"%"),
       y=paste0("PC4: ",round(varExpl[4]*100,1),"%"))


SC_RC9 <- FindNeighbors(SC_RC9, dims = 1:20, verbose = F)
SC_RC9 <- FindClusters(SC_RC9, resolution = 0.5, verbose = F)

SC_RC9 <- RunUMAP(SC_RC9, dims = 1:20, verbose = F)
gg_pca$UMAP_1 <- SC_RC9@reductions$umap@cell.embeddings[,1]
gg_pca$UMAP_2 <- SC_RC9@reductions$umap@cell.embeddings[,2]

ggplot(gg_pca)+
  geom_point(aes(x=UMAP_1,y=UMAP_2,size=4, color=time, fill=time),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x="UMAP", y="UMAP")

SC_RC9 <- RunTSNE(SC_RC9, dims = 1:20, verbose = FALSE)
gg_pca$tSNE_1 <- SC_RC9@reductions$tsne@cell.embeddings[,1]
gg_pca$tSNE_2 <- SC_RC9@reductions$tsne@cell.embeddings[,2]
## Fig 5.23 C
ggplot(gg_pca)+
  geom_point(aes(x=tSNE_1,y=tSNE_2,size=4, color=time, fill=time),
             size=2, alpha=0.3) +
  scale_color_manual(values=c("darkgoldenrod1","darkorange1",
                              "darkorange3","coral3",
                              "darkred", "black")) +
  theme_bw() +
  labs(x="tSNE", y="tSNE")


## retrieve a lookup between ensembl gene ids and mgi symbols from biomart
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

attributes <- listAttributes(ensembl)
biomart <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "mgi_symbol"),
                 mart=ensembl)


mgi_CP10K <- rep(NA,nrow(CP10K))
biotype_CP10K <- rep(NA,nrow(CP10K))
for(d in 1:dim(CP10K)[1]){
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),2])==1){
    biotype_CP10K[d] <- biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),2]
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),2])>1){
    mart_types <- biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),2]
    if(length(unique(mart_types))==1){
      biotype_CP10K[d] <- mart_types[1]
    }
    else{
      print(d)
    }
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),2])==0){
    biotype_CP10K[d] <- NA
  }
  if(length(biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),3])==TRUE){
    mgi_CP10K[d] <- biomart[which(biomart$ensembl_gene_id==rownames(CP10K)[d]),3]
  }
  else{
    mgi_CP10K[d] <- NA
  }
}
##remove all genes with NA as mgi symbol
CP10K <- CP10K[!is.na(mgi_CP10K),]
biotype_CP10K <- biotype_CP10K[!is.na(mgi_CP10K)]
mgi_CP10K <- mgi_CP10K[!is.na(mgi_CP10K)]


##remove all genes with empty mgi symbol
CP10K <- CP10K[which(mgi_CP10K!=""),]
biotype_CP10K <- biotype_CP10K[which(mgi_CP10K!="")]
mgi_CP10K <- mgi_CP10K[which(mgi_CP10K!="")]


##remove all genes with duplicate mgi symbols
CP10K <- CP10K[mgi_CP10K%in%names(table(mgi_CP10K)[which(table(mgi_CP10K)==1)]),]
biotype_CP10K <- biotype_CP10K[mgi_CP10K%in%names(table(mgi_CP10K)[which(table(mgi_CP10K)==1)])]
mgi_CP10K <- mgi_CP10K[mgi_CP10K%in%names(table(mgi_CP10K)[which(table(mgi_CP10K)==1)])]

CP10K <- CP10K[which(biotype_CP10K%in%
                                   c("protein_coding","lincRNA",
                                     "processed_transcript","antisense",
                                     "3prime_overlapping_ncRNA",
                                     "bidirectional_promoter_lncRNA",
                                     "macro_lncRNA","miRNA",
                                     "misc_RNA","lncRNA",
                                     "scaRNA","scRNA",
                                     "sense_intronic","sense_overlapping",
                                     "snoRNA","snRNA","sRNA")),]
CP10K_sct <- CP10K_sct[which(row.names(CP10K_sct)%in%
                               row.names(CP10K)),]
count_sct <- count_sct[which(row.names(count_sct)%in%
                               row.names(CP10K)),]
count_sct <- as.matrix(count_sct)

mgi_CP10K <- mgi_CP10K[which(biotype_CP10K%in%
                                   c("protein_coding","lincRNA",
                                     "processed_transcript","antisense",
                                     "3prime_overlapping_ncRNA",
                                     "bidirectional_promoter_lncRNA",
                                     "macro_lncRNA","miRNA",
                                     "misc_RNA","lncRNA",
                                     "scaRNA","scRNA",
                                     "sense_intronic","sense_overlapping",
                                     "snoRNA","snRNA","sRNA"))]
biotype_CP10K <- biotype_CP10K[which(biotype_CP10K%in%
                                         c("protein_coding","lincRNA",
                                           "processed_transcript","antisense",
                                           "3prime_overlapping_ncRNA",
                                           "bidirectional_promoter_lncRNA",
                                           "macro_lncRNA","miRNA",
                                           "misc_RNA","lncRNA",
                                           "scaRNA","scRNA",
                                           "sense_intronic","sense_overlapping",
                                           "snoRNA","snRNA","sRNA"))]

## create gain loss frame
gain_loss_SC <- data.frame(TP_0 = rep(NA,nrow(count_sct)),
                           TP_6 = rep(NA,nrow(count_sct)),
                           TP_12 = rep(NA,nrow(count_sct)),
                           TP_24 = rep(NA,nrow(count_sct)),
                           TP_48 = rep(NA,nrow(count_sct)),
                           change_sign_all = rep(NA,nrow(count_sct)),
                           log2FC_avg = rep(NA,nrow(count_sct)),
                           log2FC_6_0 = rep(NA,nrow(count_sct)),
                           log2FC_12_0 = rep(NA,nrow(count_sct)),
                           log2FC_24_0 = rep(NA,nrow(count_sct)),
                           log2FC_48_0 = rep(NA,nrow(count_sct)),
                           change_sign = rep(NA,nrow(count_sct)),
                           row.names = row.names(count_sct))

for(d in 1:nrow(gain_loss_SC)){
  ## T0
  gain_loss_SC[d,1] <- mean(as.numeric(count_sct[d,which(SC_RC9@meta.data$type=="0h")]))
  ## T6
  gain_loss_SC[d,2] <- mean(as.numeric(count_sct[d,which(SC_RC9@meta.data$type=="6h")]))
  ## T12
  gain_loss_SC[d,3] <- mean(as.numeric(count_sct[d,which(SC_RC9@meta.data$type=="12h")]))
  ## T24
  gain_loss_SC[d,4] <- mean(as.numeric(count_sct[d,which(SC_RC9@meta.data$type=="24h")]))
  ## T48
  gain_loss_SC[d,5] <- mean(as.numeric(count_sct[d,which(SC_RC9@meta.data$type=="48h")]))

  ## is there a consistent gain in expression
  if(sign(gain_loss_SC[d,2]-gain_loss_SC[d,1])%in%c(0,1) &
     sign(gain_loss_SC[d,3]-gain_loss_SC[d,1])%in%c(0,1) &
     sign(gain_loss_SC[d,4]-gain_loss_SC[d,1])%in%c(0,1) #&
     #sign(gain_loss_SC[d,5]-gain_loss_SC[d,1])%in%c(0,1)
     ){
    gain_loss_SC[d,6] <- 1
  }
  ## is there a consistent loss in expression
  if(sign(gain_loss_SC[d,2]-gain_loss_SC[d,1])%in%c(0,-1) &
     sign(gain_loss_SC[d,3]-gain_loss_SC[d,1])%in%c(0,-1) &
     sign(gain_loss_SC[d,4]-gain_loss_SC[d,1])%in%c(0,-1) #&
     #sign(gain_loss_SC[d,5]-gain_loss_SC[d,1])%in%c(0,-1)
     ){
    gain_loss_SC[d,6] <- -1
  }
  if(sign(gain_loss_SC[d,2]-gain_loss_SC[d,1])==0 &
     sign(gain_loss_SC[d,3]-gain_loss_SC[d,1])==0 &
     sign(gain_loss_SC[d,4]-gain_loss_SC[d,1])==0 #&
     #sign(gain_loss_SC[d,5]-gain_loss_SC[d,1])==0
     ){
    gain_loss_SC[d,6] <- 0
  }
  ## log2 FC all others vs T0
  gain_loss_SC[d,7] <- mean(c(log2((gain_loss_SC[d,2])/(gain_loss_SC[d,1])),
                              log2((gain_loss_SC[d,3])/(gain_loss_SC[d,1])),
                              log2((gain_loss_SC[d,4])/(gain_loss_SC[d,1]))))
  ## log2 FC all T6 vs T0
  gain_loss_SC[d,8] <- log2((gain_loss_SC[d,2])/(gain_loss_SC[d,1]))
  ## log2 FC all T12 vs T0
  gain_loss_SC[d,9] <- log2((gain_loss_SC[d,3])/(gain_loss_SC[d,1]))
  ## log2 FC all T24 vs T0
  gain_loss_SC[d,10] <- log2((gain_loss_SC[d,4])/(gain_loss_SC[d,1]))
  ## log2 FC all T48 vs T0
  gain_loss_SC[d,11] <- log2((gain_loss_SC[d,5])/(gain_loss_SC[d,1]))

}

gain_loss_ggplot <- gain_loss_SC[,6:11]
gain_loss_ggplot <- gain_loss_ggplot[rowSums(is.na(gain_loss_ggplot))!=ncol(gain_loss_ggplot),]
## plot log2 FC T48 vs T0 against average log2 FC vs T0
ggplot(gain_loss_ggplot, aes(x=log2FC_avg, y=log2FC_48_0,
                             color=as.factor(change_sign_all),
                             fill=as.factor(change_sign_all))) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              col="red") +
  labs(x="average log2 FC vs 0h",
       y="log2 FC 48h vs 0h")
## plot log2 FC T48 vs T0 against log2 FC T24 vs T0
ggplot(gain_loss_ggplot, aes(x=log2FC_24_0, y=log2FC_48_0,
                             color=as.factor(change_sign_all),
                             fill=as.factor(change_sign_all))) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0,
              col="red") +
  labs(x="log2 FC 24h vs 0h",
       y="log2 FC 48h vs 0h")

gain_loss_ggplot_sct <- gain_loss_ggplot[which(row.names(gain_loss_ggplot)%in%
                                                 row.names(count_sct)),]
mgi_CP10K_sct <- mgi_CP10K[which(row.names(CP10K)%in%
                                   row.names(CP10K_sct))]
mgi_count_sct <- mgi_CP10K[which(row.names(CP10K)%in%
                                   row.names(count_sct))]

saveRDS(data.frame(ensembl=row.names(CP10K),
                   mgi=mgi_CP10K),"RDS/log2_CP10K_Christa_lookup.rds")
saveRDS(data.frame(ensembl=row.names(CP10K_sct),
                   mgi=mgi_CP10K_sct),"RDS/log2_CP10K_sct_Christa_lookup.rds")
saveRDS(data.frame(ensembl=row.names(count_sct),
                   mgi=mgi_count_sct),"RDS/count_sct_Christa_lookup.rds")

saveRDS(log2(CP10K),"RDS/log2_CP10K_Christa.rds")
saveRDS(log2(CP10K_sct),"RDS/log2_CP10K_sct_Christa.rds")
saveRDS(count_sct,"RDS/count_sct_Christa.rds")
saveRDS(SC_RC9,"RDS/SC_RC9_Christa.rds")

saveRDS(gain_loss_ggplot,"RDS/gain_loss_ggplot_Christa.rds")
saveRDS(gain_loss_ggplot_sct,"RDS/gain_loss_ggplot_sct_Christa.rds")

