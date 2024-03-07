## In this script Gaussian Process Regression (GPR) is applied to the time course data.
## As there are no replicates from 2 to 30 hours and we assume the same variance
## accross different time points of the same gene we used GPR to smoothen the
## log2 FCs for all genes. Here information from neighbouring time points
## is used to infer a value for every time point in the time series.

library(tgp)

### load objects from previos scripts ####
log_norm_counts_short <- readRDS("RDS/log_norm_counts_short.rds")
shrunken_foldchanges_complete <- readRDS("RDS/shrunken_foldchanges_complete.rds")
tpm_table <- readRDS("RDS/tpm_table.rds")

## As we tested different workflows of applying GPR to the data we
## will run the smoothing three different times.
## 1. GPR was run on log2 normalized counts
## 2. GPR was used on shrunken log2FCs from DESeq2 (used in later analysis)
## 3. GPR was used on shrunken log2FCs from DESeq2 after
## removing the 6 hours time point as we saw some outliers here

### gaussian process regression log2 normal counts ####
## prepare table to fill with GPR outcomes
log_norm_heatmap_gpr_counts <- matrix(NA, nrow=dim(log_norm_counts_short)[1], ncol=21)
row.names(log_norm_heatmap_gpr_counts) <- row.names(log_norm_counts_short)
colnames(log_norm_heatmap_gpr_counts) <- c("RC9_2iL_1", "RC9_2iL_2",
                                           "RC9_2i_1", "RC9_2i_2", "RC9_2h",
                                           "RC9_4h", "RC9_6h", "RC9_8h",
                                           "RC9_10h", "RC9_12h", "RC9_14h",
                                           "RC9_16h", "RC9_18h", "RC9_20h",
                                           "RC9_22h", "RC9_24h", "RC9_26h",
                                           "RC9_28h", "RC9_30h",
                                           "RC9_32h_1","RC9_32h_2")

## fill 2iLIF values as they will not be included in the GPR
for(d in 1:dim(log_norm_counts_short)[1]){
  log_norm_heatmap_gpr_counts[d,c(1:2)] <- log_norm_counts_short[d,c(1:2)]
}

## prepare table for upper prediction limit
q1_gpr <- matrix(NA, nrow=dim(log_norm_counts_short)[1], ncol=19)
row.names(q1_gpr) <- row.names(log_norm_counts_short)
colnames(q1_gpr) <- c("RC9_2i_1", "RC9_2i_2", "RC9_2h",
                      "RC9_4h", "RC9_6h", "RC9_8h",
                      "RC9_10h", "RC9_12h", "RC9_14h",
                      "RC9_16h", "RC9_18h", "RC9_20h",
                      "RC9_22h", "RC9_24h", "RC9_26h",
                      "RC9_28h", "RC9_30h",
                      "RC9_32h_1","RC9_32h_2")
## prepare table for lower prediction limit
q2_gpr <- matrix(NA, nrow=dim(log_norm_counts_short)[1], ncol=19)
row.names(q2_gpr) <- row.names(log_norm_counts_short)
colnames(q2_gpr) <- c("RC9_2i_1", "RC9_2i_2", "RC9_2h",
                      "RC9_4h", "RC9_6h", "RC9_8h",
                      "RC9_10h", "RC9_12h", "RC9_14h",
                      "RC9_16h", "RC9_18h", "RC9_20h",
                      "RC9_22h", "RC9_24h", "RC9_26h",
                      "RC9_28h", "RC9_30h",
                      "RC9_32h_1","RC9_32h_2")


X <- data.frame(time=c(0,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,32))

## prepare a list where gpr models are saved
gpr_models <- vector(mode="list",length=dim(log_norm_heatmap_gpr_counts)[1])
names(gpr_models) <- row.names(log_norm_heatmap_gpr_counts)

## run through each gene, build a GPR model and save mean, upper and lower
## boundaries
for(d in 1:dim(log_norm_heatmap_gpr_counts)[1]){
  if(d %% 1000 == 0){
    message(d)
  }

  if(all(log_norm_counts_short[d,c(3:21)]==0)){
    next
  }

  Y <- data.frame(reads=log_norm_counts_short[d,c(3:21)])
  set.seed(1234)
  gpr <- bgp(X,Y,verb=0,bprior = "bmznot")
  gpr_models[[d]] <- gpr
  mean <- mean(log_norm_counts_short[d,c(3:21)])

  log_norm_heatmap_gpr_counts[d,c(3:21)] <- gpr[["Zp.mean"]]
  q1_gpr[d,] <- gpr[["Zp.q1"]]
  q2_gpr[d,] <- gpr[["Zp.q2"]]

}

## remove incomplete cases (genes)
log_norm_heatmap_gpr_counts <- log_norm_heatmap_gpr_counts[complete.cases(log_norm_heatmap_gpr_counts),]
q1_gpr <- q1_gpr[complete.cases(q1_gpr),]
q2_gpr <- q2_gpr[complete.cases(q2_gpr),]

## combine all three measuremnts in one list
gpr_list_counts <- list(log_norm_heatmap_gpr_counts,q1_gpr,q2_gpr)

### gaussian process regression normal dist deseq shrunken FCs ####
## prepare table to fill with GPR outcomes
log_norm_heatmap_gpr_shrunkenFCS <- matrix(NA, nrow=dim(shrunken_foldchanges_complete)[1], ncol=18)
row.names(log_norm_heatmap_gpr_shrunkenFCS) <- row.names(shrunken_foldchanges_complete)
colnames(log_norm_heatmap_gpr_shrunkenFCS) <- c("RC9_2iL", "RC9_2i", "RC9_2h",
                                                "RC9_4h", "RC9_6h", "RC9_8h",
                                                "RC9_10h", "RC9_12h", "RC9_14h",
                                                "RC9_16h", "RC9_18h", "RC9_20h",
                                                "RC9_22h", "RC9_24h", "RC9_26h",
                                                "RC9_28h", "RC9_30h", "RC9_32h")
## fill 2iLIF values as they will not be included in the GPR
for(d in 1:dim(shrunken_foldchanges_complete)[1]){
  log_norm_heatmap_gpr_shrunkenFCS[d,1] <- shrunken_foldchanges_complete[d,1]
}

## prepare table for upper prediction limit
q1_gpr <- matrix(NA, nrow=dim(shrunken_foldchanges_complete)[1], ncol=17)
row.names(q1_gpr) <- row.names(shrunken_foldchanges_complete)
colnames(q1_gpr) <- c("RC9_2i", "RC9_2h",
                      "RC9_4h", "RC9_6h", "RC9_8h",
                      "RC9_10h", "RC9_12h", "RC9_14h",
                      "RC9_16h", "RC9_18h", "RC9_20h",
                      "RC9_22h", "RC9_24h", "RC9_26h",
                      "RC9_28h", "RC9_30h", "RC9_32h")
## prepare table for lower prediction limit
q2_gpr <- matrix(NA, nrow=dim(shrunken_foldchanges_complete)[1], ncol=17)
row.names(q2_gpr) <- row.names(shrunken_foldchanges_complete)
colnames(q2_gpr) <- c("RC9_2i", "RC9_2h",
                      "RC9_4h", "RC9_6h", "RC9_8h",
                      "RC9_10h", "RC9_12h", "RC9_14h",
                      "RC9_16h", "RC9_18h", "RC9_20h",
                      "RC9_22h", "RC9_24h", "RC9_26h",
                      "RC9_28h", "RC9_30h", "RC9_32h")


X <- data.frame(time=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32))

## prepare a list where gpr models are saved
gpr_models <- vector(mode="list",length=dim(log_norm_heatmap_gpr_shrunkenFCS)[1])
names(gpr_models) <- row.names(log_norm_heatmap_gpr_shrunkenFCS)

## run through each gene, build a GPR model and save mean, upper and lower
## boundaries
for(d in 1:dim(log_norm_heatmap_gpr_shrunkenFCS)[1]){
  if(d %% 1000 == 0){
    message(d)
  }

  if(all(shrunken_foldchanges_complete[d,c(2:18)]==0)){
    next
  }

  Y <- data.frame(reads=shrunken_foldchanges_complete[d,c(2:18)])
  set.seed(1234)
  gpr <- bgp(X,Y,verb=0,bprior = "bmznot")
  gpr_models[[d]] <- gpr
  mean <- mean(shrunken_foldchanges_complete[d,c(2:18)])

  log_norm_heatmap_gpr_shrunkenFCS[d,c(2:18)] <- gpr[["Zp.mean"]]
  q1_gpr[d,] <- gpr[["Zp.q1"]]
  q2_gpr[d,] <- gpr[["Zp.q2"]]

}

## remove incomplete cases (genes)
log_norm_heatmap_gpr_shrunkenFCS <- log_norm_heatmap_gpr_shrunkenFCS[complete.cases(log_norm_heatmap_gpr_shrunkenFCS),]
q1_gpr <- q1_gpr[complete.cases(q1_gpr),]
q2_gpr <- q2_gpr[complete.cases(q2_gpr),]

## combine all three measuremnts in one list
gpr_list_shrunkenFCS <- list(log_norm_heatmap_gpr_shrunkenFCS,q1_gpr,q2_gpr)


### TPMs based on GPR ####
## calculate TPMs based on GPR results by using the TPMs at 2i
## TPMs for other time points are calculated by log2FCs of this time point to 2i
## and the corresponding TPMs at 2i
tpms_TC <- tpm_table
tpms_TC_gpr <- tpms_TC[which(row.names(tpms_TC)%in%row.names(gpr_list_shrunkenFCS[[1]])),]

for(d in 1:dim(tpms_TC_gpr)[1]){
  tpms_TC_gpr[d,1:21] <- (2**gpr_list_shrunkenFCS[[1]][which(row.names(gpr_list_shrunkenFCS[[1]])==row.names(tpms_TC_gpr)[d])
                                           ,c(2,2,1,1,3:17,18,18)])*mean(c(tpms_TC_gpr[d,1],tpms_TC_gpr[d,2]))
}

tpms_TC_gpr$max <- apply(tpms_TC_gpr[,1:21],1,max)
head(tpms_TC_gpr)

### >>> remove 6 hrs and test induction and repression TPs <<< ####


### export rds files in folder for other scripts ####

saveRDS(gpr_list_counts,"RDS/gpr_list_counts.rds")

saveRDS(gpr_list_shrunkenFCS,"RDS/gpr_list_shrunkenFCS.rds")
saveRDS(tpms_TC_gpr,"RDS/tpms_TC_gpr.rds")

