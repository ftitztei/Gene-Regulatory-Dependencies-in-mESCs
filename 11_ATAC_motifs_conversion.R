library(openxlsx)

setwd("/cellfile/datapublic/ftitztei/Rscripts/TCpackage/")

TF_table <- read.table("/cellfile/datapublic/ftitztei/Luis_ATAC/CisBP/TF_Information.txt",
                       sep = "\t", header = T)

bases <- c("A:", "C:", "G:", "T:")


TF_table_short <- TF_table
length(unique(TF_table_short$TF_Name))

for(d in 1:dim(TF_table_short)[1]){
  if(d>1){
    # add empty line
    write("","data/uniprobe_TF_Motifs.txt",append = T)
  }
  if(TF_table_short$MSource_Identifier[d]!="Transfac" & TF_table_short$Motif_ID[d]!="."){
    # add header for matrix
    write(paste(TF_table_short$TF_Name[d],
                TF_table_short$Motif_ID[d],
                sep="_"),
          "data/uniprobe_TF_Motifs.txt",append = T)

    # read in matrix and parse to uniprobe format
    matrix <- read.table(paste("/cellfile/datapublic/ftitztei/Luis_ATAC/CisBP/pwms_all_motifs/",
                               TF_table_short$Motif_ID[d],
                               ".txt",sep=""))

    # add matrix
    ## add each base seperately
    for(c in 1:4){
      write(paste(bases[c],paste(as.character(matrix[2:length(matrix[,c+1]),c+1]),
                                 collapse="\t"),collapse="\t"),
            "data/uniprobe_TF_Motifs.txt",append = T)
    }


  }

}

length(unique(TF_table_short$TF_Name[which(TF_table_short$MSource_Identifier!="Transfac" &
                                             TF_table_short$Motif_ID!=".")]))

length(unique(TF_table_short$Motif_ID[which(TF_table_short$MSource_Identifier!="Transfac" &
                                              TF_table_short$Motif_ID!=".")]))
