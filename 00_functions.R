### general functions ####
rowMax <- function(frame){
  rowMaxs <- apply(frame,1,max)
  return(rowMaxs)
}

rowMin <- function(frame){
  rowMins <- apply(frame,1,min)
  return(rowMins)
}

ztrans <- function(values,
                   uncenter=T){
  mean_val <- mean(values)
  sd_val <- sd(values)
  
  ztrans <- (values - mean_val) / sd_val
  if(uncenter==T){
    ztrans <- ztrans - ztrans[1]
  }
  return(ztrans)
}


### Functions timing KOs ####
distanceKOvsTP <- function(KO_wt_mat, gene_profile_24){

  distances_KOs <- matrix(NA, nrow=dim(KO_wt_mat)[2], ncol=dim(gene_profile_24)[2])
  row.names(distances_KOs) <- colnames(KO_wt_mat)
  colnames(distances_KOs) <- colnames(gene_profile_24)

  for(d in 1:dim(distances_KOs)[1]){

    row_dist <- numeric()

    for(c in 1:dim(distances_KOs)[2]){

      row_dist <- rbind(row_dist, dist(rbind(KO_wt_mat[,d],gene_profile_24[,c])))

    }

    distances_KOs[d,] <- row_dist

  }

  return(distances_KOs)

}

getRelativeDistance <- function(distances_KOs){

  distances_KOs_rel <- distances_KOs

  for(d in 1:dim(distances_KOs_rel)[1]){

    distances_KOs_rel[d,] <- (distances_KOs_rel[d,]-min(distances_KOs_rel[d,]))/
      (max(distances_KOs_rel[d,])-min(distances_KOs_rel[d,]))

  }

  return(distances_KOs_rel)

}

getGlobalTiming <- function(gene_profiles_24,distances_KOs_diff,
                            row,KO_wt_mat_diff,row_name,
                            duplicates=TRUE, lookup=lookup){

  if(duplicates==TRUE){
    if(row%in%c(1,2,3,4,20,21)){
      if(row==3){#print("test3")
        profiles <- gene_profiles_24[,c((row-2),row,(row+2))]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-2,
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]+2,
                          length.out=17)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9)[1:8],
                                seq(profiles[r,2],profiles[r,3],length.out=9))
        }
      }
      if(row==4){#print("test4")
        profiles <- gene_profiles_24[,c((row-2),row,(row+1))]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-2,
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]+2,
                          length.out=17)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9)[1:8],
                                seq(profiles[r,2],profiles[r,3],length.out=9))
        }
      }
      if(row==1|row==2){#print("test1")
        profiles <- gene_profiles_24[,c(row,(row+2))]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2],
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]+2,
                          length.out=9)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9))
        }
      }

      if(row==20){#print("test3")
        profiles <- gene_profiles_24[,c((row-1),row)]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-1,
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2],
                          length.out=9)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9))
        }
      }
      if(row==21){#print("test3")
        profiles <- gene_profiles_24[,c((row-2),row)]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-2,
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2],
                          length.out=9)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9))
        }
      }



    }

    else{
      profiles <- gene_profiles_24[,c((row-1):(row+1))]

      inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-2,
                        lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]+2,
                        length.out=17)
      ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
      colnames(ext_profiles) <- inter_time
      row.names(ext_profiles) <- row.names(profiles)
      if(dim(profiles)[2]<3){
        print(head(profiles))
      }
      for(r in 1:dim(ext_profiles)[1]){
        ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9)[1:8],
                              seq(profiles[r,2],profiles[r,3],length.out=9))
      }

    }
    #print("test")
    #print(head(ext_profiles))

    distance_row <- distanceKOvsTP(KO_wt_mat_diff,ext_profiles)
    distance_row <- distance_row[which(row.names(distance_row)==row_name),]
    best_time_inter <- inter_time[which(distance_row==min(distance_row))]
    return(list(best_time_inter,min(distance_row),ext_profiles,inter_time,distance_row))
  }
  else{
    if(row%in%c(1,18)){
      if(row==1){#print("test1")
        profiles <- gene_profiles_24[,c(row,(row+1))]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2],
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]+1,
                          length.out=9)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9))
        }
      }

      if(row==18){#print("test3")
        profiles <- gene_profiles_24[,c((row-1),row)]

        inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-1,
                          lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2],
                          length.out=9)

        ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
        colnames(ext_profiles) <- inter_time
        row.names(ext_profiles) <- row.names(profiles)
        #if(dim(profiles)[2]<3){
        #  print(head(profiles))
        #}
        for(r in 1:dim(ext_profiles)[1]){
          ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9))
        }
      }
    }

    else{
      profiles <- gene_profiles_24[,c((row-1):(row+1))]

      inter_time <- seq(lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]-2,
                        lookup[which(lookup$names==colnames(distances_KOs_diff)[row]),2]+2,
                        length.out=17)
      ext_profiles <- as.data.frame(matrix(NA,nrow=dim(profiles)[1],ncol=length(inter_time)))
      colnames(ext_profiles) <- inter_time
      row.names(ext_profiles) <- row.names(profiles)
      if(dim(profiles)[2]<3){
        print(head(profiles))
      }
      for(r in 1:dim(ext_profiles)[1]){
        ext_profiles[r,] <- c(seq(profiles[r,1],profiles[r,2],length.out=9)[1:8],
                              seq(profiles[r,2],profiles[r,3],length.out=9))
      }

    }
    #print("test")
    #print(head(ext_profiles))

    distance_row <- distanceKOvsTP(KO_wt_mat_diff,ext_profiles)
    distance_row <- distance_row[which(row.names(distance_row)==row_name),]
    best_time_inter <- inter_time[which(distance_row==min(distance_row))]
    return(list(best_time_inter,min(distance_row),ext_profiles,inter_time,distance_row))
  }

}
### time induction repression functions ####
get_regulation_from_runs <- function(gene_names_cutoff=character(),
                                     slopes_frame_TC=data.frame(),
                                     log2FC_frame_TC=data.frame(),
                                     thresh_0=0.1){
  ## idee: slopes bestimmen sign aber log2FC für thresh falls wir 0 setzten
  ## log2FC muss vom darauf folgenden Punkt genommen werden damit kann 32 kein induktions
  ## oder repressions Punkt sein

  ## set entrys to 0 that do not make the fat enough from 0 cutoff
  if(length(gene_names_cutoff)>0){
    log2FC_frame_TC_0 <- log2FC_frame_TC[gene_names_cutoff,2:17] ## omit 2i as it does not add info
    slope_sign_frame <- sign(slopes_frame_TC[gene_names_cutoff,])
    for(d in 1:dim(log2FC_frame_TC_0)[1]){
      log2FC_frame_TC_0[d,which(abs(log2FC_frame_TC_0[d,]) <= thresh_0)] <- 0
      slope_sign_frame[d,which(log2FC_frame_TC_0[d,] == 0)-1] <- 0
    }
  }

  ## set slopes to 0 that are connected to 0 entries in log2FCs

  #return(list(log2FC_frame_TC_0,slope_sign_frame))

  induct_repr_results <- list()
  for(d in 1:dim(slope_sign_frame)[1]){
    if(length(slope_sign_frame[d,]!=0)>0){

      holder <- slope_sign_frame[d,]
      names(holder) <- seq(0,30,2)
      holder <- holder[which(holder!=0)]

      keep_first_per_run <- 1
      for(c in 2:length(holder)){
        if(holder[c]==-holder[c-1]){
          keep_first_per_run <- c(keep_first_per_run,c)
        }
      }
      holder <- holder[keep_first_per_run]

      induct_repr_results <- c(induct_repr_results,
                               list(holder))
    }
    else{
      induct_repr_results <- c(induct_repr_results,
                               list(NA))
    }
  }
  names(induct_repr_results) <- row.names(slope_sign_frame)

  return(induct_repr_results)
}

get_induction_repression <- function(second_derivate=data.frame(),
                                     gene_names_cutoff=character(),
                                     slopes_frame_TC=data.frame(),
                                     MinMaxFrame=data.frame(),
                                     MinMaxType="log2FC",
                                     dynamic_thresh_perc=1.2){

  ## function will calculate induction and repression time points for time series based on
  ## slopes and changes of slopes over the time series for a set of genes

  ## The function does miss genes that are strictly monotonic in their slopes if they are changing
  ## more or less constant over the time. You can define genes that behave strictly monotonic
  ## and change over time but are not assigned induction or repression time points
  ##and choose 0 as induction or repression time point for those

  ## second_derivate should be a dataframe that has the changes in slopes of neighbouring time points
  ## as columns and genes in rows

  ## slopes_frame_TC should be a dataframe that has the slopes at given time points
  ## as columns and genes in rows
  ## colnames should have structure: slope_timepoint1_timepoint2

  ## number of colums in second_derivate should be the same as in slopes_frame_TC
  ## an additional column will be added to the beginning of slopes_frame_TC (slope_0_0)
  ## with the assumption that the gene is coming out of a steady state

  ## gene_names_cutoff is a character vector with names of genes that should be analyzed as the
  ## function will identify mainly noise if genes do not change over the time course and induction
  ## and repression points are not usefull

  ## MinMaxFrame is the data frame that was used to calculate slopes in the first place
  ## might be TPMs or log2FCs at the moment (important to define a dynamic threshold per gene)

  ## dynamic_thresh_perc changes the percentage of the number of measurements in calculating the
  ## dynamic threshold. Making it bigger than one will be more sensetive at the cost of
  ## picking up more noisy changes and vice versa.

  ## in this step the 2nd derivative and slope frames are cut to the genes that are defined by
  ## gene_names_cutoff and a slope_0_0 column is added to the beginning of the slopes frame
  holder <- trim_2nd_derivative_and_slopes(second_derivate,
                                           gene_names_cutoff,
                                           slopes_frame_TC)

  short_2nd_deriv <- holder[[1]]
  short_slopes_frame_TC <- holder[[2]]
  MinMaxFrame <- MinMaxFrame[which(row.names(MinMaxFrame)%in%
                                     gene_names_cutoff),]

  ## which columns in the second derivate frame exceed a dynamic threshold and migh be considered
  ## as an induction or repression time point
  cols_2nd_deriv <- list()

  for(d in 1:dim(short_2nd_deriv)[1]){
    gene <- row.names(short_2nd_deriv)[d]

    if(MinMaxType=="TPM"){
      dynamic_thresh <- log2(((max(MinMaxFrame[gene,])+0.01)/(min(MinMaxFrame[gene,])+0.01)))/(dim(short_slopes_frame_TC)[2]*dynamic_thresh_perc)
    }
    else if(MinMaxType=="log2FC"){
      dynamic_thresh <- dist(c(max(MinMaxFrame[gene,]),min(MinMaxFrame[gene,])))/(dim(short_slopes_frame_TC)[2]*dynamic_thresh_perc)
    }

    holder <- short_2nd_deriv[gene,]
    holder <- holder[which(abs(holder) > dynamic_thresh)]
    names(holder) <- which(abs(short_2nd_deriv[gene,]) > dynamic_thresh)

    if(length(holder)>0){
      cols_2nd_deriv <- c(cols_2nd_deriv, list(holder))
    } else {
      cols_2nd_deriv <- c(cols_2nd_deriv, list(NA))
    }

  }
  names(cols_2nd_deriv) <- row.names(short_2nd_deriv)

  ## list of cols that might be induction or repression time points and slope and 2nd
  ## derivative  frames are used to identify which of those points are actually an
  ## induction or repression point
  induction_repression_all_points <- get_induction_repression_points(cols_2nd_deriv,
                                                                     short_2nd_deriv,
                                                                     short_slopes_frame_TC)

  return(induction_repression_all_points)
}

trim_2nd_derivative_and_slopes <- function(second_derivate=data.frame(),
                                           gene_names_cutoff=character(),
                                           slopes_frame_TC=data.frame()){

  ## only selecting genes that were predefined (genes that change over time)
  ## in the 2nd derivative and slope frames
  ## also adding a slope_0_0 column to the slope frame
  short_2nd_deriv <- second_derivate[which(row.names(second_derivate)%in%
                                             gene_names_cutoff),]

  short_slopes_frame_TC <- slopes_frame_TC[which(row.names(slopes_frame_TC)%in%
                                                   gene_names_cutoff),]
  short_slopes_frame_TC <- cbind(matrix(0, nrow = dim(short_slopes_frame_TC)[1],
                                        ncol = 1), short_slopes_frame_TC)
  colnames(short_slopes_frame_TC)[1] <- "slope_0_0"

  return(list(short_2nd_deriv,short_slopes_frame_TC))
}

get_induction_repression_points <- function(cols_2nd_deriv=data.frame(),
                                            short_2nd_deriv=data.frame(),
                                            short_slopes_frame_TC=data.frame()){

  ## function that uses information of slopes and 2nd derivative from the time course as well
  ## as information on which columns surpass a certain change in slopes (2nd derivative).
  ## the output will be a dataframe that is ready to be plotted with ggplot. Each row is
  ## one induction or repression time point which has information on the type (induction or repression),
  ## the Gene (row names of dataframes) and the time point (based on col names of slopes frame)

  induction_repression_all_points <- list()
  for(d in 1:length(cols_2nd_deriv)){
    if(!is.na(cols_2nd_deriv[[d]])[1]){
      holder <- c()
      for(c in 1:length(cols_2nd_deriv[[d]])){
        if(short_slopes_frame_TC[d,(as.numeric(names(cols_2nd_deriv[[d]]))[c]+1)]<0){
          holder <- c(holder, -1)
        }
        else if(short_slopes_frame_TC[d,(as.numeric(names(cols_2nd_deriv[[d]])[c]))+1]>0){
          holder <- c(holder, 1)
        } else {
          holder <- c(holder, 0)
        }
      }
      names(holder) <- sapply(strsplit(colnames(short_slopes_frame_TC)[as.numeric(names(cols_2nd_deriv[[d]]))+1],
                                       split = "_"),function(x) x[2])

      ## we only want to define the first time point to be an induction or repression time point
      ## if it is a series of time points that make the cut off. We also don't want an induction
      ## point followed by an induction point (same for repression)
      omit <- c()
      if(length(holder)>1){
        for(p in 2:length(holder)){
          if(holder[p]==holder[p-1]){
            omit <- c(omit,p)
          }
        }
        holder <- holder[which(!1:length(holder)%in%omit)]
      }

      induction_repression_all_points <- c(induction_repression_all_points,
                                           list(holder))
    } else {
      induction_repression_all_points <- c(induction_repression_all_points,
                                           list(NA))
    }
  }
  names(induction_repression_all_points) <- row.names(short_2nd_deriv)

  return(induction_repression_all_points)
}


convert_results_to_mat <- function(induction_represson_result){
  out_frame <- data.frame(matrix(0,ncol=length
                                 (c("0","2h","4h","6h",
                                    "8h","10h","12h","14h",
                                    "16h","18h","20h","22h",
                                    "24h","26h","28h","30h")),
                                 nrow=length(induction_represson_result)))
  colnames(out_frame) <- c("0","2h","4h","6h",
                           "8h","10h","12h","14h",
                           "16h","18h","20h","22h",
                           "24h","26h","28h","30h")
  row.names(out_frame) <- names(induction_represson_result)

  for(l in 1:length(induction_represson_result)){
    if(!is.na(induction_represson_result[[l]])[1]){
      indexes <- as.numeric(names(induction_represson_result[[l]]))
      out_frame[l,indexes] <- induction_represson_result[[l]]
    }

  }

  return(out_frame)
}

### edge hierarchy for cytoscape ####
get_essential_edges <- function(cytoscape_frame,
                                root_node){

  cytoscape_frame_short <- data.frame()

  current_node <- root_node

  daughter_nodes <- as.character(cytoscape_frame[which(cytoscape_frame$out_node==root_node),2])

  while(length(daughter_nodes)>0){
    daughters_next_lvl <- daughter_nodes[which(daughter_nodes%in%as.character
                                               (cytoscape_frame[which(cytoscape_frame$out_node%in%
                                                                        daughter_nodes),2]))]
    left_daughters <- daughter_nodes[which(!daughter_nodes%in%daughters_next_lvl)]
    print(paste("current node:",current_node))

    cytoscape_frame_short <- rbind(cytoscape_frame_short,
                                   cytoscape_frame[which(cytoscape_frame$out_node==current_node&
                                                           cytoscape_frame$in_node%in%left_daughters),])
    #if(current_node == "cluster_24"){
    #  print(cytoscape_frame[which(cytoscape_frame$out_node==current_node),])
    # print(daughter_nodes)
    #  print(left_daughters)
    #  print(cytoscape_frame_short)
    #}

    for(d in 1:length(left_daughters)){
      print(paste("left daughters:",d))
      cytoscape_frame_short <- rbind(cytoscape_frame_short,
                                     get_essential_edges(cytoscape_frame,
                                                         left_daughters[d]))
    }
    daughter_nodes <- c()

  }
  return(unique(cytoscape_frame_short))
}
### plots ####
plotGprKinetics <- function(gene_name,gpr_list,mgi=TRUE,
                            type="shrunken", unfitted = NA,
                            induction_repression=NA,
                            title="gene", biomart=NA,
                            shrun){

  if(!type%in%c("counts","shrunken","shrunken 6 dropped"))stop("'type' must be one of the following: 'counts','shrunken','shrunken 6 dropped")

  induct_rep_tps <- c()
  induct_rep_groups <- c()
  if(!is.null(dim(induction_repression))){
    if(dim(induction_repression)[[1]]==0){
      stop("empty induction and repression time point frame submitted")
    }else{
      induct_rep_tps <- induction_repression[,2]
      induct_rep_groups <- as.character(induction_repression[,3])
      induct_rep_groups[which(induct_rep_groups=="induction")] <- "red"
      induct_rep_groups[which(induct_rep_groups=="repression")] <- "blue"
    }
  }

  if(title=="gene"){
    if(mgi==TRUE){
      title <- gene_name
    }else{
      title <- biomart[which(biomart$ensembl_gene_id==gene_name),3]
    }

  }

  if(type=="counts"){
    time <- c(0,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,32)
    if(mgi==TRUE){
      ensmug <- biomart[which(biomart$mgi_symbol==gene_name),1]
    }
    else{
      ensmug <- gene_name
      gene_name <- biomart[which(biomart$ensembl_gene_id==ensmug),3]
    }

    if(is.null(unfitted)){
      unfitted <- log_norm_counts
    }
    orig <- as.numeric(unfitted[which(row.names(unfitted)==ensmug),c(3:21)])
    smoothed <- as.numeric(gpr_list[[1]][which(row.names(gpr_list[[1]])==ensmug),c(3:21)])

    q1 <- as.numeric(gpr_list[[2]][which(row.names(gpr_list[[2]])==ensmug),])
    q2 <- as.numeric(gpr_list[[3]][which(row.names(gpr_list[[3]])==ensmug),])

    #fit_mean <- loess(smoothed~time)
    points <- data.frame(x=time,y=orig,group="points")
    fit_mean <- data.frame(x=time,y=smoothed,group="mean")
    fit_95 <- data.frame(x=time,y=q1,group="q1")
    fit_5 <- data.frame(x=time,y=q2,group="q2")

    ggframe <- rbind(points, fit_mean,
                     fit_95, fit_5)
    #print(head(subset(ggframe,group== "points")))

    plot <- ggplot(ggframe,aes(x=x, y=y))+
      geom_point(data=subset(ggframe,group== "points"),
                 aes(x=x, y=y), shape=21, size=2) +
      geom_line(data=subset(ggframe,group== "mean"),
                aes(x=x, y=y))+
      geom_line(data=subset(ggframe,group== "q1"),
                aes(x=x, y=y),col="red")+
      geom_line(data=subset(ggframe,group== "q2"),
                aes(x=x, y=y),col="red")+
      geom_hline(yintercept = mean(smoothed),lty=2,col="grey 70")+
      xlab("time in hrs")+
      ylab("log2fc expression vs 2i expression")+
      ggtitle(title)+
      geom_vline(xintercept = induct_rep_tps,
                 col = as.factor(induct_rep_groups),
                 lty = 2, alpha = 0.5)+
      scale_x_continuous(breaks = seq(0,32,4))+
      theme(plot.title=element_text(hjust=0.5),
            panel.background = element_rect(fill=NA),
            panel.grid.major = element_line(colour = "grey90"))
    return(plot)

    #plot(time,orig,xlab="time in hrs",ylab="log2fc expression vs 2i expression",main=gene_name,xaxt="n")
    #abline(h=mean(smoothed),lty=2,col="blue")
    #lines(fit_mean)
    #lines(fit_95,col="red")
    #lines(fit_5,col="red")
    #axis(side=1,at=seq(0,32,4))

  }

  if(type=="shrunken"){
    time <- c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)
    if(mgi==TRUE){
      ensmug <- biomart[which(biomart$mgi_symbol==gene_name),1]
    }
    else{
      ensmug <- gene_name
      gene_name <- biomart[which(biomart$ensembl_gene_id==ensmug),3]
    }

    if(is.null(unfitted)){
      unfitted <- shrunken_foldchanges_complete
    }
    orig <- as.numeric(unfitted[which(row.names(unfitted)==ensmug),c(2:18)])
    smoothed <- as.numeric(gpr_list[[1]][which(row.names(gpr_list[[1]])==ensmug),c(2:18)])

    q1 <- as.numeric(gpr_list[[2]][which(row.names(gpr_list[[2]])==ensmug),])
    q2 <- as.numeric(gpr_list[[3]][which(row.names(gpr_list[[3]])==ensmug),])

    #plot(time,orig,xlab="time in hrs",ylab="log2fc expression vs 2i expression",main=gene_name,xaxt="n")
    #fit_mean <- loess(smoothed~time)
    points <- data.frame(x=time,y=orig,group="points")
    fit_mean <- data.frame(x=time,y=smoothed,group="mean")
    fit_95 <- data.frame(x=time,y=q1,group="q1")
    fit_5 <- data.frame(x=time,y=q2,group="q2")


    ggframe <- rbind(points, fit_mean,
                     fit_95, fit_5)
    #print(head(subset(ggframe,group== "points")))

    plot <- ggplot(ggframe,aes(x=x, y=y))+
      geom_point(data=subset(ggframe,group== "points"),
                 aes(x=x, y=y), shape=21, size=2) +
      geom_line(data=subset(ggframe,group== "mean"),
                aes(x=x, y=y))+
      geom_line(data=subset(ggframe,group== "q1"),
                aes(x=x, y=y),col="red")+
      geom_line(data=subset(ggframe,group== "q2"),
                aes(x=x, y=y),col="red")+
      geom_hline(yintercept = mean(smoothed),lty=2,col="grey 70")+
      xlab("time in hrs")+
      ylab("log2fc expression vs 2i expression")+
      ggtitle(title)+
      geom_vline(xintercept = induct_rep_tps,
                 col = as.factor(induct_rep_groups),
                 lty = 2, alpha = 0.5)+
      scale_x_continuous(breaks = seq(0,32,4))+
      theme(plot.title=element_text(hjust=0.5),
            panel.background = element_rect(fill=NA),
            panel.grid.major = element_line(colour = "grey90"))
    return(plot)
    #abline(h=mean(smoothed),lty=2,col="blue")
    #lines(fit_mean)
    #lines(fit_95,col="red")
    #lines(fit_5,col="red")
    #axis(side=1,at=seq(0,32,4))
  }

  if(type=="shrunken 6 dropped"){
    time <- c(0,2,4,8,10,12,14,16,18,20,22,24,26,28,30,32)
    if(mgi==TRUE){
      ensmug <- biomart[which(biomart$mgi_symbol==gene_name),1]
    }
    else{
      ensmug <- gene_name
      gene_name <- biomart[which(biomart$ensembl_gene_id==ensmug),3]
    }

    if(is.null(unfitted)){
      unfitted <- shrunken_foldchanges_complete
    }
    orig <- as.numeric(shrunken_foldchanges_complete[which(row.names(shrunken_foldchanges_complete)==ensmug),c(2:4,6:18)])
    print(orig)
    smoothed <- as.numeric(gpr_list[[1]][which(row.names(gpr_list[[1]])==ensmug),c(2:17)])
    print(smoothed)
    q1 <- as.numeric(gpr_list[[2]][which(row.names(gpr_list[[2]])==ensmug),])
    q2 <- as.numeric(gpr_list[[3]][which(row.names(gpr_list[[3]])==ensmug),])

    #plot(time,orig,xlab="time in hrs",ylab="log2fc expression vs 2i expression",main=gene_name,xaxt="n")
    #fit_mean <- loess(smoothed~time)
    points <- data.frame(x=time,y=orig,group="points")
    fit_mean <- data.frame(x=time,y=smoothed,group="mean")
    fit_95 <- data.frame(x=time,y=q1,group="q1")
    fit_5 <- data.frame(x=time,y=q2,group="q2")


    ggframe <- rbind(points, fit_mean,
                     fit_95, fit_5)
    #print(head(subset(ggframe,group== "points")))

    plot <- ggplot(ggframe,aes(x=x, y=y))+
      geom_point(data=subset(ggframe,group== "points"),
                 aes(x=x, y=y), shape=21, size=2) +
      geom_line(data=subset(ggframe,group== "mean"),
                aes(x=x, y=y))+
      geom_line(data=subset(ggframe,group== "q1"),
                aes(x=x, y=y),col="red")+
      geom_line(data=subset(ggframe,group== "q2"),
                aes(x=x, y=y),col="red")+
      geom_hline(yintercept = mean(smoothed),lty=2,col="grey 70")+
      xlab("time in hrs")+
      ylab("log2fc expression vs 2i expression")+
      ggtitle(title)+
      geom_vline(xintercept = induct_rep_tps,
                 col = as.factor(induct_rep_groups),
                 lty = 2, alpha = 0.5)+
      scale_x_continuous(breaks = seq(0,32,4))+
      theme(plot.title=element_text(hjust=0.5),
            panel.background = element_rect(fill=NA),
            panel.grid.major = element_line(colour = "grey90"))
    return(plot)
    #abline(h=mean(smoothed),lty=2,col="blue")
    #lines(fit_mean)
    #lines(fit_95,col="red")
    #lines(fit_5,col="red")
    #axis(side=1,at=seq(0,32,4))
  }

}


### comparing SC to TC ####
gene_gene_dependencies <- function(Gene_A, Gene_B,
                                   SC_data,
                                   gpr_TPMs_short,
                                   NA_cut=numeric(),
                                   early = numeric(),
                                   late = numeric(),
                                   return_dependency=F,
                                   return_extreme_mean=F,
                                   extreme_cut=0.2,
                                   counts=F,
                                   buffer_genes=500){

  holder_A <- SC_data[Gene_A,]
  holder_B <- SC_data[Gene_B,]

  signs <- c(gpr_TPMs_short[Gene_A,21],
             gpr_TPMs_short[Gene_B,21])
  time_points <- sapply(colnames(SC_data), function(x) strsplit(x,"_")[[1]][1])

  if(length(NA_cut)==1){
    early_A <- early[which(SC_data[Gene_A,c(early)]>NA_cut)]
    late_A <- late[which(SC_data[Gene_A,c(late)]>NA_cut)]

    early_B <- early[which(SC_data[Gene_B,c(early)]>NA_cut)]
    late_B <- late[which(SC_data[Gene_B,c(late)]>NA_cut)]

    early <- early_A[which(early_A%in%early_B)]
    late <- late_A[which(late_A%in%late_B)]

    cells_non_NA <- which(holder_A>NA_cut & holder_B>NA_cut)
    holder_A <- holder_A[cells_non_NA]
    holder_B <- holder_B[cells_non_NA]
    SC_data <- SC_data[,cells_non_NA]
  }

  if(length(holder_A) < 40){
    warning("not enough cells in either early and late time point")
    return(NULL)
  }
  if(counts == F){
    A_scaled <- scale(as.numeric(holder_A), center = T, scale = T)
    B_scaled <- scale(as.numeric(holder_B), center = T, scale = T)

    A_scaled_sign <- signs[1]*A_scaled

    B_scaled_sign <- signs[2]*B_scaled
  }else{
    ## no scaling on counts
    A_scaled <- as.numeric(holder_A)

    B_scaled <- as.numeric(holder_B)

    A_scaled_sign <- signs[1]*A_scaled

    B_scaled_sign <- signs[2]*B_scaled
  }


  if(return_dependency == F){
    return(list(as.numeric(A_scaled_sign),
                as.numeric(B_scaled_sign),
                time_points,
                colMeans(rbind(holder_A,holder_B))))
  }else{
    rank_A <- rank(A_scaled_sign, ties.method = "first")
    rank_B <- rank(B_scaled_sign, ties.method = "first")
    rank_diff <- rank_A - rank_B

    if(counts == F){
      if(!is.na(wilcox.test(rank_diff)[3])){
      wilc_test <- wilcox.test(rank_diff)

      if(return_extreme_mean == T){
        rank_diff_extr <- rank_diff[order(abs(rank_diff),decreasing = T)]
        rank_diff_extr <- rank_diff_extr[1:round(length(rank_diff_extr)*extreme_cut,1)]

        return(list(mean=mean(rank_diff_extr)/length(rank_A),
                    pval=wilc_test$p.value))
        }else{
          return(list(median=-median(rank_diff)/length(rank_A),
                      pval=wilc_test$p.value))
        }
      }
    }else{
      if(signs[1]==1){
        max_A <- max(abs(rank_diff[which(holder_A==0.1 &
                                           holder_B==0.1 &
                                           time_points=="24h")]))
      }else{
        max_A <- max(abs(rank_diff[which(holder_A==0.1 &
                                           holder_B==0.1 &
                                           time_points=="0h")]))
      }
      if(signs[2]==1){
        max_B <- max(abs(rank_diff[which(holder_B==0.1 &
                                           holder_A==0.1 &
                                           time_points=="24h")]))
      }else{
        max_B <- max(abs(rank_diff[which(holder_B==0.1 &
                                           holder_A==0.1 &
                                           time_points=="0h")]))
      }

      buffer <- max(max_A,max_B)
      buffer <- buffer * 1.1

      rank_diff_buff <- rank_diff[which(abs(rank_diff)>buffer)]

      if(length(rank_diff_buff)>=buffer_genes){
        binom_p <- binom.test(x = sum(rank_diff_buff<0),
                            n = length(rank_diff_buff),
                            p = 0.50,
                            alternative = "two.sided")$p.value

        rank_diff <- rank_diff_buff

        if(return_extreme_mean == T){
          rank_diff_extr <- rank_diff[order(abs(rank_diff),decreasing = T)]
          rank_diff_extr <- rank_diff_extr[1:round(length(rank_diff_extr)*extreme_cut,1)]

          return(list(mean=mean(rank_diff_extr)/length(rank_A),
                      pval=binom_p))
          }else{
            return(list(mean=mean(rank_diff)/length(rank_A),
                        pval=binom_p))
          }
        }

      }
   }
}


get_relations_from_SC <- function(SC_relations,
                                  SC_data,
                                  gpr_TPMs_short,
                                  NA_cut=numeric(),
                                  early = numeric(),
                                  late = numeric(),
                                  return_dependency=F,
                                  return_extreme_mean=F,
                                  extreme_cut=0.2,
                                  contingency_out=F,
                                  early_late_cont=F,
                                  padj_cut=0.05,
                                  counts=F,
                                  buffer_genes=500){

  index_matrix <- matrix(NA,nrow=nrow(SC_relations),
                                ncol=ncol(SC_relations))
  row.names(index_matrix) <- row.names(SC_relations)
  colnames(index_matrix) <- row.names(SC_relations)

  index <- 1

  ## direction of relation between genes
  A_relation_B <- list()
  ## mean relation between genes
  A_mean_B <- list()
  ## adjusted pval is the distribution A-B signigicantly different from 0
  A_padj_B <- list()
  ## pval is the distribution A minus B signigicantly different from 0
  A_pval_B <- list()
  ## contingency of 0 count cells
  contingency <- list()
  ## contingency of 0 count cells
  contingency_early_late <- list()


  ## run through each field in the possible dependency matrix
  for(r in 1:nrow(index_matrix)){
    print(r)
    for(c in 1:ncol(index_matrix)){
      if(r>c){
        if(SC_relations[r,c]!=0){
          holder <- gene_gene_dependencies(Gene_A = row.names(index_matrix)[r],
                                           Gene_B = colnames(index_matrix)[c],
                                           SC_data = SC_data,
                                           gpr_TPMs_short=gpr_TPMs_short,
                                           NA_cut = NA_cut,
                                           early = T0_cells, late = T24_cells,
                                           return_dependency = return_dependency,
                                           return_extreme_mean=return_extreme_mean,
                                           extreme_cut=extreme_cut,
                                           counts=counts,
                                           buffer_genes=buffer_genes)
          if(contingency_out==T){
            holder_cont <- gene_gene_zero_contingency(Gene_A = row.names(index_matrix)[r],
                                                    Gene_B = colnames(index_matrix)[c],
                                                    SC_data = SC_data,
                                                    NA_cut = NA_cut,
                                                    early_late=early_late_cont,
                                                    gpr_TPMs_short=gpr_TPMs_short)
            }
          }else{
            holder <- NA
            }

          if(is.numeric(holder[[1]])){
            A_mean_B <- c(A_mean_B,
                                       holder[[1]])
            A_pval_B <- c(A_pval_B,
                                 holder[[2]])
            if(contingency_out==T){
              contingency <- c(contingency,
                                            list(holder_cont))
            }



            index_matrix[r,c] <- index
            index <- index + 1
          }
        }
      }
    }


  A_padj_B <- lapply(A_pval_B, FUN = function(x)
    p.adjust(x,"BH",length(A_pval_B)))

  A_relation_B <- mapply(x=A_padj_B, y=A_mean_B, SIMPLIFY = F,
                         FUN = function(x, y)
    if(x < padj_cut){y}else{0})

  if(contingency_out==T){
    return(list(index_matrix,
              A_relation_B,
              A_mean_B,
              A_padj_B,
              A_pval_B,
              contingency))
  }else{
    return(list(index_matrix,
                A_relation_B,
                A_mean_B,
                A_padj_B,
                A_pval_B))
  }

}


### relations based on percentage ####

## function to calculate relation based on percentage of WT table
## the percentage table carries the percentage of WT completion
## of different processes in all KOs

get_relations_from_percentage_genewise <- function(percentage_table,
                                                   sig_diff_0=TRUE,
                                                   sig_diff_0_p=0.05){
  ## construct tables to be filled
  ## pval
  A_pval_B_matrix <- matrix(NA,nrow=dim(percentage_table)[1],
                            ncol=dim(percentage_table)[1])
  row.names(A_pval_B_matrix) <- row.names(percentage_table)
  colnames(A_pval_B_matrix) <- row.names(percentage_table)
  ## mean
  A_mean_B_matrix <- matrix(NA,nrow=dim(percentage_table)[1],
                            ncol=dim(percentage_table)[1])
  row.names(A_mean_B_matrix) <- row.names(percentage_table)
  colnames(A_mean_B_matrix) <- row.names(percentage_table)
  ## relation
  A_relation_B_matrix <- matrix(NA,nrow=dim(percentage_table)[1],
                                ncol=dim(percentage_table)[1])
  row.names(A_relation_B_matrix) <- row.names(percentage_table)
  colnames(A_relation_B_matrix) <- row.names(percentage_table)
  ## consistency
  #consistency_frame <- data.frame(A=character(),
  #                               B=character(),
  #                              relation=numeric(),
  #                             mean=numeric())

  process_mat <- lapply(1:as.numeric(nrow(percentage_table)^2),
                        function(x) matrix(NA, ncol=2,
                                           nrow=ncol(percentage_table)))
  relation_mat <- lapply(1:as.numeric(nrow(percentage_table)^2),
                         function(x) matrix(NA, ncol=2,
                                            nrow=ncol(percentage_table)))

  percentage_mat <- lapply(1:as.numeric(nrow(percentage_table)^2),
                           function(x) matrix(NA, ncol=2,
                                              nrow=1))
  pair_vec <- lapply(1:as.numeric(nrow(percentage_table)^2),
                           function(x) as.character(NA))

  for(r in 1:dim(A_pval_B_matrix)[1]){
    if(r%%50==0){
      print(r)
    }
    for(c in 1:dim(A_pval_B_matrix)[1]){
      Gene_A <- percentage_table[r,]
      Gene_B <- percentage_table[c,]
      diff_A_B <- as.numeric(Gene_A-Gene_B)
      ## all KOs are used, previously this was not the case
      #index_KOs <- 1:dim(percentage_table)[2]
      ## if the difference between two processes in all KOs is 0 fill the matrices
      ## with corresponding values
      if(length(which((diff_A_B)!=0))==0){
        A_pval_B_matrix[r,c] <- 1
        A_mean_B_matrix[r,c] <- 0
        A_relation_B_matrix[r,c] <- 0
      }
      ## else calculate mean relation between the two processes and a pvalue
      else{
        A_mean_B_matrix[r,c] <- mean(diff_A_B)
        A_pval_B_matrix[r,c] <- t.test(-(diff_A_B),mu = 0,alternative = "two.sided")$p.value

        ## add consistency frame entries to the output frame
        process_mat[[as.numeric(((r-1)*nrow(A_mean_B_matrix)+c))]] <-
          cbind(rep(row.names(A_mean_B_matrix)[r],
                    ncol(percentage_table)),
                rep(colnames(A_mean_B_matrix)[c],
                    ncol(percentage_table)))
        relation_mat[[as.numeric(((r-1)*nrow(A_mean_B_matrix)+c))]] <-
          cbind(diff_A_B,
                rep(mean(diff_A_B),
                    ncol(percentage_table)))

        pair_vec[[as.numeric(((r-1)*nrow(A_mean_B_matrix)+c))]] <-
          paste(row.names(A_mean_B_matrix)[r],
                colnames(A_mean_B_matrix)[c],
                sep=";")
        percentage_mat[[as.numeric(((r-1)*nrow(A_mean_B_matrix)+c))]] <-
          cbind(sum(sign(diff_A_B)==sign(mean(diff_A_B)))/length(diff_A_B),
          mean(diff_A_B))

        #consistency_frame <- rbind(consistency_frame,
        #                          data.frame(A=row.names(A_mean_B_matrix)[r],
        #                                    B=colnames(A_mean_B_matrix)[c],
        #                                   relation=diff_A_B,
        #                                  mean=mean(diff_A_B)))

      }
    }
  }
  process_mat <- do.call(rbind, process_mat)
  relation_mat <- do.call(rbind, relation_mat)
  consistency_frame <- cbind(as.data.frame(process_mat),
                             as.data.frame(relation_mat))
  colnames(consistency_frame) <- c("A", "B", "relation", "mean")
  consistency_frame <- na.omit(consistency_frame)
  row.names(consistency_frame) <- 1:nrow(consistency_frame)

  percentage_mat <- do.call(rbind, percentage_mat)
  pair_vec <- do.call(rbind, pair_vec)
  summary_consistency <- cbind(as.data.frame(percentage_mat),
                             as.data.frame(pair_vec))
  colnames(summary_consistency) <- c("fraction_same_side", "mean", "pair")
  summary_consistency <- na.omit(summary_consistency)
  row.names(summary_consistency) <- 1:nrow(summary_consistency)

  ##adjust pvalues of previous results
  A_padj_B_matrix <- A_pval_B_matrix
  for(d in 1:dim(A_padj_B_matrix)[1]){
    A_padj_B_matrix[d,] <- p.adjust(A_padj_B_matrix[d,],"BH",length(A_padj_B_matrix[d,]))
  }
  ## if sorting out non significant differences between two processes
  ## calculate a direction for the relationship matrix carrying only directionality
  if(sig_diff_0==TRUE){
    for(r in 1:dim(A_relation_B_matrix)[1]){
      for(c in 1:dim(A_relation_B_matrix)[2]){
        if(A_pval_B_matrix[r,c] < sig_diff_0_p){
          if(A_mean_B_matrix[r,c] < 0){
            A_relation_B_matrix[r,c] <- -1
          }
          if(A_mean_B_matrix[r,c] > 0){
            A_relation_B_matrix[r,c] <- 1
          }
        }
        else{
          A_relation_B_matrix[r,c] <- 0
        }
      }
    }
  }


  return(list(relation_A_B=A_relation_B_matrix,
              mean_A_B=A_mean_B_matrix,
              padj_A_B=A_padj_B_matrix,
              pval_A_B=A_pval_B_matrix,
              consistency=consistency_frame,
              summary_consist=summary_consistency))
}


percentage_plot <- function(percentage_table,
                            Gene_A=character(),
                            Gene_B=character(),
                            nbins=20,
                            plot_out=TRUE){

  weak_pheno <- c("Msi2", "Rps6ka1", "Arih2", "Cabin1",
                  "Hprt", "Tet1", "Igf2bp1", "Nes",
                  "Ssr2", "Irak3", "Etl4", "Pum1",
                  "Myc", "Dido1", "Hnrnph1")

  diff_A_B <- percentage_table[Gene_A,]-percentage_table[Gene_B,]

  hist_obj <- hist(diff_A_B, breaks = nbins, plot = FALSE)

  hist_table <- data.frame(mids=hist_obj$mids,
                           counts=hist_obj$counts,
                           fraction_strong=NA)

  for(d in 1:nrow(hist_table)){
    if(hist_table[d,2] != 0){
      hist_table[d,3] <- 1-(sum(diff_A_B[weak_pheno]>=hist_obj$breaks[d] &
                                  diff_A_B[weak_pheno]<hist_obj$breaks[d+1])/hist_table[d,2])
    }
  }

  ggpl <- ggplot(hist_table,
                 aes(x=mids ,y=counts, fill=fraction_strong))+
    scale_fill_gradient(low="white",high="red", name="percentage strong phenotypes") +
    geom_bar(stat="identity", col="black",
             width=dist(c(hist_obj$breaks[1],hist_obj$breaks[2]))) +
    ggtitle(paste(Gene_A, "-", Gene_B, sep = " ")) +
    xlab(paste(Gene_A, "-", Gene_B, sep = " ")) +
    ylab("Frequency")+
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5))
  if(plot_out){
    plot(ggpl)
  }else{
    return(hist_table)
  }

}

pheno_hist_test <- function(percentage_table,
                            dependency_mat){

  weak_pheno <- c("Msi2", "Rps6ka1", "Arih2", "Cabin1",
                  "Hprt", "Tet1", "Igf2bp1", "Nes",
                  "Ssr2", "Irak3", "Etl4", "Pum1",
                  "Myc", "Dido1", "Hnrnph1")

  dependencies <- sum(dependency_mat!=0)/2

  out_frame <- data.frame(position=rep(NA,dependencies*ncol(percentage_table)),
                          KO=rep(NA,dependencies*ncol(percentage_table)),
                          type=rep(NA,dependencies*ncol(percentage_table)),
                          dependency=rep(NA,dependencies*ncol(percentage_table)))

  range <- c(1,ncol(percentage_table))

  for(r in 1:nrow(dependency_mat)){
    if(r%%10==0){
      print(r)
    }
    for(c in 1:ncol(dependency_mat)){
      if(dependency_mat[r,c] > 0){

        diff_A_B <- percentage_table[row.names(dependency_mat)[r],]-
          percentage_table[colnames(dependency_mat)[c],]

        diff_A_B_mean <- diff_A_B-mean(diff_A_B)
        diff_A_B_mean <- diff_A_B_mean/max(abs(diff_A_B_mean))

        out_frame[range[1]:range[2], 1] <- diff_A_B_mean
        out_frame[range[1]:range[2], 2] <- names(diff_A_B_mean)
        out_frame[range[1]:range[2], 3] <- as.numeric(names(diff_A_B_mean)%in%
                                                        weak_pheno)
        out_frame[range[1]:range[2], 4] <- paste(row.names(dependency_mat)[r],
                                                 colnames(dependency_mat)[c],
                                                 collapse = "_")


        ## update indexes for next non zero dependency
        range[1] <- range[1]+ncol(percentage_table)
        range[2] <- range[2]+ncol(percentage_table)
      }
    }
  }
  return(out_frame)

}


