vectorSubset = function(vec, mat){
    # This function is from Shila Ghazanfar as used in the spatial organogenesis paper (Nat. Biotech 2021)
  # vec is a named vector
  # mat is a matrix containing the names or indices for which you want
  # to get the entries of vec
  
  vmat = c(mat)
  vvec = vec[vmat]
  
  vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
  colnames(vecmat) <- colnames(mat)
  rownames(vecmat) <- rownames(mat)
  
  return(vecmat)
}
 
vectorMatch = function(vec, mat, vecnames){
    # This function is from Shila Ghazanfar as used in the spatial organogenesis paper (Nat. Biotech 2021)
  # vec is an unnamed vector
  # vecnames is the names of vec
  # mat is a matrix containing the names or indices for which you want
  # to get the entries of vec, matching vecnames
  
  vmat = c(mat)
  
  vecind = match(vmat,vecnames)
  
  vvec = vec[vecind]
  
  vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
  colnames(vecmat) <- colnames(mat)
  rownames(vecmat) <- rownames(mat)
  
  return(vecmat)
}
 
getmode <- function(v, dist) {
    # This function is from Shila Ghazanfar as used in the spatial organogenesis paper (Nat. Biotech 2021)
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}
 
 
get_MVS_vals = function(predmat) {
    # This function is from Shila Ghazanfar as used in the spatial organogenesis paper (Nat. Biotech 2021)
  MVS_vals_svm = t(apply(predmat,1,function(x) {
    # first column is the mapped celltype
    # second column is the runner-up celltype
    # mode_prop is the proportion of nearest neighbours with either mode1 or mode2
    # d is the nodepath distance for those celltypes
    
    if (all(is.na(x))) return(c(NA,NA,NA,NA))
    
    mode1 = getmode(x, 1:length(x))
    xx = x[x != mode1]
    if (length(xx) == 0) {
      mode2 = NA
      
    } else {
      mode2 = getmode(xx, 1:length(xx))
      if (is.null(mode2)) {
        mode2 <- NA
      }
      
    }
    return(c(mode1, mode2, mean(x %in% mode1), mean(x %in% c(mode1,mode2))))
  }))
  colnames(MVS_vals_svm) <- c("celltype_mapped",
                              "celltype_alternative",
                              "mapping_score",
                              "mapping_alternative_score"
  )
  return(MVS_vals_svm)
  
}
bubblePlot <- function(m, markers, grps, cluster_col=TRUE, cluster_row=TRUE, angled=TRUE, zscores=FALSE) { 
    out <- data.frame(numeric(length(markers)))
    colnames(out) <- levels(factor(grps))[1]
    out.freq <- out
    for (cond in levels(factor(grps))) {
	expr <- rowMeans(m[markers,grps==cond])
	expr.freq <- rowMeans(m[markers,grps==cond]>0)
	colname <- cond
	out[,colname] <- expr
	out.freq[,colname] <- expr.freq
    }
    out <- as.matrix(out)
    out.freq <- as.matrix(out.freq)
    rownames(out.freq) <- rownames(out) <- markers
    if(zscores) {
	out <- t(scale(t(out))) 
    } else {
	out <- out/rowMax(out)
    }
    out.long <- reshape2::melt(out,value.name="Mean")
    out.freq <- reshape2::melt(out.freq,value.name="Frequency")
    out.long$Frequency <- out.freq$Frequency * 100

    if (cluster_row) {
    dis <- dist(out)
    hclst <- hclust(dis,method="ward.D2")
    lvlsVar1 <- rownames(out)[hclst$order]
    out.long$Var1 <- factor(out.long$Var1, levels=lvlsVar1)
    }
    if (cluster_col) {
    dis <- dist(t(out))
    hclst <- hclust(dis,method="ward.D2")
    lvlsVar2 <- colnames(out)[hclst$order]
    out.long$Var2 <- factor(out.long$Var2, levels=lvlsVar2)
    }
    p <- ggplot(out.long, aes(x=Var2, y=Var1, color=Mean, size=Frequency)) +
    geom_point() +
    scale_color_gradient2(low="grey90",high="black") + ##,mid="orange") +
   # scale_color_distiller(palette="Spectral") +
    scale_size(range=c(0,3)) +
    theme(panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
	  axis.text=element_text(size=7),
	  legend.position="bottom",
	  legend.direction="horizontal") +
    xlab("") +
    ylab("")
    if (angled) {
    p <- p + theme(axis.text.x=element_text(size=9, angle=45, hjust=1, vjust=1))
    }
    return(p)
}
#Cluster stuff
mergeCluster <- function(x, clusters, min.DE=10, maxRep=20, removeGenes=NULL, merge=TRUE, ...)
{
    # This function iteratively merges clusters with a number of DE genes less or equal than min.DE
    # Clusters is a factor corresponding with a length equal to ncol(x)
    # x is the log-transformed normalized gene expression matrix (genes x cells)
    # DE is defined as lfc>1 and log FDR <= -2 
    require(scran)
    counter <- c(1:maxRep)
    # if uninteresting features are to be removed
    if (!is.null(removeGenes)) {
	x <- x[!(rownames(x) %in% removeGenes),]
    }
    # Iterate through computing DE genes and merging clusters (if merge=TRUE)
    for (i in counter) {
	# Compute DE Genes
	clust.vals <- levels(clusters)
	marker.list <- findMarkers(x,clusters,subset.row=rowMeans(x)>0.01, lfc=1, full.stats=TRUE, ...)

	## Create Cluster x Cluster matrix with each entry corresponding to number of DE genes (log.FDR <=-2)
	# Setup matrix
	out <- matrix(nrow=length(clust.vals), ncol=length(clust.vals))
	colnames(out) <- clust.vals
	rownames(out) <- clust.vals
	# loop through each cluster and get DE genes against every other cluster
	for (cl in names(marker.list)) {
	    sb <- marker.list[[cl]]
	    sb <- data.frame(sb)
	    sb <- as.data.frame(sb[,grepl(".log.FDR",colnames(sb))]) # as.data.frame in case it's a single column
	    if(ncol(sb)==1) {
		colnames(sb) <- setdiff(names(marker.list),cl)
	    } else {
	    colnames(sb) <- gsub("stats.|.log.FDR","",colnames(sb))
	    }
	    sb <- colSums(as.matrix(sb) <= -3) # -3 is the log FDR threshold for DE
	    sb[cl] <- 0
	    sb <- sb[match(colnames(out),names(sb))]
	    out[,cl] <- sb
	}

	if (merge) {

	# If a test does not have enough d.f. then findMarkers will output NA
	# This usually happens if you use say a blocking factor and 
	# some clusters only appear in some level of that factor
	# In that case I don't have information to merge them
	if (any(is.na(out))) { 
	    print("A test did not have enough D.F.") 
	    out[is.na(out)] <- min.DE+1 
	}
	#Compute new groups
	if (any(is.na(out))) { print("A test did not have enough D.F.") }
	out[is.na(out)] <- min.DE+1 # This is to prevent hclust from crashing, NAs are produced when their are not enough d.f. for testing. In this case I do not wnat the clusters to be merged because I formally cant test for DE.
	hc <- hclust(as.dist(out))
	#get min number of DE genes
	min.height <- min(hc$height)
	if (min.height <= min.DE) {
	    print(paste0(i,". joining of clusters with ",min.height," DE genes"))
	    cs <- cutree(hc,h=min.height)
	    newgrps <- sapply(c(1:max(cs)), function(i) paste(names(cs)[cs==i],collapse="."))
	    csnew <- plyr::mapvalues(cs,c(1:max(cs)), newgrps)
	    clusters <- plyr::mapvalues(clusters,names(csnew),csnew)
	} else {
	    break
	}
    } else {
	break}

	if (length(unique(clusters))==1) {
	    break
	}
    }
	output <- list("Marker"=marker.list,
		       "NumberDE"=out,
		       "NewCluster"=clusters)
	return(output)
}

pointsAndLabels <- function(sce, DimRed, Labels,sz=1,label=FALSE,colour=TRUE) {
    c25 <- c(
	      "dodgerblue2", "#E31A1C", # red
	      "green4",
	      "#6A3D9A", # purple
	      "#FF7F00", # orange
	      "black", "gold1",
	      "skyblue2", "#FB9A99", # lt pink
	      "palegreen2",
	      "#CAB2D6", # lt purple
	      "#FDBF6F", # lt orange
	      "gray70", "khaki2",
	      "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
	      "darkturquoise", "green1", "yellow4", "yellow3",
	      "darkorange4", "brown"
	    )
    require(ggplot2)
    require(scater)
    require(ggrepel)
    pD <- data.frame(colData(sce))
    dr <- reducedDim(sce,DimRed)
    pD[,paste0(DimRed,".1")] <- dr[,1]
    pD[,paste0(DimRed,".2")] <- dr[,2]
    X <- model.matrix(~0+pD[,Labels])
    colnames(X) <- levels(factor(pD[,Labels]))
    X <- t(t(X)/colSums(X))
    celltyp_ctrs <- data.frame(crossprod(X, as.matrix(pD[,c(paste0(DimRed,".1"),paste0(DimRed,".2"))])))
    celltyp_ctrs$Lbls <- rownames(celltyp_ctrs)
    if(!colour){
    pout <- ggcells(sce, aes_string(x=paste0(DimRed,".1"), y=paste0(DimRed,".2"))) +
	geom_point(size=sz) +
	theme_void() +
	theme(legend.position="none",
	      text=element_text(size=18))
    } else {
    pout <- ggcells(sce, aes_string(x=paste0(DimRed,".1"), y=paste0(DimRed,".2"), color=Labels)) +
	geom_point(size=sz) +
	scale_color_manual(values=c25) +
	theme_void() +
	theme(legend.position="none",
	      text=element_text(size=18))
    }
    if(label) {
	pout <- pout + geom_label_repel(data=celltyp_ctrs, aes_string(x=paste0(DimRed,".1"), y=paste0(DimRed,".2"), label="Lbls", fill=NULL),fontface="bold",color="black") 
    } else {
	pout <- pout + geom_text_repel(data=celltyp_ctrs, aes_string(x=paste0(DimRed,".1"), y=paste0(DimRed,".2"), label="Lbls", fill=NULL),fontface="bold",color="black")
    }

    return(pout)
}
theme_pub <- function(base_size=16) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
  + theme(plot.title = element_text(face = "bold"),
          text = element_text(size=base_size),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1.2)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(size=rel(1)), 
          axis.line = element_line(colour="black",size=rel(1.5)),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
	  #           legend.direction = "horizontal",
	  #           legend.key.width= unit(1.0, "cm"),
	  #           legend.key.height = unit(0.5, "cm"),
          legend.title = element_blank(),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
  ))
}

# I use the umap funciton to also be able to extract the underlying knn graph
get_umap_graph <- function(umap_output) {
    require(Matrix)
    require(igraph)
    # extracts the graph information and outputs a weighted igraph
    dists.knn <- umap_output$knn[["distances"]]
    indx.knn <- umap_output$knn[["indexes"]]
    m.adj <- Matrix(0, nrow=nrow(indx.knn), ncol=nrow(indx.knn), sparse=TRUE) 
    rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn)

    for (i in seq_len(nrow(m.adj))) {
	m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- dists.knn[i,] 
    }

    igr <- graph_from_adjacency_matrix(m.adj, weighted=TRUE)
    return(igr)
}
