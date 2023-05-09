#!/usr/bin/env Rscript
#@author : Thomas NEFF

library(pheatmap)
library(DESeq2)
library(RColorBrewer)

# data --------------------------------------------------------------------

save_pheatmap <- function(x,file,n_samples,n_genes,c_width=80,c_height=15){
  stopifnot(!missing(x))
  stopifnot(!missing(file))
  png(filename = file,width = (c_width*n_samples)+(c_width/2),height=n_genes*c_height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

heatmap_dge <- function(y=stop("DGEList object missing"),
                        genes=NULL,pval=0.05,name.pval="FDR",
                        vars_samples_sort=list(names=NULL,decreasing=NULL),
                        vars_samples=NULL, label_samples=NULL,angle_col=45,
                        scale="row",graph=T,save_png=NULL,
                        color=colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(200),...){
  
  samples <- data.frame(eval(y)$samples)
  
  if (is.null(genes)) { genes <- rownames(eval(y)$tt[eval(y)$tt[eval(name.pval)] <= pval,]) }
  
  if (length(genes)>1 & is.character(genes)) { counts <- edgeR::cpm(y = eval(y),log=TRUE)[genes,]
  } else { stop("WARNING: the p-value threshold may be too low") }
  
  #order samples
  if (!is.null(vars_samples_sort$names)) {
    samples <- samples[order(samples[,vars_samples_sort$names],decreasing = vars_samples_sort$decreasing),]
    counts <- counts[,rownames(samples)]
    cluster_cols <- F
  } else { cluster_cols <- T}
  
  if (!is.null(label_samples)) {
    labels_col <- paste(unlist(samples[label_samples],use.names = F))
  } else { labels_col <- NULL }
  
  if (!is.null(vars_samples)) {
    annotation_col <- samples[paste(unlist(vars_samples,use.names = F))]
    colnames(annotation_col) <- names(vars_samples)
  } else { annotation_col <- NULL }
  
  mat_main <- as.matrix(counts)
  
  p <- pheatmap::pheatmap(mat = mat_main,color = color,scale = scale, angle_col = angle_col,
                          labels_col = labels_col,cluster_cols = cluster_cols,
                          annotation_col = annotation_col,...)
  
  if (!is.null(save_png)) {
    save_pheatmap(x = p,file = save_png,n_genes = dim(mat_main)[1],n_samples = dim(mat_main)[2])
  }
  
  if (graph) { print(p) }
  
  return(p)
}

# dds heatmap -------------------------------------------------------------


heatmap_dds <- function(dds=stop("DESeqDataSet object is missing"),
                        results=NULL,genes=NULL,pval=0.05,name.pval="padj",rlog.blind=F,
                        vars_samples_sort=list(names=NULL,decreasing=F),
                        vars_samples=NULL, label_samples=NULL,angle_col=45, scale="row",
                        color=colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(200),
                        save_png=NULL,graph=T,...){
  
  samples <- data.frame(eval(dds)@colData)
  
  if ("sizeFactor" %in% colnames(samples)) { samples <- samples[!colnames(samples) %in% "sizeFactor"]
  } else { stop("The dataset can be normalized")}

  if (is.null(genes)) {

    stopifnot(!is.null(results))

    results <- data.frame(eval(results))

    results <- results[complete.cases(results),]

    genes <- rownames(results[results[name.pval] <= pval,])
    }


  if (length(genes)>1 & is.character(genes)) { counts <- DESeq2::rlog(object = dds,blind = rlog.blind)[genes,]
  } else {stop("WARNING: the p-value threshold may be too low")}

  #order samples
  if (!is.null(vars_samples_sort$names)) {

    if (is.null(vars_samples_sort$decreasing)) {  vars_samples_sort$decreasing <- F}

    samples <- samples[order(samples[,vars_samples_sort$names],decreasing = vars_samples_sort$decreasing),]
    counts <- counts[,rownames(samples)]
    cluster_cols <- F
  } else {cluster_cols <- T}

   if (!is.null(label_samples)) {
     labels_col <- paste(unlist(samples[label_samples],use.names = F))
   } else { labels_col <- NULL }


  if (!is.null(vars_samples)) {
    annotation_col <- samples[paste(unlist(vars_samples,use.names = F))]
    if (!is.null(names(vars_samples))) { colnames(annotation_col) <- names(vars_samples) }
  } else { annotation_col <- NULL }

   mat_main <- as.matrix(assay(counts))

  p <- pheatmap::pheatmap(mat = mat_main,color = color,scale = scale, angle_col = angle_col,
                          labels_col = labels_col,cluster_cols = cluster_cols,
                          annotation_col = annotation_col,...)

  if (!is.null(save_png)) {
    save_pheatmap(x = p,file = save_png,n_genes = dim(mat_main)[1],n_samples = dim(mat_main)[2])
  }

  if (graph) { print(p) }

  return(p)
}




