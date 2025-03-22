#' Prepare Design matrix and Cross-subject Variance for MuSiC Deconvolution
#'
#' This function is used for generating cell type specific cross-subject mean and variance for each gene. Cell type specific library size is also calcualted.
#'
#' @param x SingleCellExperiment, single cell dataset
#' @param non.zero logical, default as TRUE. If true, remove all gene with zero expression.
#' @param markers vector or list of gene names. Default as NULL. If NULL, then use all genes provided.
#' @param clusters character, the colData used as clusters;
#' @param samples character,the colData used as samples;
#' @param select.ct vector of cell types. Default as NULL. If NULL, then use all cell types provided.
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data.
#' @param ct.cov logical. If TRUE, use the covariance across cell types.
#' @param verbose logical, default as TRUE.
#' @return a list of
#'  \itemize{
#'     \item {gene by cell type matrix of Design matrix;}
#'     \item {subject by celltype matrix of Library size;}
#'     \item {vector of average library size for each cell type;}
#'     \item {gene by celltype matrix of average relative abundance;}
#'     \item {gene by celltype matrix of cross-subject variation.}
#'     }
#' @importFrom Matrix rowSums
#'
#' @export
music_basis = function(x, non.zero = TRUE, markers = NULL, clusters, samples, select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE){
  x <- x[rowSums(Biobase::exprs(x))>0, ]
  clusters <- x$celltype
  samples <- x$sampleID
  M.theta <- sapply(unique(clusters), function(ct){
    MuSiC::my.rowMeans(sapply(unique(samples), function(sid){
      y = Biobase::exprs(x)[,clusters %in% ct & samples %in% sid]
      if(is.null(dim(y))){
        return(y/sum(y))
      }else{
        return(rowSums(y)/sum(y))
      }
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Relative Abudance Matrix...")}
  if(ct.cov){
    nGenes = nrow(x);
    n.ct = length(unique(clusters));
    nSubs = length(unique(samples))
    
    Theta <- sapply(unique(clusters), function(ct){
      sapply(unique(samples), function(sid){
        y = Biobase::exprs(x)[,clusters %in% ct & samples %in% sid]
        if(is.null(dim(y))){
          return(y/sum(y))
        }else{
          return( rowSums(y)/sum(y) )
        }
      })
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Theta))
      Theta = Theta[, m.ct]
    }
    
    Sigma.ct = sapply(1:nGenes, function(g){
      sigma.temp = Theta[nGenes*(0:(nSubs - 1)) + g, ];
      Cov.temp = cov(sigma.temp)
      Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes*(0:(nSubs - 1)) + 1, ])) == 0, ])
      Cov.temp[which(colSums(is.na(sigma.temp))>0), ] = Cov.temp1[which(colSums(is.na(sigma.temp))>0), ]
      Cov.temp[, which(colSums(is.na(sigma.temp))>0)] = Cov.temp1[, which(colSums(is.na(sigma.temp))>0)]
      return(Cov.temp)
    })
    colnames(Sigma.ct) = rownames(x);
    
    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma.ct <- Sigma.ct[ , m.ids]
    }
    if(verbose){message("Creating Covariance Matrix...")}
  }else{
    Sigma <- sapply(unique(clusters), function(ct){
      apply(sapply(unique(samples), function(sid){
        y = Biobase::exprs(x)[,clusters %in% ct & samples %in% sid]
        if(is.null(dim(y))){
          return(y/sum(y))
        }else{
          return(rowSums(y)/sum(y))
        }
      }), 1, var, na.rm = TRUE)
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Sigma))
      Sigma = Sigma[, m.ct]
    }
    
    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma <- Sigma[m.ids, ]
    }
    if(verbose){message("Creating Variance Matrix...")}
  }
  
  S <- sapply(unique(clusters), function(ct){
    MuSiC::my.rowMeans(sapply(unique(samples), function(sid){
      y = Biobase::exprs(x)[, clusters %in% ct & samples %in% sid]
      if(is.null(dim(y))){
        return(sum(y))
      }else{
        return(sum(y)/ncol(y))
      }
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Library Size Matrix...")}
  
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  #S.ra = relative.ab(S, by.col = FALSE)
  #S.ra[S.ra == 0] = NA
  #S[S == 0] = NA
  #M.S = mean(S, na.rm = TRUE)*ncol(S)*colMeans(S.ra, na.rm = TRUE)
  
  if(!is.null(cell_size)){
    if(!is.data.frame(cell_size)){
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }else if(sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))){
      stop("Cell type names in cell_size must match clusters")
    }else if (any(is.na(as.numeric(cell_size[, 2])))){
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  
  D <- t(t(M.theta)*M.S)
  
  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }
  
  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }
  
  if(ct.cov){
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma.ct = Sigma.ct))
  }else{
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma = Sigma))
  }
}
