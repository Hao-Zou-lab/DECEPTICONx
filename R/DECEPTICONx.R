#########bayesprism############

#' @title Creat BayesPrism base
#'
#' @param bk.dat is a matrix with rows for gene names and columns for sample names
#' @param sc.dat is a matrix with rows for gene names and columns for cell names
#' @param cell.type.labels is an n*2 matrix, the first column is the cell name of sc.dat, and the second column is the cell subtype
#'
#' @return BayesPrism base
#' @export
#'
#' @examples bayes_base(bulk.samples,sc.dat,subtype)
bayes_base <- function(bk.dat,sc.dat,cell.type.labels,species="hs",input.type="count.matrix",gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"),gene.type = "protein_coding"){
  library(BayesPrism)
  sc.dat = t(sc.dat)
  bk.dat = t(bk.dat)
  cell.type.labels <- t(cell.type.labels)
  if (species == "hs" || species == "mm") {
    sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                      input.type=input.type,
                                      species=species,
                                      gene.group=gene.group ,
                                      exp.cells=1)
    sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                            gene.type = gene.type)
  } else {
    # do something else
    sc.dat.filtered.pc <- sc.dat
  }
  myPrism <- new.prism(
    reference=sc.dat.filtered.pc,
    mixture=bk.dat,
    input.type=input.type,
    cell.type.labels = cell.type.labels,
    cell.state.labels = NULL,
    key=NULL,
    outlier.cut=0.01,
    outlier.fraction=0.1,
  )
  write.table(t(myPrism@phi_cellState@phi),file = "./custom_signature_matrix/Bayes_base.txt",row.names = TRUE,col.names = TRUE,sep='\t')
  message("Creat BayesPrism base sucessful!")
  }

#########monocle3############

#' @title Creat monocle3 base
#'
#' @param expression_matrix is a matrix with rows for gene names and columns for cell names
#' @param cell_metadata is an n*2 matrix, the first column is the cell name of sc.dat, and the second column is the cell subtype
#'
#' @return monocle3 base
#' @export
#'
#' @examples monocle3_base <- monocle3_base(sc.dat,subtype)
monocle3_base <- function(expression_matrix,cell_metadata){
  library(monocle3)
  row.names(cell_metadata) <- colnames(expression_matrix)
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
  )
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds, resolution=1e-5)

  rowname <- cell_group <- marker_score <- cell_id <- mean_expression <- NULL
  fraction_expressing <- specificity <- pseudo_R2 <- NULL
  lrtest_p_value <- lrtest_q_value <- gene_short_name <- NULL
  cell_group_df <- data.frame(row.names = row.names(colData(cds)),
                              cell_id = row.names(colData(cds)))

  cell_group_df$cell_group <- as.character(cell_metadata[,1])
  cluster_mean_exprs = as.matrix(aggregate_gene_expression(cds,
                                                           cell_group_df = cell_group_df, norm_method = "size_only",
                                                           scale_agg_values = FALSE))
  write.table(cluster_mean_exprs, './custom_signature_matrix/monocle3_mean_base.txt', sep='\t', row.names= T, col.names= T, quote=F)
  message("Creat monocle3 base sucessful!")
  }
#####################################################################
#' @title Creat c
#'
#' @param sc.sce is a matrix with rows for gene names and columns for cell names
#' @param bulk is a matrix with rows for gene names and columns for sample names
#' @param cell.type.labels is an n*2 matrix, the first column is the cell name of sc.dat, and the second column is the cell subtype
#'
#' @return MuSiC2 base
#'
MuSiC2_base <- function(sc.sce, bulk, cell.type.labels,
                        prop_r = 0.1, eps_c = 0.05, expr_low = 20, cutoff_expr = 0.05,
                        maxiter = 200, markers = NULL, cell_size = NULL, centered = FALSE,
                        ct_cov = FALSE, normalize = FALSE, cutoff_fc = 2, eps_r = 0.01, n_resample = 20, sample_prop = 0.5,
                        cutoff_c = 0.05, cutoff_r = 0.01, iter_max = 1000, nu = 1e-04,
                        eps = 0.01,  non.zero = FALSE) {
  library(MuSiC)
  library(MuSiC2)
  select.ct = unique(cell.type.labels)[,1]
  bulk.control.mtx = bulk
  bulk.case.mtx = bulk
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if (length(gene.bulk) < 0.1 * min(nrow(bulk.control.mtx),
                                    nrow(bulk.case.mtx))) {
    stop("Not enough genes for bulk data! Please check gene annotations.")
  }
  bulk.mtx = bulk.control.mtx[gene.bulk, ]
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if (length(gene_all) < 0.2 * min(length(gene.bulk), nrow(sc.sce))) {
    stop("Not enough genes between bulk and single-cell data! Please check gene annotations.")
  }
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr >= expr_low])
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)]
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)]
  bulk.mtx = bulk.control
  bulk.gene = rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
  bulk.mtx = bulk.mtx[bulk.gene, ]
  if (is.null(markers)) {
    sc.markers = bulk.gene
  } else {
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  markers = sc.markers
  if (non.zero) {
    all_zero_rows <- apply(sc.sce, 1, function(row) all(row == 0))
    sc.sce<- sc.sce[!all_zero_rows, ]
  }
  clusters <- cell.type.labels[,1]
  samples <- cell.type.labels[,1]
  M.theta <- sapply(unique(clusters), function(ct) {
    my.rowMeans(sapply(unique(samples), function(sid) {
      y = sc.sce[, clusters %in% ct & samples %in%
                   sid]
      if (is.null(dim(y))) {
        return(y/sum(y))
      }else {
        return(rowSums(y)/sum(y))
      }
    }), na.rm = TRUE)
  })
  S <- sapply(unique(clusters), function(ct) {
    my.rowMeans(sapply(unique(samples), function(sid) {
      y = sc.sce[, clusters %in% ct & samples %in%
                   sid]
      if (is.null(dim(y))) {
        return(sum(y))
      } else {
        return(sum(y)/ncol(y))
      }
    }), na.rm = TRUE)
  })
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  if (!is.null(cell_size)) {
    if (!is.data.frame(cell_size)) {
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
      stop("Cell type names in cell_size must match clusters")
    } else if (any(is.na(as.numeric(cell_size[, 2])))) {
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[,
                                                      1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),
    ]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  D <- t(t(M.theta) * M.S)
  if (!is.null(select.ct)) {
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }
  if (!is.null(markers)) {
    ids <- intersect(unlist(markers), rownames(sc.sce))
    m.ids = match(ids, rownames(sc.sce))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }
  write.table(D,file = "./custom_signature_matrix/MuSiC2_base.txt",row.names = TRUE,col.names = TRUE,sep='\t')
  message("Creat MuSiC2 base sucessful!")
}
##################################################################

#' @title generate_signature_matrix
#
#' @param cell_types can enter in nk, b, cd4, cd8, fibroblast, five kinds of macrophage cell types
#'
#' @return
#' @export
#'
#' @examples
generate_signature_matrix <- function(cell_types) {
  # 创建一个字典来存储每种细胞类型对应的signature.matrix文件名
  cell_type_files <- list(
    "nk" = c("NK_monocle3.txt", "NK_BayesPrism.txt", "NK_MuSiC2.txt"),
    "b" = c("B_monocle3.txt", "B_BayesPrism.txt", "B_MuSiC2.txt"),
    "cd4" = c("CD4_monocle3.txt", "CD4_BayesPrism.txt", "CD4_MuSiC2.txt"),
    "cd8" = c("CD8_monocle3.txt", "CD8_BayesPrism.txt", "CD8_MuSiC2.txt"),
    "fibroblast" = c("Fibroblast_monocle3.txt", "Fibroblast_BayesPrism.txt", "Fibroblast_MuSiC2.txt"),
    "macrophage" = c("Macrophage_monocle3.txt", "Macrophage_BayesPrism.txt", "Macrophage_MuSiC2.txt")
  )

  # 将输入的cell_types转换为小写，并去除空格
  cell_types <- tolower(trimws(cell_types))

  # 检查输入的细胞类型是否有效
  valid_cell_types <- names(cell_type_files)
  invalid_types <- cell_types[!cell_types %in% valid_cell_types]
  if (length(invalid_types) > 0) {
    warning(paste("无效的细胞类型:", paste(invalid_types, collapse = ", ")))
  }

  # 为有效的细胞类型生成signature.matrix
  signature_files <- list()
  for (cell_type in cell_types) {
    if (cell_type %in% valid_cell_types) {
      signature_files[[cell_type]] <- file.path("./signature_matrix", cell_type_files[[cell_type]])
    }
  }

  # 合并所有有效的signature.matrix文件路径
  signature.matrix <- unlist(signature_files)

  return(signature.matrix)
}



######################################################

#' @title run_DECEPTICONx
#'
#' @param bulk.samples is a matrix with rows for gene names and columns for sample names
#' @param sc.dat is a matrix with rows for gene names and columns for cell names
#' @param subtype is an n*2 matrix, the first column is the cell name of sc.dat, and the second column is the cell subtype
#' @param light set TRUE for light version of DECEPTICON, running time reduced.
#' @param custom.template
#' @param custom.signature set TRUE for custom signature matrix.
#' @param signature.matrix When custom.signature is set to TRUE, need to enter expression template(s).
#' @param cell_types can enter in nk, b, cd4, cd8, fibroblast, five kinds of macrophage cell types
#'
#' @return
#' @export
#'
#' @examples
run_DECEPTICONx <- function (bulk.samples,sc.dat,subtype,cell_types=c("nk","b","cd4","cd8","fibroblast","macrophage"), RUNpath, light = TRUE, custom.template = FALSE ,custom.signature = FALSE,
                             single.cell = FALSE, single.cell.data = NULL)
{
  library(DECEPTICON)
  library(pkgmaker)
  library(xbioc)
  light = light
  single.cell.data = single.cell.data
  setwd(RUNpath)
  if (custom.template == FALSE) {
    dir.create("res")
    signature.matrix = generate_signature_matrix(cell_types)
    run_DECEPTICON(bulk.samples = bulk.samples,
                   RUNpath = RUNpath,
                   custom.signature = TRUE,
                   signature.matrix = signature.matrix,
                   light = TRUE)
  }
  else {
    bayes_base(bulk.samples,sc.dat,subtype)
    monocle3_base <- monocle3_base(sc.dat,subtype)
    MuSiC2_base(sc.dat,bulk.samples, subtype)
    run_DECEPTICON(bulk.samples = bulk.samples,
                   RUNpath = RUNpath,
                   custom.signature = TRUE,
                   signature.matrix = c('./custom_signature_matrix/monocle3_mean_base.txt',
                                        './custom_signature_matrix/Bayes_base.txt',
                                        './custom_signature_matrix/MuSiC2_base.txt'),
                   light = TRUE)
  }

  message("Run DECEPTICONx sucessful!")
}
