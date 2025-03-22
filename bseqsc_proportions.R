bseqsc_proportions <- function(x, reference = NULL, log = NULL, ..., verbose = TRUE){
  
  message <- if( verbose ) message else function(...) invisible(NULL)
  
  # load pancreas basis matrix if necessary
  if( is.null(reference) ){
    message("* Using pancreatic islet reference basis matrix: ", appendLF = FALSE)
    reference <- ldata('PancreasIslet')
  }
  
  # load CIBERSORT
  y <- x
  if( is(y, 'ExpressionSet') ) y <- exprs(y)
  
  # use common features only
  ids <- intersect(rownames(y), rownames(reference))
  message(sprintf("* Data features: %s", pkgmaker::str_out(rownames(y), total = TRUE)))
  message(sprintf("* Basis features: %s", pkgmaker::str_out(rownames(reference), total = TRUE)))
  message(sprintf("* Common features: %s", pkgmaker::str_out(ids, total = TRUE)))
  islog <- xbioc::is_logscale(y)  
  if( islog ){
    message("* Converting to linear scale")
    y <- expb(y, 2)
  }
  
  # setup run directory and files
  tdir <- tempfile('CIBERSORT')
  dir.create(tdir)
  owd <- setwd( tdir )
  on.exit({
    unlink(tdir, recursive = TRUE)
    setwd(owd)
  }) 
  message("* Writing input files ... ", appendLF = FALSE)
  write.table(reference, file = xf <- 'reference.tsv', sep = "\t", row.names = TRUE, col.names = NA)
  write.table(y, file = yf <- 'mixture.tsv', sep = "\t", row.names = TRUE, col.names = NA)
  message('OK')
  
  # run
  message("* Running CIBERSORT ... ", appendLF = FALSE)
  res <- CIBERSORT(xf, yf, QN=FALSE,...)
  message('OK')
  
  stats <- setdiff(colnames(res), colnames(reference))
  list(coefficients = t(res[, !colnames(res) %in% stats]), stats = res[, stats])
}
