#!/usr/bin/env Rscript


baseDir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
  script.basename <- dirname(script.name)
  return(script.basename)
})()

ParseArgs <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  # read in args, if exists --batch_input, the script will ignore other args.
  tmp <- grepl("--batch-input=", args)
  batch.file <- gsub("--batch-input=", "", args[tmp])
  print(batch.file)

  if (length(batch.file) == 0) {
    cancer.expression <- args[1]
    cancer.category <- args[2]
  } else {
    cancer.expression <- NULL
    cancer.category <- NULL
  }
  return(list(batch = batch.file, expression = cancer.expression, category = cancer.category))
}


## Evaluate the code immediately so that error can be detected as early as possible
cancers <- (function() {
  args <- ParseArgs()
  cancers.available <- c('kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg',
                         'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct',
                         'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca',
                         'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol')

  if (length(args$batch) != 0) {
    cat("Enter batch mode\n")
    cancers <- as.matrix(read.table(args$batch, sep=","))
  } else {
    cancers<- c(args$expression, args$category)
    dim(cancers) <- c(1, 2)
  }
  print(cancers)
  for (i in seq(nrow(cancers))) {
    cancer.category <- cancers[i, 2]
    if (!(cancer.category %in% cancers.available)) {
      stop(paste('unknown cancers:', cancer.category))
    }
  }
  return(cancers)
})()


source(paste(baseDir, '/utils.r', sep=''))

##----- Constrained regression method implemented in Abbas et al., 2009 -----##
GetFractions.Abbas <- function(XX, YY, w=NA){
  ## XX is immune expression data
  ## YY is cancer expression data
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0, ncol(XX)))
      tmp.XX=tmp.XX[, -ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0, ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX, YY, intercept=F) else tmp=lsfit(tmp.XX, YY, w, intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0, ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}

ConvertRownameToLoci <- function(cancerGeneExpression) {
  ## Extract only the loci information for row name

  ## Example of origin row name is 'LOC389332|389332'
  ## Coverted row name is 'LOC389332'

  ## Args:
  ##   geneExpression: the orginal geneExpression load from .Rdata file
  ##
  ## Returns:
  ##   Modified geneExpression

  tmp <- strsplit(rownames(cancerGeneExpression), '\\|')
  tmp <- sapply(tmp, function(x) x[[1]])
  tmp.vv <- which(nchar(tmp) > 1)
  rownames(cancerGeneExpression) <- tmp
  extracted <- as.matrix(cancerGeneExpression[tmp.vv, ])
  colnames(extracted) <- colnames(cancerGeneExpression)
  return(extracted)
}

ParseInputExpression <- function(path) {
  ret <- read.csv(path, sep='\t', row.names=1)
  ret <- as.matrix(ret)
  mode(ret) <- 'numeric'
  ret <- ConvertRownameToLoci(ret)
  return(ret)
}


DrawQQPlot <- function(cancer.exp, immune.exp, name='') {
  ## q-q plot by sample should look like a straight line.
  ## Extreme values may saturate for Affy array data, but most of the data should align well.
  qqplot(cancer.exp[,1], immune.exp[,1], xlab='Tumor Expression', ylab='Ref Expression',
         main='Sample-Sample Q-Q plot')
  mtext(name, col="gray11")
}


main <- function() {

#   help_msg = 'Usage：
#   For single run: Rscript regression.R expFile cancer_catlog
#   For batch run: Rscript regression.R --batch_input=table.txt
# '
  # cat(help_msg)

  TimerINFO('Loading immune gene expression')
  immune <- LoadImmuneGeneExpression()
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes
  if (!dir.exists(paste(baseDir, '/results', sep=''))) {
    dir.create(paste(baseDir, '/results', sep=''))
  }
  pdf(paste(baseDir, '/results/output.pdf', sep=''))

  for (i in seq(nrow(cancers))) {
    cancer.expFile <- cancers[i, 1]
    cancer.category <- cancers[i, 2]
    gene.selected.marker.path <- paste(baseDir,
                                       '/data/precalculated/genes_', cancer.category, '.RData',
                                       sep='')
    gene.selected.marker <- get(load(gene.selected.marker.path))
    cancer.expression <- ParseInputExpression(cancer.expFile)

    TimerINFO(paste("Removing the batch effect of", cancer.expFile))
    DrawQQPlot(cancer.expression, immune.geneExpression, name=cancer.expFile)

    tmp <- RemoveBatchEffect(cancer.expression, immune.geneExpression, immune.cellTypes)
    cancer.expNorm <- tmp[[1]]
    immune.expNormMedian <- tmp[[3]]

    DrawQQPlot(cancer.expNorm, immune.expNormMedian,
               name=paste("After batch removing and aggregating for", cancer.expFile))


    XX = immune.expNormMedian[gene.selected.marker, c(-4)]
    YY = cancer.expNorm[gene.selected.marker, ]


    fractions <- GetFractions.Abbas(XX, YY)
    print(fractions)
  }

  dev.off()


}

main()

