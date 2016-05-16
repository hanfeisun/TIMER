#!/usr/bin/env Rscript

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



help_msg = 'Usage：
  For single run: Rscript regression2.R expFile cancer_catlogo
  For batch run: Rscript regression2.R --batch_input=table.txt

'
cat(help_msg)

cancers.available <- c('kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg', 'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct', 'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca', 'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol')

args <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
script.basename <- dirname(script.name)

baseDir = script.basename
#print(args)

# read in args, if exists --batch_input, the script will ignore other args.
args <- commandArgs(trailingOnly = TRUE)
batchf_index <- grepl("--batch_input=", args)
batch_file <- gsub("--batch_input=", "", args[batchf_index])

if (length(batch_file) != 0) {
  cancers <- as.matrix(read.table(batch_file, sep="\t"))
} else {
  cancers<- c(args[1], args[2])
  dim(cancers) <- c(1, 2)
}

for (i in seq(nrow(cancers))) {
  cancer.expFile <- cancers[i, 1]
  cancer.category <- cancers[i, 2]

  if (!(cancer.category %in% cancers.available)) {
    stop(paste('unknown cancers:', cancer.category))
  }

  gene.selected.marker.path <- paste(baseDir,
                          '/data/precalculated/genes_', cancer.category, '.RData',
                          sep='')
  gene.selected.marker <- get(load(gene.selected.marker.path))
  immune.agg.median <- get(load(paste(baseDir,
                                      '/data/precalculated/immune_median.RData',
                                      sep='')))

  cancer.expression <- ParseInputExpression(cancer.expFile)

  XX = immune.agg.median[gene.selected.marker, c(-4)]
  YY = cancer.expression[gene.selected.marker, ]

  fractions <- GetFractions.Abbas(XX, YY)
  print(fractions)

}



