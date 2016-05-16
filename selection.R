## License GPL2.0
## Developed by Bo Li, bli@jimmy.harvard.edu, 2016

## Select genes whose expression levels are negatively correlated with tumor purity,
## Select them because they are informative to deconvolve immune calls in the tumor tissue


options(warn=-1)
library(sqldf)
library(sva)
library(crayon)

args <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
script.basename <- dirname(script.name)
baseDir = script.basename


cancers <- c('kich','blca','brca','cesc','gbm','hnsc','kirp','lgg','lihc','luad','lusc','prad','sarc','pcpg','paad','tgct','ucec','ov','skcm','dlbc','kirc','acc','meso','thca','uvm','ucs','thym','esca','stad','read','coad','chol')


## A hack for assigning from multiple return type
## Such as:  list[a, b] <- functionReturningTwoValues()
## Pasted from https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
     a <- args[[i]]
     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}


LoadCancerGeneExpression <- function(cancer) {
  ## Load and process gene expression data

  ## Args:
  ##   cancer: String of cancer type
  ##
  ## Returns:
  ##   The data frame of gene expression for the input cancer category 
  ##   (col for sample, row for gene)



  if (cancer %in% c('gbm', 'ov')) {
    geneExpression <- get(load(paste(baseDir,'/data/RNAseq/',cancer,'Affy.Rdata',
                                     sep='')))
  } else {
    geneExpression <- get(load(paste(baseDir,'/data/RNAseq/',cancer,'RNAseq.Rdata',
                                     sep='')))
  }

  geneExpression <- as.matrix(geneExpression)
  mode(geneExpression) <- 'numeric'
  if(!cancer %in% c('gbm', 'ov', 'esca', 'stad')) {
    ## rsem scaled estimates needs multiply 1e6, Array or RPKM does not need.
    geneExpression <- geneExpression * 1e6   
  }

  geneExpression <- ConvertRownameToLoci(geneExpression)
  colnames(geneExpression) <- GetCancerSampleID(colnames(geneExpression), 4)

  return(geneExpression)
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
  return(cancerGeneExpression[tmp.vv,])
}

immuneCuratedData <- paste(baseDir,
                           '/data/precalculated/immune.expression.curated.RData',
                           sep='')

LoadImmuneGeneExpression <- function() {
  ## Load gene expression data for immune cells

  ## Returns:
  ##   A data frame of expression data for immune cells
  ##   (cols for immune cell sample, rows for "gene name;probe ID")
  exp <- get(load(paste(baseDir,'/data/immune_datasets/HPCTimmune.Rdata',sep='')))

  ##----- Select single reference samples of pre-selected immune cell types -----##
  B_cell <- 362:385
  T_cell.CD4 <- grep('T_cell.CD4',colnames(exp))
  T_cell.CD8 <- grep('T_cell.CD8',colnames(exp))
  NK <- 328:331
  Neutrophil <- 344:361
  Macrophage <- 66:80
  DC <- 151:238

  curated.ref <- exp[, c(B_cell, T_cell.CD4, T_cell.CD8,
                         NK, Neutrophil, Macrophage, DC)]

  curated.cell.types <- colnames(curated.ref)
  names(curated.cell.types) <- c(rep('B_cell', length(B_cell)),
                                 rep('T_cell.CD4', length(T_cell.CD4)),
                                 rep('T_cell.CD8', length(T_cell.CD8)),
                                 rep('NK', length(NK)),
                                 rep('Neutrophil', length(Neutrophil)),
                                 rep('Macrophage', length(Macrophage)),
                                 rep('DC', length(DC)))

  curated.ref.genes <- ConvertImmuneProbeToRefgene(curated.ref)
  return(list(genes=curated.ref.genes, celltypes=curated.cell.types))
}


ConvertImmuneProbeToRefgene <- function(curated.ref){
  ##----- function to preprocess the reference dataset, not necessary if the processed data "curated.ref.genes.Rdata" is available -----##

  if (file.exists(immuneCuratedData)) {
    curated.ref.genes <- get(load(immuneCuratedData))
    return(curated.ref.genes)
  }

  tmpDD <- data.frame(curated.ref)
  tmpDD <- tmpDD[order(rownames(tmpDD)), ]
  ## sort the immune expression data by rownames

  colnames(tmpDD) <- gsub('\\.', '_', colnames(tmpDD))
  genes <- sapply(strsplit(rownames(tmpDD), ';'), function(x) x[[1]])
  ## remove the probe ID, only keep gene names

  tmpDD <- cbind(genes, tmpDD)
  tmpDD <- tmpDD[order(genes), ]
  ## sort by gene names

  tmp0 <- c()
  cnt <- 0
  for(i in colnames(tmpDD)[2:ncol(tmpDD)]){
    ## start from the second column (the first column is gene information)

    print(cnt)
    print(i)
    tmp <- sqldf(paste('select max(', i, ') from tmpDD group by genes', sep <- ''))
    ## select the maximum probe expression level when a gene has multiple probes

    if(length(tmp0) == 0) tmp0 <- tmp else tmp0 <- cbind(tmp0, tmp)
    cnt <- cnt + 1
  }
  colnames(tmp0) <- colnames(tmpDD)[2:ncol(tmpDD)]
  rownames(tmp0) <- unique(tmpDD[, 1])
  curated.ref.genes <- tmp0
  save(curated.ref.genes, file=immuneCuratedata)
  return(curated.ref.genes)
}

LoadImmuneMarkerProbe <- function () {
  ## Load immune marker probe from Abbas et al., 2005

  ## Return:
  ##   Vector of probeID for immule marker genes
  ##   (name: cancer category, value: probe ID)

  tmp <- read.csv(paste(baseDir, 'data/immune_datasets/IRIS-marker-gene.txt', sep='/'),
                  header=T, sep='\t', stringsAsFactors=F)
  marker.list <- tmp[, 1]
  names(marker.list) <- tmp[, 7]
  names(marker.list) <- gsub(' ','_',tmp[,7])
  names(marker.list) <- gsub('Dendritic_Cell', 'DC', names(marker.list))
  names(marker.list) <- gsub('Neutrophil', 'Neutrophils', names(marker.list))
  gsub('Cell', 'cell', names(marker.list)) -> names(marker.list)
  return(marker.list)
}

LoadImmuneMarkerGene <- function() {
  ## Load immune marker genes from Abbas et al., 2005
  ## Convert the result of LoadImmuneMarkerProbe() to gene names with the help of immune gene expression data

  ## Args:
  ##    ref: The immune gene expression data which is generated by LoadImmuneGeneExpression()
  ## Returns:
  ##    Vector of immune marker gene names


  ref <- get(load(paste(baseDir,'/data/immune_datasets/HPCTimmune.Rdata',sep='')))
  tmp <- strsplit(rownames(ref), ';')
  AffyIDtoGenes <- sapply(tmp, function(x) x[[1]])
  names(AffyIDtoGenes) <- sapply(tmp, function(x) x[[2]])

  marker.list <- LoadImmuneMarkerProbe()
  marker.list.genes <- AffyIDtoGenes[marker.list]
  return(marker.list.genes)
}

LoadCancerTumorPurity <- function(cancer) {
  ## Load and process tumor purity data

  ## Args:
  ##   cancer: String of cancer type
  ##
  ## Return:
  ##    A dataframe of tumor purity data

  AGP <- read.table(paste(baseDir, '/data/AGP/AGP-', cancer, '.txt', sep=''),
                    sep='\t', header=T)
  AGP <- AGP[which(AGP[, 'PoP'] > 0.01), ]
  rownames(AGP) <- GetCancerSampleID(AGP[, 1], 4)
  return(AGP)
}

GetCancerSampleID <- function(sID, num.res=3){
  ## Edit TCGA ID, with the option of keeping the first num.res fields

  ## Args:
  ##    sID: vector of original TCGA sample IDs
  ##    num.res: how many fields to keep

  ## Return:
  ##    Vector of sample IDs with the extracted fields


  mm <- c()
  for (id in sID) {
    tmp <- unlist(strsplit(id, '-'))
    if (length(tmp) == 1) {
      tmp <- unlist(strsplit(id, '\\.'))
    }
    ll <- 'TCGA'
    for (j in 2:num.res) {
      ll <- paste(ll, tmp[j], sep='-')
    }
    mm <- c(mm, ll)
  }
  return(mm)
}


RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {

  ## intersect the gene names of cancer.exp and immune.exp
  tmp.dd <- as.matrix(cancer.exp)
  tmp <- sapply(strsplit(rownames(cancer.exp), '\\|'), function(x) x[[1]])
  rownames(tmp.dd) <- tmp
  tmp.dd <- tmp.dd[which(nchar(tmp)>1), ]
  tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))

  ## bind cancer and immune expression data into one dataframe
  N1 <- ncol(tmp.dd)

  tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
  tmp.dd <- as.matrix(tmp.dd)
  mode(tmp.dd) <- 'numeric'

  ## remove batch effects
  N2 <- ncol(immune.exp)
  tmp.batch <- c(rep(1, N1), rep(2, N2))
  tmp.dd0 <- ComBat(tmp.dd, tmp.batch, c())

  ## separate cancer and immune expression data after batch effect removing
  dd.br <- tmp.dd0[, 1:N1]
  immune.exp.br <- tmp.dd0[, (N1+1):(N1+N2)]

  ## a immune category has multiple samples, use the median expression level for a gene
  tmp0 <- c()
  for(kk in unique(names(immune.cellType))){
    tmp.vv <- which(names(immune.cellType)==kk)
    tmp0 <- cbind(tmp0, apply(immune.exp.br[, tmp.vv], 1, median, na.rm=T))
  }


  immune.exp.agg.br <- tmp0
  colnames(immune.exp.agg.br) <- unique(names(immune.cellType))
  return(list(dd.br, immune.exp.br, immune.exp.agg.br))
}


##----- function to select genes with expression values negatively correlated with tumor purity -----##
GetPurityGenes <- function(dd, AGP, thr.p=0.05, thr.c=0, mode='env'){
  tmp.ss <- intersect(colnames(dd), rownames(AGP))
  if(length(tmp.ss)==0){
    ## TODO: check whether this paragraph is necessary
    colnames(dd) <- getID(colnames(dd))
    tmp.ss <- intersect(colnames(dd), rownames(AGP))
  }

  tmp.dd <- dd[, tmp.ss]


  ## 's' means Spearman Correlation
  tmp <- lapply(rownames(tmp.dd),
                function(x) {
                  cor.test(tmp.dd[x, ],
                           as.numeric(AGP[colnames(tmp.dd), 2]),
                           method='s')})

  tmp.pp <- sapply(tmp, function(x)x$p.value)
  tmp.cor <- sapply(tmp, function(x)x$estimate)
  names(tmp.pp) <- names(tmp.cor) <- rownames(dd)

  if (mode == 'env') {
    geneNames <- names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
  } else if (mode == 'tumor') {
    geneNames <- names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
  } else {
    stop("invalid mode for purity-based gene selection")
  }
  return(geneNames)
}

##----- remove outlier genes whose expression may drive the colinearity of similar covariates in the regression -----##
RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
  ## removes upper thr.q quantile for every reference feature
  remove.vv=c()
  for(i in 1:ncol(ref.dd)){
    tmp=quantile(ref.dd[vv,i],thr.q)[1]
    tmp.vv=which(ref.dd[vv,i]>tmp)
    remove.vv=c(remove.vv,tmp.vv)
  }
  remove.vv=unique(remove.vv)
  return(vv[-remove.vv])
}




##----- setup parameters and establish the output file -----##



for(cc in cancers) {
  if (cc == 'skcm') {
    cc.type <- '06A'
  } else {
    cc.type <- '01A'
  }

  geneCalculated <- paste(baseDir,
                        '/data/precalculated/genes_', cc, '.RData',
                        sep='')

  
  ccDir = paste(baseDir, 'results', cc, sep='/')
  dir.create(ccDir, showWarnings=FALSE, recursive=TRUE)



  options(scipen=6)


  write(paste(cc, ' output\n'),
        file=paste(baseDir, '/results/', cc, '/output-statistics.txt', sep=''))


  cancer.geneExpression <- LoadCancerGeneExpression(cancer=cc)
  cancer.tumorPurity <- LoadCancerTumorPurity(cancer=cc)

  immune <- LoadImmuneGeneExpression()
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes


  immune.markerGene <- LoadImmuneMarkerGene()

  cat(green(paste("## Removing the batch effect of", cc, '\n')))
  list[cancer.expNorm, immune.expNorm, immune.expNormMedian] <-
    RemoveBatchEffect(cancer.geneExpression, immune.geneExpression, immune.cellTypes)


  gene.purity.neg <- GetPurityGenes(cancer.geneExpression, cancer.tumorPurity,
                                    thr.p = 0.05, thr.c = -0.2)
  gene.purity.neg <- intersect(gene.purity.neg, rownames(immune.expNormMedian))

  gene.selected.marker <- intersect(gene.purity.neg, immune.markerGene)


  ##---- calculate differences between the correlations of reference immune cells using Pearson's or Spearman's correlations -----##
  tmp.diff <-
    sum(sum(
      abs(cor(immune.expNormMedian[gene.selected.marker, ], method='p') -
          cor(immune.expNormMedian[gene.selected.marker, ], method='s'))))

  if (tmp.diff >= -10000) {
    gene.selected.marker <- RemoveOutliers(gene.selected.marker,
                                           immune.expNormMedian[, -4])
  }

  cat("Number of genes inversely correlated with purity is ",
      length(gene.purity.neg), '\n\n',
      sep='', file='output-statistics.txt', append=T)

  cat("Number of immune genes inversely correlated with purity is ",
      length(gene.selected.marker), '\n\n',
      sep='', file='output-statistics.txt', append=T)

  ##----- calculate the significance of enrichment for purity selected genes to immune marker genes -----##
  gene.all=intersect(rownames(immune.expNormMedian), rownames(cancer.expNorm))
  n.immune=length(intersect(immune.markerGene, gene.all))

  cat("Test if immune genes are enriched for inverse correlation with purity: \n\n", file='output-statistics.txt', append=T)

  save(gene.selected.marker, file=geneCalculated, compress='gzip')

  sink(file='output-statistics.txt', append=T)
  print(
    fisher.test(
      matrix(
        c(
          length(gene.selected.marker),
          length(gene.purity.neg)-length(gene.selected.marker),
          n.immune,
          length(gene.all)-n.immune),
        2, 2)))
  sink()

}
