# Created: 2018.07.31 - adapted from differential_expression_analysis_express.R and TNIP_master_readin.covariate.R
# By: Stephanie R Hilz
# Usage: Perform diff expression analysis in mass for a set of specified comparisons while
#   controlling for a covariate (in this case pool)

source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db", ask=FALSE, suppressUpdates=TRUE)
biocLite("edgeR", ask=FALSE, suppressUpdates=TRUE)
biocLite("goseq", ask=FALSE, suppressUpdates=TRUE)
library(org.Hs.eg.db)
library(edgeR)#call edgeR
library(goseq)

if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

makeCountsTable <- function(configFile, headerFile){
  config <- read.table(configFile, sep='\t',header=T, stringsAsFactors=F)
  config <- config[,c("file","name")]
  rsemheader <- read.table(headerFile,sep='\t', header=F, stringsAsFactors = F)
  rsemheader <- as.character(unlist(unname(rsemheader[1,])))
  rsemheader[which(rsemheader==' txStart')] <- 'txStart'
  for (i in seq_along(config$file)){
    fileName <- paste0('/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/20180813_ker_pool1_2_covariate_analysis/rsem_20180706_20180323/',config$file[i],'.rsem.txt')
    data <- read.table(fileName,sep='\t', header=F, stringsAsFactors = F)
    print(config$file[i])
    print(config$name[i])
    colnames(data) <- rsemheader
    data <- data[,c('Hugo','gene_expected_count')]
    data <- data[!duplicated(data),]
    print(dim(data))
    data$gene_expected_count <- round(data$gene_expected_count)
    colnames(data)[which(colnames(data)=='gene_expected_count')] <- config$name[i]
    if (i==1){
      merged <- data
    } else{
      merged <- merge(merged, data, by="Hugo")
    }
  }
  return(merged)
}

annotateGO <- function(GO.wall,invertedGetgo){
  GO.wall$foreground_genes <- 'NA'
  for (cat in GO.wall$category){
    GO.wall[which(GO.wall$category==cat),]$foreground_genes <- paste(invertedGetgo[[cat]],collapse=';')
  }
  return(GO.wall)
}

invertGetgo <- function(getgolist){
  inversion <- list()
  getgolist <- getgolist[!is.na(names(getgolist))]
  for (gene in names(getgolist)){
    for (goterm in getgolist[[gene]]){
      if (!goterm %in% names(inversion)){
        inversion[[goterm]] = c()
      }
      inversion[[goterm]] <- append(inversion[[goterm]], gene)
    }
  }
  return(inversion)
}

performGO <- function(binaryList, outfile){
  print("Table of input values")
  print(table(binaryList))
  pwf=nullp(binaryList,'hg19',"geneSymbol")
  GO.wall=goseq(pwf,"hg19","geneSymbol")
  print("Top 20 most significant GO terms")
  top <- GO.wall[,c(6,2)]
  colnames(top) <- c("term","pvalue")
  rownames(top) <- NULL
  print(head(top,20))
  getgolist <- getgo(names(binaryList[which(binaryList==1)]), 'hg19','geneSymbol')
  getGeneList <- invertGetgo(getgolist)
  GO.wall.anno <- annotateGO(GO.wall,getGeneList)
  write.csv(GO.wall.anno,outfile, quote=F)
}

## Define file paths
tag <- "20180921_fpool_1_2_"
configFile <- '/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/20180813_ker_pool1_2_covariate_analysis/configfile_20180813.txt' 
headerFile <- '/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/pooled_keratinocyte_orf/rsem_header.txt' 
CPMOutputFile <- paste0('/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/20180813_ker_pool1_2_covariate_analysis/20180921_with_GO/',tag,'CPMs.csv') 
boxFile <- paste0('/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/20180813_ker_pool1_2_covariate_analysis/20180921_with_GO/',tag,'boxplot.pdf') 
logFile <- paste0('/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/20180813_ker_pool1_2_covariate_analysis/20180921_with_GO/',tag,'log.csv')
DAOutputPath <- paste0('/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/20180813_ker_pool1_2_covariate_analysis/20180921_with_GO/',tag,"DE_anal.csv")
covariateToUse <- c('pool')
contrastToUse <- c('group')

## Inspect config and confirm groupings - can probably delete
config <- read.table(configFile, sep='\t',header=T, stringsAsFactors=F, comment.char='&')
config$descriptive <- paste0(config$orf.vector,'_',config$treatment,'_',config$time,'h')

## Merge all files by geneID to get counts file
merged <- makeCountsTable(configFile, headerFile)

## Format for edgeR
rownames(merged) <- merged$Hugo 
sampleTable_edgeR <- merged[,2:dim(merged)[2]]

## Set group variable after confirming order is the same as in edgeR
config$name == colnames(sampleTable_edgeR)
group <- config$group

## Build DGEList object
d<-DGEList(counts=sampleTable_edgeR,group=group)
dim(d)

## Filter out genes that aren't expressed in at least n samples (here n=3, because have 3 replicates for each grouping)
keep <- rowSums(cpm(d)>1) >= 6 
d<- d[keep,]
dim(d)

## Check CPM distribution across libraries before normalization
boxplot(cpm(d), ylim=c(0,100), col='red')
boxplot(log(cpm(d)+.001,2), col='red', ylab = "log2 CPM")

## Perform nomralizations
d<-calcNormFactors(d)#normalizing by dfault
d<-estimateCommonDisp(d)
d<-estimateTagwiseDisp(d)
sqrt(d$common.disp)

## Inspect normalization with boxplots
pdf(boxFile)
boxplot(log(cpm(d)+.001,2), col='red', ylab = "log2 CPM")
dev.off()

## Write CPMs to file
write.csv(cpm(d),CPMOutputFile)

## Define all possible pairwise comparisons
comparisons_forward <- combn(unique(group),2)
comparisons_reverse <- rbind(comparisons_forward[2,],comparisons_forward[1,])
comparisons <- cbind(comparisons_forward, comparisons_reverse)
print(dim(comparisons))

comparisonsToUse = c('egfp_none_24h_vs_egfp_IL17_24h',
                     'egfp_none_1h_vs_egfp_TNF_1h',
                     'egfp_none_24h_vs_tnip1_none_24h',
                     'egfp_none_24h_vs_tnfaip3_none_24h',
                     'egfp_IL17_24h_vs_tnip1_IL17_24h',
                     'egfp_IL17_24h_vs_tnfaip3_IL17_24h',
                     'egfp_TNF_1h_vs_tnip1_TNF_1h',
                     'egfp_TNF_1h_vs_tnfaip3_TNF_1h',
                     'egfp_none_24h_vs_egfp_TNF_24h',
                     'egfp_TNF_24h_vs_tnip1_TNF_24h',
                     'egfp_TNF_24h_vs_tnfaip3_TNF_24h',
                     'egfp_none_1h_vs_egfp_IL17_1h',
                     'egfp_IL17_1h_vs_tnip1_IL17_1h',
                     'egfp_IL17_1h_vs_tnfaip3_IL17_1h')

## Loop through each comparison and perform differential expression analysis
log <- rbind(comparisons, rep(NA,dim(comparisons)[2]))
log <- rbind(log, rep(NA,dim(comparisons)[2])) #will store counts of ups and downs per comparison
for (c in seq_len(dim(comparisons)[2])){
  ## Define the groups for this comparison
  group_a <- comparisons[1,c]
  group_b <- comparisons[2,c]
  description_a <- config[which(config$group == group_a),]$descriptive[1]
  description_b <- config[which(config$group == group_b),]$descriptive[1]
  descriptionComparison <- paste0(description_a,"_vs_",description_b)
  
  if (descriptionComparison %in% comparisonsToUse){#this is where the express part happens :)
    print(c)
    print(descriptionComparison)
    
    ## Create a subdirectory for all output
    print(paste0("Comparing group ",group_a," vs ",group_b,"(",description_a," vs ",description_b,")"))  
    subfolderPath <- paste0(DAOutputPath,tag,descriptionComparison)
    dir.create(subfolderPath, showWarnings <- TRUE)
    
    ## Build design matrix
    info <- config
    samples_a <- as.character(info[which(info$group == group_a),]$name)
    samples_b <- as.character(info[which(info$group == group_b),]$name)
    rownames(info) <- info$name
    info$name <- NULL
    info <- info[c(samples_a, samples_b),c(contrastToUse,covariateToUse)]
    covariate <- factor(info[,covariateToUse], levels=unique(info[,covariateToUse]))
    contrast_order <- c(group_a, group_b)
    contrast <- factor(info[,contrastToUse], levels=contrast_order)
    design <- model.matrix(~covariate+contrast)
    
    ## Do comparison-specific cpm filter and output CPMs (added 2018.04.19 because noticed some genes with low exp were being called sig DE)
    comparisonCPMOutputFile <- paste0(subfolderPath,'/',tag,descriptionComparison,"_CPMs.csv")
    keep <- rowSums(cpm(d[,c(config[which(config$group==group_a),]$name,
                             config[which(config$group==group_b),]$name)])>1) >= 6 #filters out anything without CPM>1 in at least six libraries #edit prior to every run
    y <- d[keep,c(samples_a, samples_b)]
    dim(y)
    write.csv(cpm(y),comparisonCPMOutputFile)
    
    ## Differential expression analysis controlling for covariate
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)
    
    ## Format results and grab top tags
    results <- topTags(qlf, n = nrow(qlf$table ) )$table
    head(results)
    detags <- rownames(results)[results$FDR < 0.05]
    deUp <- sum(results$FDR < 0.05 & results$logFC > 0)
    deDown <- sum(results$FDR < 0.05 & results$logFC < 0)
    
    ## Generate MA file
    MAFile <- paste0(subfolderPath,'/',tag,descriptionComparison,"_MA.pdf")
    pdf(MAFile)
    plotSmear(qlf, de.tags=detags) 
    dev.off()
    
    ## Generate output file
    DAOutputFile <- paste0(subfolderPath,'/',tag,descriptionComparison,"_diffexp_analysis.csv")
    write.csv(results,DAOutputFile)
   
    ## Save info into log
    log[3,c] <- deDown #downregulated genes
    log[4,c] <- deUp #upregulated genes
    
    ## Perform GO on significantly upregulated genes
    GOUpFile <- paste0(subfolderPath,'/',tag,descriptionComparison, '_Up_GeneOntology.csv')
    genes_up <- as.integer(results$logFC > 0 & 
                                 rownames(results) %in% detags)
    names(genes_up) <- rownames(results)
    performGO(genes_up, GOUpFile)
    
    # perform GO on significantly downregulated genes
    GODownFile <- paste0(subfolderPath,'/',tag,descriptionComparison, '_Down_GeneOntology.csv')
    genes_down=as.integer(results$logFC < 0 & 
                                  rownames(results) %in% detags)
    names(genes_down) <- rownames(results)
    performGO(genes_down, GODownFile)
  }
}

## Format and output log, which shows num up and down for each comp
log <- t(log)
colnames(log) <- c('group_a','group_b','down_b','up_b')
write.csv(log, file=logFile)

