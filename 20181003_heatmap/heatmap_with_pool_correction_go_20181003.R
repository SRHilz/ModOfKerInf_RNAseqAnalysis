# Created: 2018.05.15
# By: Stephanie R Hilz
# Usage: Create a polished heatmap from Paymann's cpm data for a list of genes of interest 

#source("http://bioconductor.org/biocLite.R")
library(RColorBrewer)
library(goseq)
library(gplots)

## Functions
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
  GO.wall.anno$p.BH.adjust <- p.adjust(GO.wall.anno$over_represented_pvalue, method="BH")
  write.csv(GO.wall.anno,outfile, quote=F)
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

computeZscore <- function(cpmList){
  zscores <- c()
  mean <- mean(cpmList)
  sd <- sd(cpmList)
  for (i in cpmList){
    z <- (i - mean)/sd
    zscores <- append(zscores, z)
  }
  names(zscores) <- names(cpmList)
  return(zscores)
}

## Define file paths and name components
tag <- '20181003_poolker_combined_covariate'
pathtodir <- "/Users/paym01/Documents/UCSF/Bioinformatics/TNIP1/heatmap_pool1_2_covariate/20181003/"
CPMFile <- '20180813_fpool_1_2_covariate_CPMs.csv' 
geneFile <- 'Pool1_2_covariate_861_heatmap_genes.txt'
experiments_OverExpGrouped <- c(
                  'fpool_eGFPorf_unt_24h',
                  'fpool2_eGFPorf_unt_24h',
                  'fpool_eGFPorf_IL17_24h',
                  'fpool2_eGFPorf_IL17_24h',
                 'fpool_eGFPorf_TNF_24h',
                 'fpool2_eGFPorf_TNF_24h',
                 'fpool_TNFAIP3orf_unt_24h',
                 'fpool2_TNFAIP3orf_unt_24h',
                 'fpool_TNFAIP3orf_IL17_24h',
                 'fpool2_TNFAIP3orf_IL17_24h',
                 'fpool_TNFAIP3orf_TNF_24h',
                 'fpool2_TNFAIP3orf_TNF_24h',
                 'fpool_TNIP1orf_unt_24h',
                 'fpool2_TNIP1orf_unt_24h',
                 'fpool_TNIP1orf_IL17_24h',
                 'fpool2_TNIP1orf_IL17_24h',
                 'fpool_TNIP1orf_TNF_24h',
                 'fpool2_TNIP1orf_TNF_24h')

experiments_TreatmentGrouped <- c(
  'fpool2_eGFPorf_IL17_24h',
  'fpool_eGFPorf_IL17_24h',
  'fpool2_eGFPorf_unt_24h',
  'fpool_eGFPorf_unt_24h',
  'fpool2_eGFPorf_TNF_24h',
  'fpool_eGFPorf_TNF_24h',
  'fpool2_TNFAIP3orf_IL17_24h',
  'fpool_TNFAIP3orf_IL17_24h',
  'fpool2_TNFAIP3orf_unt_24h',
  'fpool_TNFAIP3orf_unt_24h',
  'fpool2_TNFAIP3orf_TNF_24h',
  'fpool_TNFAIP3orf_TNF_24h',
  'fpool2_TNIP1orf_IL17_24h',
  'fpool_TNIP1orf_IL17_24h',
  'fpool2_TNIP1orf_unt_24h',
  'fpool_TNIP1orf_unt_24h',
  'fpool2_TNIP1orf_TNF_24h',
  'fpool_TNIP1orf_TNF_24h')

experiments <- experiments_OverExpGrouped # set this to one of the above lists; it is 
# what will then be used to both pick the samples you use and order your samples on
# the x axis.

## Define heatmap color scheme
my_palette <- colorRampPalette(c("black", "pink", "red"))(n = 1000)

## Read in cpm file
cpmData <- read.csv(paste0(pathtodir,CPMFile), header=T, row.names=1)

## Define samples of interest that correspond to experiements defined above
samplesOfInterest <- c()
sampleSuffix <- c(1,2,3)#these are the 3 replicates
for (exp in experiments){
  samplesOfInterest <- append(samplesOfInterest, c(paste0(exp, sampleSuffix)))
}

## Read in genes of interest
genesOfInterest <- scan(paste0(pathtodir,geneFile), what=character())

## Subset the data you will put into the heatmap
genesOfInterest[!genesOfInterest %in% rownames(cpmData)] # these are genes that are not in the CPM file
subset <- cpmData[genesOfInterest,samplesOfInterest]
subset <- subset[complete.cases(subset),]
dim(subset) #this is the final size of your input matrix - may not be all 880

## Set up your heatmap input matrix (can change whether you do log + psuedocount or not)
m <- as.matrix(log(subset+.01,2))

## Normalize first within each sample (can skip this part if don't need to)
pool1Subset <- m[,grepl('pool_',colnames(m))]
pool2Subset <- m[,grepl('pool2_',colnames(m))]
pool1SubsetNormed <- t(apply(pool1Subset,1,computeZscore))
dim(pool1SubsetNormed)
pool2SubsetNormed <- t(apply(pool2Subset,1,computeZscore))
dim(pool2SubsetNormed)
m <- cbind(pool1SubsetNormed, pool2SubsetNormed)
dim(m)
m <- m[,samplesOfInterest] #fixes order back to original

## And now we get to making lots of heatmaps, then doing GO on one of them
## A. Make heatmap where samples appear in original order and are scaled by row Z-score
dist.pear <- function(x) as.dist(1-cor(t(x)))
h <- heatmap.2(m, 
          trace="none", 
          dendrogram="row",
          distfun=dist.pear,
          Colv = FALSE,
          density.info="none",
          cexCol=.3,
          cexRow=.06,
          scale='row',
          col = my_palette)
str(h$rowDendrogram, max.level=5)# checks dendro structure, helpful for

## B. Make heatmap where samples are also clustered (here I am also using average linkage clustering method )
# hclust.ave <- function(x) hclust(x, method="average")
# h <- heatmap.2(m, 
#           trace="none", 
#           dendrogram="row",
#           distfun=dist.pear,
#           hclust=hclust.ave,
#           density.info="none",
#           cexCol=.3,
#           cexRow=.06,
#           scale='row',
#           col = my_palette)
# str(h$rowDendrogram, max.level=6)# checks dendro structure

# C. here we now cut the tree from A based on a desired number of groups, k; for this analysis
#  I am using the k determined by a line Paymann drew to cut his tree
set.seed(1)
k <- 14
hm <- hclust(dist.pear(m)) # we are using the same settings importantly as in the heatmap
groups <- cutree(hm, k=k)
length(groups[which(groups==1)]) #something to run manually to see which k group matches the clusters of the dendro

## Perform GO on each gene group, using all expressed genes as background
for (i in 1:k){
  print(i)
  GOFile <- paste0(pathtodir,'/',tag,'_group_',i, '_GeneOntology.csv')
  genes_group <- as.integer(rownames(cpmData) %in% names(groups[which(groups==i)]))
  names(genes_group) <- rownames(cpmData)
  performGO(genes_group, GOFile)
}

## Create final heatmap in which clusters are annotated
mycol <- c(brewer.pal(8, "Dark2"), brewer.pal(6, "Set2"))
groups.annotation <- factor(groups)
h <- heatmap.2(m, 
               trace="none", 
               dendrogram="row",
               distfun=dist.pear,
               density.info="none",
               Colv = FALSE,
               cexCol=.3,
               cexRow=.06,
               scale='row',
               col = my_palette,
               RowSideColors=mycol[groups.annotation])
legend(x="topright", legend=levels(groups.annotation), col=unique(mycol), pch=15, cex=.9)
str(h$rowDendrogram, max.level=6)# checks dendro structure


# Number of genes per group (group# -> # of genes) for k=14
# 1 -> 23: [NS] single organism signalingresponse to stimulus
# 2 -> 124: immune response term, response to cytokines
# 3 -> 6: [NS] extracellular space, bone regeneration, response to vitamin K
# 4 -> 20: keratinization, cell death
# 5 -> 9: intermediate filament, keratinization
# 6 -> 11: [NS] MAP kinase phosphatase activity, negative regulation of MAPK cascade
# 7 -> 154: response to molecule of bacterial origin, response to lipopolysaccharide
# 8 -> 29: [NS] keratinocyte apoptotic process
# 9 -> 74: antimicrobial humoral response, immune system process 
# 10 -> 107: DNA packaging complex, chromatin, nucleosome
# 11 -> 54: [NS] sequestering of TGFbeta in extracellular matrix
# 12 -> 223: cell adhesion
# 13 -> 12: gastrulation, extracellular matrix organization, angiogenesis
# 14 -> 15 [NS] response to exogenous dsRNA, insulin-like growth factor binding


