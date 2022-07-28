##############################################################################

##  description: Processes mouse datasets including quantile normalization, log2-transformation and other steps.

##  usage: Processdata_Mm(counts)

##  arguments: 
##  counts: The scRNA-seq data matrix with rows labeling genes and columns labeling single cells.

##  value:
##  exp: The scRNA-seq data matrix after preprocessing.

##############################################################################


Processdata_Mm<-
  function(counts){
    
    require("preprocessCore");
    require("AnnotationDbi");
    require("org.Mm.eg.db");
    
    counts <- counts[!duplicated(counts[,1]),];
    rownames(counts) <- counts[,1];
    counts <- counts[,-1];
    ncounts <- normalize.quantiles(as.matrix(counts),copy=FALSE);
    ncounts[ncounts < 1] <- 1;
    ncounts <- log2(ncounts+0.1);
    rownames(ncounts) <- rownames(counts);
    colnames(ncounts) <- colnames(counts);
     
    return(ncounts);
  }

