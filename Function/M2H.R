##############################################################################

##  description: Gene ID conversion in processing mouse datasets.

##  usage: M2H(mdata)

##  arguments: 
##  mdata: The mouse scRNA-seq data matrix after preprocessing.

##  value:
##  hdata: The human homologene scRNA-seq data matrix.

##############################################################################


M2H<-	
  function(mdata){
    
    library(homologene)
    library(dplyr)
    library(org.Mm.eg.db)
    
    genelist <- mdata[, 1]
    
    id <- homologene(genelist, inTax=10090, outTax=9606)
    id1 <- id[!duplicated(id[,1]),]
    
    colnames(mdata)[1] <- '10090'
    
    x <- dplyr::inner_join(mdata, id1, by='10090')
    hdata <- x[!duplicated(x[,1]),]
    
    hdata[,1] <- id1["9606_ID"]
    
    return(hdata)
  } 

