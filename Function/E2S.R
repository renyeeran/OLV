##############################################################################

##  description: Gene ID conversion in processing mouse datasets.

##  usage: E2S(exp,curve)

##  arguments: 
##  exp: The output value of function Processdata.
##  curve: The curvature matrix transformed from the output value of CompORC.py

##  value:
##  exp: The scRNA-seq data matrix after gene ID conversion.
##  curve: The curvature matrix after gene ID conversion.

##############################################################################


E2S<-	
  function(exp,curve){
    require("preprocessCore");
    require("AnnotationDbi");
    require("org.Hs.eg.db");

    r = as.character(exp[,1])
    geneIDselect1.0 <- select(org.Hs.eg.db, keys=r, columns="SYMBOL", keytype="ENTREZID");
    geneIDselect1.1 <- geneIDselect1.0[!duplicated(geneIDselect1.0[,2]),]
    geneIDselect1.2 <- geneIDselect1.1[!is.na(geneIDselect1.1[,2]),]
    
    #e = exp[,-1]
    e = e[!duplicated(geneIDselect1.0[,2]),]
    e = e[!is.na(geneIDselect1.1[,2]),]
    rownames(e) = geneIDselect1.2[,2]
    
    c = curve[,-1]
    c = c[!duplicated(geneIDselect1.0[,2]),!duplicated(geneIDselect1.0[,2])]
    c = c[!is.na(geneIDselect1.1[,2]),!is.na(geneIDselect1.1[,2])]
    rownames(c) = geneIDselect1.2[,2]
    colnames(c) = geneIDselect1.2[,2]
    
    return(list(exp=e,curve=c))
  }

