##############################################################################

##  description: This function computes the OLV value for each cell with the scRNA-Seq expression matrix, curvature matrix and gene-to-gene functional similarity matrix. 

##  usage: CompOLV(exp,curve,km)

##  arguments: 
##  exp: The output value of function Process_Hs or E2S.
##  curve: The curvature matrix transformed from the output value of CompORC.py or of function E2S.
##  km: The pre-compiled pairwise Kappa similarity matrix on Gene Ontology of human genes.

##  value:
##  OLV_nor: The normalized single-cell potency measure computing by scRNA-seq, Ollivier Curvature and Gene Ontology similarity scores.

##############################################################################


CompOLV<-	
  function(exp,curve,km){

    commonID1.v <- intersect(rownames(km),rownames(curve));
    commonID.v <- intersect(commonID1.v,rownames(exp));
    match(commonID.v,rownames(curve)) -> map1.idx;
    cint <- curve[map1.idx,map1.idx];
    match(commonID.v,rownames(km)) -> map2.idx;
    kmint <- km[map2.idx,map2.idx];
    match(commonID.v,rownames(exp)) -> map3.idx;
    expint <- exp[map3.idx,];
    
    cr <- cint*kmint;
    ORG <- vector();
    for(i in 1:ncol(cr)){
      ORG[i] <- sum(cr[i,]);
    }
    OLV <- as.vector(t(expint)%*%ORG);
    OLV_nor <- OLV/max(OLV);  #normalization
    return(OLV_nor);
  }

