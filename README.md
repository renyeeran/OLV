# OLV

OLV is a topology-based method to estimate the differentiation potency of single-cells.

OLV is performed as follows:  

(1) Obtainment of expression matrix and adjacency matrix from scRNA-Seq data and PPI network;

(2) Computation of Ollivier Ricci Curvature(ORC) based on adjacency matrix;

(3) Assignment of functional similarities between genes to ORC as weights;

(4) Combination of weighted ORC with expression matrix.

Note: Conversion of human homologous genes is involved when dealing with mouse data.
