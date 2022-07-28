
counts = read.csv("GSE75748_sc_cell_type_ec.csv", header=TRUE)
counts = counts[-1,]

load("net13Jun12.rda")
load("net17Jan16.rda")
load("hs_km.Rda")

res0 <- Processdata_Hs(counts, net17Jan16.m)

res1 <- DoIntegPPI(res0[["exp"]], res0[["adj"]])

exp <- res1[["expMC"]]
adj <- res1[["adjMC"]]

## running CompORC.py
curve = read.csv("Chu1_curve.csv", header=TRUE)

rownames(curve) <- curve[,1]
curve <- curve[,-1]

score = CompOLV(exp,curve,km)

write.csv(score, "Chu1_score.csv")
