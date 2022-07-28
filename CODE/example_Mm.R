
count <- read.csv("GSE97390_direct_programming.filtered_normalized_counts.csv", header=TRUE)
count <- count[-1,]
count <- M2H(count)
write.csv(count, "Briggs_trans1.csv")

load("net13Jun12.rda")
load("net17Jan16.rda")
load("hs_km.Rda")

counts <- read.csv("Briggs_trans2.csv")

res0 <- Processdata_Mm(counts)

res1 <- DoIntegPPI(res0, net17Jan16.m)

exp <- res1[["expMC"]]
adj <- res1[["adjMC"]]

## running CompORC.py
curve <- read.csv("Briggs_curve.csv", header=TRUE)

res2 <- E2S(exp, curve)

score = CompOLV(res2[["exp"]], res2[["curve"]], km)

write.csv(score, "Briggs_score.csv")

