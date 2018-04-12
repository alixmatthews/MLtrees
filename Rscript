# Building easy maximum likelihood trees with R
# example using CAD data from ants from GenBank

setwd("...")
library(ape)
library(phangorn)

# import aligned fasta file
cad <- read.dna("CAD_aligned.fas", "f")

# set up data
mlCad <- as.phyDat(cad) # changing fasta file into a phyDat class (need for later)
tr.Cad <- nj(dist.dna(cad)) # neighbor-joining tree
m0.Cad <- pml(tr.Cad, mlCad, k = 4) # use as base tree

# modelTest for best nucleotide substitution model
mltest.Cad <- modelTest(mlCad) # test for best model from your fasta file
mltest.Cad # results as table
which.min(mltest.Cad$AIC) # which is the lowest AIC?
which.min(mltest.Cad$AICc) # which is the lowest AICc? (small sample size)

# in this case, HKY, build tree under best model
m1.hky.Cad <- optim.pml(m0.Cad, model="HKY", optGamma = FALSE, optInv = FALSE) # build off of base tree; if you get +G or +I in best model, put TRUE instead
plot(m1.hky.Cad) # unrooted, best model tree

# optimize tree toplogies with best model
m2.hky.Cad <- optim.pml(m1.hky.Cad, optNni = TRUE)
plot(m2.hky.Cad) # unrooted, optimized best model tree

# root tree
rooted.Cad <- root(m2.hky.Cad$tree, "...") # your outgroup here
plot(rooted.Cad)

# easy bootstraps
ctr <- pml.control(trace=0) # helps "hide"/control results as they're running
m3.hky.Cad.bs <- bootstrap.pml(m2.hky.Cad, optNni = TRUE, bs = 100, control = ctr) # boostraps (100)
bs.rooted.Cad <- plotBS(rooted.Cad, m3.hky.Cad.bs, p = 80, type = "phylogram") # p = ... for your cut-off for bootstrap values
plot(bs.rooted.Cad)
write.tree(bs.rooted.Cad, "Cad_final.tre") # you can open this final .tre file in FigTree now

############################################
#### another bootstrap option ####
# bootstraps
# ctr <- pml.control(trace=0) # helps "hide"/control results as they're running
# m3.hky.Cad.bs <- bootstrap.pml(m2.hky.Cad, optNni = TRUE, bs = 100, control = ctr) # boostraps (100)
# prop.clades.bp.Cad <- prop.clades(m2.hky.Cad$tree, m3.hky.Cad.bs) # connects ($tree) to your bootstraps
# plot(m2.hky.Cad$tree)
# drawSupportOnEdges(prop.clades.bp.Cad)
# 
# # root and redraw
# rooted.Cad <- root(m2.hky.Cad$tree, "...") # your outgroup here
# plot(rooted.Cad)
# drawSupportOnEdges(prop.clades.bp.Cad) # final tree
# 
# # save to a .tre file (and open/edit in FigTree)
# write.tree(rooted.Cad, "Cad_final.tre")
############################################
