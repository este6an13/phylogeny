library(ape)
library(phangorn)
library(phytools)
library(seqinr)

alignment <- read.phyDat("alignment.fas", format="fasta", type="DNA")

random_tree <- rtree(n=12, tip.label=names(alignment))

verosim_random_tree_obj <- pml(tree=random_tree, data=alignment)

verosim_random_tree <- verosim_random_tree_obj$tree

model_test <- modelTest(alignment, verosim_random_tree)

best_model_tree_obj <- optim.pml(object=verosim_random_tree_obj, model="GTR", optGamma=TRUE, optInv=TRUE, rearrangement="ratchet")

best_model_tree <- best_model_tree_obj$tree

bootstrap <- bootstrap.pml(verosim_random_tree_obj, bs=1000, optNni=TRUE)

plotBS(tree=ladderize(root(best_model_tree, outgroup='Homo_sapiens'), FALSE), 
       align.tip.label=TRUE, cex=0.8)
