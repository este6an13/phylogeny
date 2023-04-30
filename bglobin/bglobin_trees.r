library(ape)
library(phangorn)
library(seqinr)
library(phytools)

alignment = read.phyDat("bglobin.fas", format="fasta", type="DNA")
alignment

dist_matrix = dist.dna(as.DNAbin(alignment))
dist_matrix

layout(matrix(c(1,2)), height=c(1,1.25))
par(mar=c(.1,.1,.1,.1))

UPGMA_tree = upgma(dist_matrix)
plot(UPGMA_tree)

UPGMA_tree

NJ_tree = nj(dist_matrix)
plot(NJ_tree)

NJ_tree

UPGMA_tree_root <- root(UPGMA_tree, "human", r=TRUE)
plot(UPGMA_tree_root)

NJ_tree_root <- root(NJ_tree, "human", r=TRUE)
plot(NJ_tree_root)

plot(ladderize(UPGMA_tree_root), align.tip.label=TRUE)
plot(ladderize(NJ_tree_root), align.tip.label=TRUE)

writeNexus(UPGMA_tree_root, file="bglobin_dna_upgma_tree.nex")
writeNexus(NJ_tree_root, file="bglobin_dna_nj_tree.nex")

boot_trees <- boot.phylo(NJ_tree, as.DNAbin(alignment), FUN = function(xx) nj(dist.dna(xx, model = "JC69")), B = 10000)

parsimony(UPGMA_tree, alignment)
parsimony(UPGMA_tree_root, alignment)

parsimony(NJ_tree, alignment)
parsimony(NJ_tree_root, alignment)

best_parsimony_UPGMA = optim.parsimony(UPGMA_tree, alignment)
plot(best_parsimony_UPGMA)

best_parsimony_NJ = optim.parsimony(NJ_tree, alignment)
plot(best_parsimony_NJ)

parsimony_pratchet <- pratchet(alignment, all=TRUE)
parsimony_pratchet

parsimony_pratchet_root <- root(parsimony_pratchet, "human")
parsimony_pratchet_root

plot(parsimony_pratchet_root, cex=0.6)

strict_100 <- ape::consensus(parsimony_pratchet_root, p=1)
plot(strict_100)

strict_50 <- ape::consensus(parsimony_pratchet_root, p=0.5)
plot(strict_50)

boot_pratchet_trees <- bootstrap.phyDat(alignment, FUN = pratchet, bs = 1000)
boot_pratchet_trees
plot(boot_pratchet_trees)

strict_100 <- ape::consensus(boot_pratchet_trees, p=1)
plot(strict_100)

strict_50 <- ape::consensus(boot_pratchet_trees, p=0.5)
plot(strict_50)

write.tree(strict_50, file="parsimony_pratchet_strict_50.nex")

parsimony_pratchet_strict_50 <- read.tree(file="parsimony_pratchet_strict_50.nex")
parsimony_pratchet_strict_50
parsimony_pratchet_strict_50_root <- root(parsimony_pratchet_strict_50, "human", r=TRUE)

plot(ladderize(parsimony_pratchet_strict_50_root), align.tip.label=TRUE)