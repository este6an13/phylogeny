library(ape)
library(phangorn)
library(seqinr)
library(phytools)

alignment = read.phyDat("alignment.fas", format="fasta", type="DNA")
alignment

dist_matrix = dist.dna(as.DNAbin(alignment))
dist_matrix

UPGMA_tree = upgma(dist_matrix)
plot(UPGMA_tree)

NJ_tree = nj(dist_matrix)
plot(NJ_tree)

UPGMA_tree_root <- root(UPGMA_tree, "Lemur_catta", r=TRUE)
plot(UPGMA_tree_root)

NJ_tree_root <- root(NJ_tree, "Lemur_catta", r=TRUE)
plot(NJ_tree_root)

plot(ladderize(UPGMA_tree_root), align.tip.label=TRUE, cex=0.8)
plot(ladderize(NJ_tree_root), align.tip.label=TRUE, cex=0.8)

writeNexus(UPGMA_tree_root, file="primates_dna_upgma_tree.nex")
writeNexus(NJ_tree_root, file="primates_dna_nj_tree.nex")

boot_trees <- boot.phylo(NJ_tree, as.DNAbin(alignment), FUN = function(xx) nj(dist.dna(xx, model = "JC69")), B = 10000)
boot_trees