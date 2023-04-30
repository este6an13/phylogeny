library(ape)
library(phangorn)
library(seqinr)
library(phytools)

alignment = read.phyDat("alignment.fas", format="fasta", type="DNA")
alignment

dist_matrix = dist.dna(as.DNAbin(alignment))
dist_matrix

parsimony_pratchet <- pratchet(alignment, all=TRUE)
parsimony_pratchet

boot_pratchet_trees <- bootstrap.phyDat(alignment, FUN = pratchet, bs = 100)
boot_pratchet_trees

strict_50 <- ape::consensus(boot_pratchet_trees, p=0.5)

write.tree(strict_50, file="primates_dna_parsimony_boot_pratchet_strict_50.nex")

strict_50_root <- root(strict_50, "Lemur_catta", r=TRUE)

plot(ladderize(strict_50_root), align.tip.label=TRUE, cex=0.8)