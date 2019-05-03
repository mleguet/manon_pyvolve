# install.packages('ape')
# install.packages('gtools')
# install.packages('seqinr')
# install.packages('plotrix')

setwd("~/Documents/Manon_Zombi")

library('vhica')
library('ape')
library('gtools')
library('seqinr')
library('plotrix')

gene = dir("gene_ali/")
#TE= "MAX.fasta"

TE= "deug-dsim.fasta"
vc <- read.vhica(gene.fasta = gene, target.fasta = TE, coding = FALSE)
image(vc, "deug-dsim", treefile="Input/tree.nwk", skip.void=TRUE)
plot(vc,'dsim','dyak')


vc <- read.vhica(gene.fasta = gene, target.fasta = TE)
image(vc, "deug-dsim", treefile="Input/tree.nwk", skip.void=TRUE)
plot(vc,'dsim','dyak')
