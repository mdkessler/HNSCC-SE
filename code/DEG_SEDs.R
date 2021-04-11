# Exploratory analysis and figure generation for relationships between
# (differential) gene expression and (differential) SEDs (super enhancer
# domains). Aims to prepare data to allow for exploration of
# gene expression patterns for genes around differential
# SEDs, association between differential expression and
# differential SEs, etc.

######################################################################
# Set up working environment
######################################################################

# packages
library(DESeq2) # need to load this first to import rds below
library(rtracklayer)
library(dplyr)
library(plyranges)
library(GenomicRanges)
library(differential.coverage)
library(stringr)

#load dSEDs
SEDs<-readRDS('../data/SE_by_ChipSeq.rds')
# calculate the midpoint of each SE
# this will be used to measure distance to nearest gene
SEDmids <- makeGRangesFromDataFrame(SEDs,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames",
                                                     "seqname",
                                                     "chromosome",
                                                     "chrom",
                                                     "chr",
                                                     "chromosome_name",
                                                     "seqid"),
                                    start.field="start",
                                    end.field=c("end", "stop"),
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)# copy object as genomic ranges
# Note: in data, common_SE is how Sasha calls SED (super enhancer domain)
start(SEDmids) <- (start(SEDmids) + end(SEDmids)) / 2 # reassign start to midpoint
width(SEDmids) <- 1 # set width to 0, which should make end == start
# set all strands to *, as SEDs seem to be on both strands
strand(SEDmids) <- "*"

# Use Sasha's differential.coverage pacakge to get hg19 gene annotations
#gene.anno <- differential.coverage::get.Known.Gene.List(genome.annotation.id = 'gencode19')
gene.anno <- differential.coverage::gencode_hs19_genes
# Note: genes exist in chrX, chrY, chrM. We will remove chrY and chrM
# remove duplicate genes 
gene_freq <- table(gene.anno$ensembl)
gene.dups <- names(which(gene_freq > 1))
# remove duplicates and chrY/chrM
gene.anno.u <- gene.anno %>%
  filter((!ensembl %in% gene.dups) &
           seqnames != "chrY" &
           seqnames != "chrM")

# set the ranges of each gene to only be the start site (presumably transcription start site)
# note: account for strand - if strand == "-", set start to end
start(gene.anno.u) <- ifelse(as.character(strand(gene.anno.u)) == '-',
                             end(gene.anno.u),
                             start(gene.anno.u))
width(gene.anno.u) <- 1

# use distanceToNearest function to calculate the closest SED for each gene
gene_SED_hits <- distanceToNearest(gene.anno.u, SEDmids, ignore.strand=TRUE)
# error check that all entries had a nearest SED hit
stopifnot(length(queryHits(gene_SED_hits)) == length(gene.anno.u))
# subset out hits
genehits <- gene.anno.u[queryHits(gene_SED_hits)]
SEDhits <- SEDmids[subjectHits(gene_SED_hits)]
dists <- mcols(gene_SED_hits)$distance # conveniently access below

# adjust dist direction (sign) depending on whether gene TSS is upstream or
# downstream of SED location
gene_strand_minus <- as.character(strand(genehits))=='-'
SED_bigger <- start(genehits) < start(SEDhits)
# update dists - dists including direction from SED to gene TSS
# if minus and SED bigger, than its upstream, and +dist to gene
# if plus and SED smaller (FALSE == FALSE in comparison below), than its upstream, and +dist to gene
## otherwise, its the reverse
dists <- ifelse(gene_strand_minus == SED_bigger, dists, -dists)

# make df with results
SED.df <- data.frame(
  distance = dists,
  SE_ID = SEDhits$ID,
  SE_logFC = SEDhits$logFC,
  SE_logCPM = SEDhits$logCPM,
  SE_F = SEDhits$F,
  SE_PValue = SEDhits$PValue,
  SE_FDR = SEDhits$FDR,
  stringsAsFactors=FALSE)

# add genes as rownames
rownames(SED.df) <- genehits$ensembl

# get differential expression results for 77 samples from RNA-seq
DEG <- readRDS('../data/T_to_N_DE_results_hg19_ensembl_77_samples.rds')
# rename columns
colnames(DEG) <- str_c("GENE_",colnames(DEG))
colnames(DEG)[2] <- c("GENE_log2FC")
  
# find gene names with differential expression values that are also
# in the SED table I just made
genes.int <- intersect(rownames(DEG), rownames(SED.df))

# combine DEG with SED
DEG.SED.df <- as.data.frame(cbind(SED.df[genes.int,], DEG[genes.int,]))
saveRDS(DEG.SED.df,'../data/DEG.SED.df.rds')
