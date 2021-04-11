# The following analysis will be centered around the TFs
# that were identified as enriched around HNSCC dSEDs in our
# Cistrome analysis.
# It will evaluate the following:
# 1) The trancriptional levels of the TFs and their target
# genes (as annotated in the trrust database)
# 2) The transcriptional levels of these TFs + targets in
# TCGA, split by tissue and by HPV status (clustering/network)
# 3) When looking at clustering in the previous analyses, 
# consider clustering by only TFs, TFs + targets, and TFs + 
# targets + the targets of TFs that are targets of enriched TFs
# 4) Evaulate transcriptional changes post JQ1 treatment
# within the aforementioned analysis frameworks
# 5) See what pathways are upregulated within the TFs + target
# genes, and TFs + targets + target TF targets

######################################################################
# Set up working environment
######################################################################

# packages
library(dplyr)
library(org.Hs.eg.db)
library(stringr)
library(RCurl)
library(fgsea)
library(igraph)
library(paletteer)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(qgraph)
library(gplots)
library(plot3D)
library(ggrepel)
library(tidyverse)
library(pals)
library(sva)
library(ClassDiscovery)
library(simpleaffy)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(ggupset)

# set dir
setwd("/Users/michaelkessler/Desktop/Workspace/POSTDOC/enhancers/HNSCC-SEs-genes/MDK")

# source functions
source("TFTargetJaccard.R")

# load HNSCC data
df.tn <- as.data.frame(readRDS("DEG.SED.df.rds"))
# Load JQ1 data
df.047 <- readRDS("../RNAseq-gencode19/JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
df.090 <- readRDS("../RNAseq-gencode19/JQ1_vs_DMSO_DE_results_O90_JQ1.rds")

####################################################################
# Preprocessing
####################################################################

# rename cols for 047 and 090 data
names(df.047) <- str_c(names(df.047), "_047")
names(df.090) <- str_c(names(df.090), "_090")

# match gene names just in case the order differs (even though I believe the gene names are in the same order)
mtch_047 <- match(rownames(df.tn), rownames(df.047))
mtch_090 <- match(rownames(df.tn), rownames(df.090))

# combine dfs
df <- as.data.frame(cbind(df.tn, df.047[mtch_047,], df.090[mtch_047,]))

# annotate gene symbols
symbols <- mapIds(org.Hs.eg.db, rownames(df), 'SYMBOL', 'ENSEMBL')
entrez <- mapIds(org.Hs.eg.db, rownames(df), 'ENTREZID', 'ENSEMBL')

# Now add gene symbols to df
df$gene_symbols <- symbols
df$gene_entrez <- entrez
# remove any rows with missing gene symbols or entrez IDs
df <- subset(df, !is.na(gene_symbols) & !is.na(gene_entrez))
# remove duplicates
df <- df[!duplicated(df$gene_symbols),]

####################################################################
# Make custom gene lists representing TF target genes
####################################################################

# read in TF enrichment data
TF.df <- readRDS("TF.enrichment_analysis.rds") # get TFs whose target genes you want to analyze for overlap with genes comprising leading edges of enriched pathways

# select the TFs for whose target genes you will look for enrichments
# these will be those that are significantly (padj <= 0.05)
# enriched in tumor dSED peaks (i.e. delta > 0)
TFs <- unique(subset(TF.df, padj <= 0.05 & delta > 0)$TF)

# get TF gene sets
# TF gene sets from https://www.grnpedia.org/trrust/downloadnetwork.php
trustt_URL <- "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
myfile <- getURL(trustt_URL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
TF_paths.df <- read.table(textConnection(myfile), header=F)

# make list with TF as name and gene sets as values of character vectors
TF_gene_list <- list()
for (TF in TFs){
  TF_genes <- unique(as.character(subset(TF_paths.df,
      V1 == TF)$V2))
  #TF_gene_regulation <- as.character(subset(TF_paths.df,
      # V1 == TF)$V3)
  TF_gene_list[[TF]] <- data.frame(TF_genes = TF_genes)
      # TF_gene_regulation = TF_gene_regulation)
}

# now expand this TF_gene_list to include any genes that are
# targets of the TFs that are themselves targets of the
# enriched TFs
expanded_TF_gene_list <- TF_gene_list

# first initialize an expanded TF_gene_list from the TF_gene_list with only the TF from our enrichment analysis
# now expand this list by seeing whether any of the genes in the TF lists are TFs themselves that are not already in the TF list. If they are, add them to the TF list, along with the genes they target.
TF_in_trustt <- unique(as.character(TF_paths.df$V1))
for (TF in names(TF_gene_list)){
  TF_df <- TF_gene_list[[TF]]
  genes_are_TFs <- TF_df$TF_genes %in% TF_in_trustt &
    !TF_df$TF_genes %in% names(TF_gene_list)
  if (sum(genes_are_TFs) > 0){ # at least one gene is annotated TF
    for (TF2 in TF_df$TF_genes[genes_are_TFs]){
      TF2_genes <- as.character(subset(TF_paths.df,
            V1 == TF2)$V2)
      # TF2_gene_regulation <- as.character(subset(TF_paths.df,
            # V1 == TF2)$V3)
      expanded_TF_gene_list[[TF2]] <- data.frame(TF_genes = TF2_genes)
            # TF_gene_regulation = TF2_gene_regulation)
    }
  }
}

# now make a list of only TFs
# that is, TF as list names, and their values are character
# vectors of any target genes that are also TFs. So any
# non-TF target genes will be removed from this list.
only_TFs_list <- lapply(TF_gene_list, function(x){
    hit <- x$TF_genes %in% names(TF_gene_list) 
    return(data.frame(TF_genes = as.character(x$TF_genes[hit])))
  }
  )

#############################################################
# Prepare nodes and edges for network plots
#############################################################

# make edges df using dplyr bind_rows function
# Note: the warnings produced by the following command are fine - they just corce factors into character vectors as thigns are combined
edges <- dplyr::bind_rows(only_TFs_list, .id = "column_label")
# rename cols
names(edges) <- c("label", "pair")

# get rid of self connected nodes
edges <- edges[!edges$label == edges$pair,]

# node should be all TFs and Gene
nodes <- as.data.frame(unique(c(names(only_TFs_list),
                        as.character(edges$label),
                        as.character(edges$pair))))

# name col
names(nodes) <- "label"

# now add meta data to nodes df
# first whether the node is a TF or a gene
nodes$status <- "Gene"
nodes$status[nodes$label %in% names(TF_gene_list)] <- "TF"

# add transcriptomic data as meta data
mtch <- match(nodes$label, df$gene_symbols)
nodes$log2FC <- df$GENE_log2FC[mtch]
nodes$log2FC_JQ1047 <- df$log2FoldChange_047[mtch]
nodes$log2FC_JQ1090 <- df$log2FoldChange_090[mtch]
# mean of cell line changes post JQ1
nodes$log2FC_JQ1 <- rowMeans(cbind(nodes$log2FC_JQ1047, nodes$log2FC_JQ1090), na.rm=TRUE)

# add pvals
nodes$padj <- df$GENE_padj[mtch]
nodes$padj_JQ1 <- pmin(df$padj_047[mtch], df$padj_090[mtch])

# generate a color ramp over log2FC
# Create a color palette
nColor <- 14
# colors <- paletteer_c(palette = "viridis::viridis", n = nColor)
#colors <- brewer.pal(n = nColor, name = "Blues")
pal1 <- colorRampPalette(c("blue", "purple",
                          "white", "orange", "red"))
pal2 <- colorRampPalette(c("darkgreen",
                           "white", "hotpink"))

colors1 <- pal1(nColor)
colors2 <- pal2(nColor)

# bin the data but account for extreme
cuts <- apply(nodes[,3:6], 2, cut,
              c(-Inf,seq(-3, 3, 0.5), Inf), labels=1:14)
# set colors using these bins
nodes$log2FC_col <- colors1[as.numeric(cuts[,1])]
nodes$log2FC_JQ1047_col <- colors2[as.numeric(cuts[,2])]
nodes$log2FC_JQ1090_col <- colors2[as.numeric(cuts[,3])]
nodes$log2FC_JQ1_col <- colors2[as.numeric(cuts[,4])]
# now assign NA values to a color
nodes$log2FC_col[which(is.na(nodes$log2FC))] <- "darkgrey"
nodes$log2FC_JQ1047_col[which(is.na(nodes$log2FC_JQ1047))] <- "darkgrey"
nodes$log2FC_JQ1090_col[which(is.na(nodes$log2FC_JQ1090))] <- "darkgrey"
nodes$log2FC_JQ1_col[which(is.na(nodes$log2FC_JQ1))] <- "darkgrey"

#############################################################
# Make network plots
#############################################################

# plot with iRanges package

# make graphs with only TFs and with the coords determined
# by dist and mds

# first make the graph from the data I prepared
ig <- graph_from_data_frame(d = edges,
          vertices = nodes, directed = T)

coords <- layout_in_circle(ig)
#coords <- layout_with_kk(ig)
#coords <- layout_with_fr(ig)# * 0.2#, dim = 3,
            #weights = rep(0.01, nrow(edges)))#, grid = "nogrid")
#coords <- layout_with_drl(ig) * 3#, options=list(edge.cut=0))

pdf("TF_network.pdf", height = 15, width = 15)
plot(ig,
     layout = coords,
     rescale = T,
     #vertex.label = ifelse(V(ig)$padj <= 0.05, V(ig)$name, NA),
     vertex.label = V(ig)$name,
     #vertex.label.cex = ifelse(V(ig)$status == "TF", 1.5, NA),
     vertex.label.cex = 2,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_col,
     #vertex.size = ifelse(V(ig)$status == "TF", 6, 1.2),
     vertex.size = 10,
     edge.width = 1,
     edge.color = "grey60",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)

text(x = coords[,1]*1.2,
     y = coords[,2]*1.2,
     labels = ifelse(V(ig)$padj <= 0.05, "*", NA),
     col = "red",
     cex = 7)

dev.off()

# make another network colored by post JQ1 (mean across both cell lines)

pdf("TF_network.JQ1.pdf", height = 15, width = 15)
plot(ig,
     layout = coords,
     rescale = T,
     #vertex.label = ifelse(V(ig)$padj <= 0.05, V(ig)$name, NA),
     vertex.label = V(ig)$name,
     #vertex.label.cex = ifelse(V(ig)$status == "TF", 1.5, NA),
     vertex.label.cex = 2,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_JQ1_col,
     #vertex.size = ifelse(V(ig)$status == "TF", 6, 1.2),
     vertex.size = 10,
     edge.width = 1,
     edge.color = "grey60",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)

text(x = coords[,1]*1.2,
     y = coords[,2]*1.2,
     labels = ifelse(V(ig)$padj_JQ1 <= 0.05, "*", NA),
     col = "red",
     cex = 7)

dev.off()

# make color legends
pdf("TF_network.log2_legend.pdf")
plot(rep(1, 16), 1:16, col = c("white", colors1, "white"),
     pch = 19, cex = 3)
text(x = rep(1, 15), 1:15,
     labels = as.character(c(-Inf,seq(-3, 3, 0.5), Inf)),
     pos  = 4, offset = 2, cex = 1.5)
     #adj = c(0, -1), cex = 0.5)
dev.off()

pdf("TF_network.log2_JQ1_legend.pdf")
plot(rep(1, 16), 1:16, col = c("white", colors2, "white"),
     pch = 19, cex = 3)
text(x = rep(1, 15), 1:15,
     labels = as.character(c(-Inf,seq(-3, 3, 0.5), Inf)),
     pos  = 4, offset = 2, cex = 1.5)
dev.off()

# plot TF distance matrix
pdf("TF.jaccard_distance_matrix.heatmap.pdf")
par(mar = c(6,6,2,2))
image(1:ncol(jac_mat), 1:ncol(jac_mat),
      jac_mat, axes = FALSE, xlab="", ylab="", )

axis(1, 1:ncol(jac_mat), colnames(jac_mat),
                  cex.axis = 1, las=3)
axis(2, 1:ncol(jac_mat), colnames(jac_mat),
     cex.axis = 1, las=2)

dev.off()
#######################################################################
# Make plot of TFs relative to one another based on target gene overlap
#######################################################################

# calculate x and y coords for TFs using jaccard dist and mds
jac_mat <- TFTargetJaccard(TF_gene_list, dist = T)
jac_mat[is.nan(jac_mat)] <- 1 # deal with Nan
# now perform mds
mds <- cmdscale(d = jac_mat, eig = TRUE, k = 2) # k is the number of dim
coords <- as.data.frame(mds$points)
# x <- mds$points[,1]
# y <- mds$points[,2]
# z <- mds$points[,3]
# coords <- cbind(x, y, z)
# # make sure row order is correct
# coords <- coords[match(nodes$label, rownames(coords)),]
# # rename cols
# colnames(coords) <- c("V1", "V2")
# error check
stopifnot(rownames(mds$points) == names(V(ig)))
coords$cols <- V(ig)$log2FC_col
coords$JQ1_cols <- V(ig)$log2FC_JQ1_col

# 2d plot of mds coords
pdf("TF_mds.2D.pdf")
ggplot(coords, aes(x= V1, y = V2)) + 
  geom_point(color = coords$col, size = 4) +
  geom_text_repel(aes(label = rownames(coords)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   size = 3,
                   segment.color = 'grey50') +
  theme_classic() +
  theme(#axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
dev.off()

# repeat plot using JQ1 colors
# 2d plot of mds coords
pdf("TF_mds.JQ1.2D.pdf")
ggplot(coords, aes(x= V1, y = V2)) + 
  geom_point(color = coords$JQ1_cols, size = 4) +
  geom_text_repel(aes(label = rownames(coords)),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  size = 3,
                  segment.color = 'grey50') +
  theme_classic() +
  theme(#axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
dev.off()

####################################################################
##### Run customGSEA analysis (see customGSEA script and url below)
####################################################################

# Here you will run GSEA on ranked gene list and post JQ1 gene list
# Use hallmark gene sets and oncogenic gene sets

# first make ranked gene list - this will be subset below
# for each analysis
ranked_gene_list = df$GENE_log2FC / df$GENE_lfcSE # use logfoldchange normalized by standard eror, which should recreate the gene level stat

# set names to gene symbols
names(ranked_gene_list) = df$gene_entrez
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# now, get go term gene pathways
# pathways downloaded from msigdb at http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C5

# pathway files
c5bp_file <- "c5.bp.v6.2.symbols.gmt.txt"
c6onc_file <- "c6.all.v7.0.symbols.gmt.txt"
h_all_file <- "h.all.v7.1.symbols.gmt"

# pathways
myPaths <- fgsea::gmtPathways(h_all_file)

# run fgsea
fgseaRes <- fgsea(pathways = myPaths, 
                  stats = ranked_gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 1000)

# collapsedPathways <- collapsePathways(fgseaRes[order(fgseaRes$pval),],
#                                       myPaths.filt,
#                                       ranked_gene_list.TF_target_genes)
# 
# mainPathways <- fgseaRes[fgseaRes$pathway %in% collapsedPathways$mainPathways,]
# upPaths <- mainPathways[mainPathways$NES > 0,]
# upPaths <- upPaths[head(order(upPaths$NES, decreasing = T),10),]
# downPaths <- mainPathways[mainPathways$NES < 0,]
# downPaths <- downPaths[head(order(downPaths$NES, decreasing = F),10),]
# topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

kk2 <- gseKEGG(geneList     = ranked_gene_list,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

#hall <- read.gmt(h_all_file)
egmt2 <- GSEA(ranked_gene_list, TERM2GENE=hall, verbose=FALSE)
dotplot(egmt2, showCategory=100, x = 'NES', size = 'GeneRatio')
heatplot(egmt2, foldChange=ranked_gene_list)
upsetplot(egmt2, foldChange=ranked_gene_list)

sigPaths <- subset(fgseaRes, padj <= 0.01)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

pdf("GSEA.H-all.pdf", height=6, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths[topPaths$pathway],
              ranked_gene_list, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# now repeat but with oncogenic pathways

# pathways
myPaths <- fgsea::gmtPathways(c6onc_file)

# run fgsea
fgseaRes <- fgsea(pathways = myPaths, 
                  stats = ranked_gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 1000)

# collapsedPathways <- collapsePathways(fgseaRes[order(fgseaRes$pval),],
#                                       myPaths.filt,
#                                       ranked_gene_list.TF_target_genes)
# 
# mainPathways <- fgseaRes[fgseaRes$pathway %in% collapsedPathways$mainPathways,]
# upPaths <- mainPathways[mainPathways$NES > 0,]
# upPaths <- upPaths[head(order(upPaths$NES, decreasing = T),10),]
# downPaths <- mainPathways[mainPathways$NES < 0,]
# downPaths <- downPaths[head(order(downPaths$NES, decreasing = F),10),]
# topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

sigPaths <- subset(fgseaRes, padj <= 0.01)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

pdf("GSEA.c6Onc.pdf", height=6, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths[topPaths$pathway],
              ranked_gene_list, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# JQ1 treatment

# pathways
myPaths <- fgsea::gmtPathways(h_all_file)

# calculate means across both cell lines treated with JQ1
df$log2FoldChange_JQ1 <- rowMeans(cbind(
  df$log2FoldChange_047, df$log2FoldChange_090), na.rm=TRUE)
df$lfcSE_JQ1 <- rowMeans(cbind(
  df$lfcSE_047, df$lfcSE_090), na.rm=TRUE)

# make new ranked gene list using JQ1 gene expression levels
ranked_gene_list = df$log2FoldChange_JQ1 / df$lfcSE_JQ1 # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_symbols
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# run fgsea
fgseaRes <- fgsea(pathways = myPaths, 
                  stats = ranked_gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 1000)
sigPaths <- subset(fgseaRes, padj <= 0.05)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

# output figure showing positive pathways enrichments
pdf("GSEA.JQ1.H-all.pdf", height=13, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths[topPaths$pathway],
              ranked_gene_list, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# oncogenic pathways JQ1
# pathways
myPaths <- fgsea::gmtPathways(c6onc_file)

# calculate means across both cell lines treated with JQ1
df$log2FoldChange_JQ1 <- rowMeans(cbind(
  df$log2FoldChange_047, df$log2FoldChange_090), na.rm=TRUE)
df$lfcSE_JQ1 <- rowMeans(cbind(
  df$lfcSE_047, df$lfcSE_090), na.rm=TRUE)

# make new ranked gene list using JQ1 gene expression levels
ranked_gene_list = df$log2FoldChange_JQ1 / df$lfcSE_JQ1 # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_symbols
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# run fgsea
fgseaRes <- fgsea(pathways = myPaths, 
                  stats = ranked_gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 1000)
sigPaths <- subset(fgseaRes, padj <= 0.05)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

# output figure showing positive pathways enrichments
pdf("GSEA.JQ1.c6Onc.pdf", height=13, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths[topPaths$pathway],
              ranked_gene_list, fgseaRes, 
              gseaParam = 0.5)
dev.off()


###################################################################
# GSEA using prepare analysis workflow
# https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#heatmap-like-functional-classification
###################################################################

# make ranked gene list
ranked_gene_list = df$GENE_log2FC / df$GENE_lfcSE # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_entrez
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# try enrichment analysis (overlap analysis) using disease gene networks
edoDO <- gseDO(ranked_gene_list, nPerm = 1000)
edoNCG <- gseNCG(ranked_gene_list, nPerm = 1000)
edoDGN <- gseDGN(ranked_gene_list, nPerm = 1000)
pltDO <- dotplot(edoDO, showCategory=20) + ggtitle("dotplot for GSE - DO")
pltNCG <- dotplot(edoNCG, showCategory=20) + ggtitle("dotplot for GSE - NCG")
pltDGN <- dotplot(edoDGN, showCategory=20) + ggtitle("dotplot for GSE - DGN")
pltDGN

# now JQ1
# make new ranked gene list using JQ1 gene expression levels
ranked_gene_list = df$log2FoldChange_JQ1 / df$lfcSE_JQ1 # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_entrez
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

edoDO <- gseDO(ranked_gene_list, nPerm = 1000)
edoNCG <- gseNCG(ranked_gene_list, nPerm = 1000)
edoDGN <- gseDGN(ranked_gene_list, nPerm = 1000)
pltDO <- dotplot(edoDO, showCategory=20) + ggtitle("dotplot for GSE - DO")
pltNCG <- dotplot(edoNCG, showCategory=20) + ggtitle("dotplot for GSE - NCG")
pltDGN <- dotplot(edoDGN, showCategory=20) + ggtitle("dotplot for GSE - DGN")

# plot gene enrichment networks
edoNCGx <- setReadable(edoNCG, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edoNCGx, foldChange=ranked_gene_list)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edoNCGx, categorySize="pvalue", foldChange=ranked_gene_list)
p2
p3 <- cnetplot(edoNCGx, foldChange=ranked_gene_list, circular = TRUE, colorEdge = TRUE)
p3
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

###################################################################
# ORA conditioned on SEs
###################################################################
SE.df <- subset(df, abs(distance) < 1.5e6 &
                       SE_logFC > 0 &
                       SE_PValue <= 0.05)
# make ranked gene list
ranked_gene_list = SE.df$GENE_log2FC / SE.df$GENE_lfcSE # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = SE.df$gene_entrez
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# try enrichment analysis (overlap analysis) using disease gene networks
edoDO <- enrichDO(names(ranked_gene_list))
edoNCG <- enrichNCG(names(ranked_gene_list))
edoDGN <- enrichDGN(names(ranked_gene_list))
# barplot(edo, showCategory=20)
# edo2 <- gseNCG(ranked_gene_list, nPerm=10000)
# edo2 <- gseDO(ranked_gene_list, nPerm=10000)
# edo2 <- gseDGN(ranked_gene_list, nPerm=10000)
pltDO <- dotplot(edoDO, showCategory=30) + ggtitle("dotplot for ORA - DO")
pltNCG <- dotplot(edoNCG, showCategory=30) + ggtitle("dotplot for ORA - NCG")
pltDGN <- dotplot(edoDGN, showCategory=1000) + ggtitle("dotplot for ORA - DGN")
pltDGN
pdf("test.DGN.pdf", width = 30, height = 30)
pltDGN
dev.off()


# # test plot
# # remove NAs
# View(edo2@result)
# edo2@result <- edo2@result[!is.na(edo2@result$Description),]
# # reorder factors in order of NES score
# edo2@result$Description <- factor(edo2@result$Description, levels = edo2@result$Description[order(edo2@result$NES)])
# ggplot(edo2@result,
#        aes(x = NES, y = Description, size = )) +
#     geom_point()
# 
# length(ranked_gene_list)

###################################################################
# Make heatmaps using TFs + TF_Target_genes
# Do this using our HNSCC RNA-seq data and TCGA data
###################################################################

# read in HNSCC data
RNAseq.all <- readRDS('../RNAseq-gencode19/txi_gencode_19_ensembl_DGAY_77_samples.rds')
# set RNAseq dataframe to abundance values
HNSCC.RNAseq <- RNAseq.all$abundance # set 
# clean up RNAseq colnames
colnames(HNSCC.RNAseq) = gsub('_trimed','',colnames(HNSCC.RNAseq))
colnames(HNSCC.RNAseq) = gsub('_2','',colnames(HNSCC.RNAseq))
# replace ensembl genes names with gene symbols
symbols <- mapIds(org.Hs.eg.db, rownames(HNSCC.RNAseq), 'SYMBOL', 'ENSEMBL')
# Now set rownames to gene symbol
rownames(HNSCC.RNAseq) <- symbols
# remove any rows with missing gene symbols
HNSCC.RNAseq <- subset(HNSCC.RNAseq, !is.na(rownames(HNSCC.RNAseq)))
# filter HNSCC data (77 samples) data to exclude genes that aren't expressed or that have minimal fold change
HNSCC.RNAseq <- HNSCC.RNAseq[apply(HNSCC.RNAseq,1,max) > 0, ]
HNSCC.RNAseq <- HNSCC.RNAseq[apply(HNSCC.RNAseq,1,function(x){max(x)-min(x)})>1,]
# remove any duplicated genes
HNSCC.RNAseq <- HNSCC.RNAseq[-which(rownames(HNSCC.RNAseq) %in% rownames(HNSCC.RNAseq)[duplicated(rownames(HNSCC.RNAseq))]),]
# read in labels for tumor vs normal HNSCC samples
HNSCC_anno <- read.table(file='../RNAseq-gencode19/DGAY_77_sampleID_mapping.txt',
                         header = TRUE, sep = "\t", stringsAsFactors = F)
# build named vector where ID is name and tumor class is label
HNSCC_labels <- HNSCC_anno %>% filter(class %in% c('T','N')) %>% column_to_rownames('ID') %>% dplyr::select(class)
HNSCC_labels <- rbind(HNSCC_labels, 'T') # add manually as its missing from meta data
rownames(HNSCC_labels)[nrow(HNSCC_labels)] <- 'DGay12-25742'

# read in TCGA data
expTCGA <- readRDS("/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI/TCGA_RDS/TCGA.legacy.expression.rds")
metaTCGA <- readRDS("/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI/TCGA_RDS/TCGA.legacy.meta.rds")

# remove duplicated genes
expTCGA <- expTCGA[-which(rownames(expTCGA) %in% rownames(expTCGA)[duplicated(rownames(expTCGA))]),]

# identify the intersect of genes, and subset and merges datasets
genes.intersect <- intersect(rownames(HNSCC.RNAseq), rownames(expTCGA))

# subset genes to intersected genes in both HNSCC and TCGA data sets
expTCGA.its <- expTCGA[genes.intersect,]
HNSCC.RNAseq.its <-  HNSCC.RNAseq[genes.intersect,]



# SET HEATMAP COLORS

# HNSCC col

mycols <- sample(polychrome())

# make vector of cols with samples as names - these cols will be used to batch correct as well
# sampleCols <- rep('grey',ncol(HNSCC.TCGA.df))
# names(sampleCols) <- colnames(HNSCC.TCGA.df)
HNSCCCols <- rep(mycols[1],ncol(HNSCC.RNAseq.its))
names(HNSCCCols) <- colnames(HNSCC.RNAseq.its)
# use HNSCC_labels vector created above to set HNSCC cols
HNSCCCols[rownames(HNSCC_labels)[HNSCC_labels == "T"]] <- mycols[12]
HNSCCCols[rownames(HNSCC_labels)[HNSCC_labels == "N"]] <- mycols[2]

#plot(1,1,col=mycols[12],cex=10) # test color
# TCGA cols
# ERROR CHECK - first make sure the TCGA samples in meta data are in the same order as the rows in my expression df
stopifnot(metaTCGA$barcode == colnames(expTCGA.its))
TCGACols <- mycols[3:35][as.factor(metaTCGA$project_id)]

# heatmap of enriched TFs
# subset genes to those in the TF network I'm working with
expanded_genes <- names(TF_gene_list)
# 
subdf <- df[df$gene_symbols %in% expanded_genes,]
expanded_genes.sig <- subdf[subdf$GENE_padj <= 0.05,]$gene_symbols


# heatmap 
pdf("TF_heatmap.HNSCC.TFs.pdf")
heatmap.2(HNSCC.RNAseq.its[rownames(HNSCC.RNAseq.its) %in% expanded_genes.sig,],
          scale='row',trace='none',
          dendrogram = 'col',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = HNSCCCols,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='kendall'))/2),
          main="Expanded TF Network")
dev.off()

# heatmap of expanded TFs list
# subset genes to those in the TF network I'm working with
expanded_genes <- names(expanded_TF_gene_list)

# heatmap 
pdf("TF_heatmap.HNSCC.expanded_TFs.pdf")
heatmap.2(HNSCC.RNAseq.its[rownames(HNSCC.RNAseq.its) %in% expanded_genes,],
          scale='row',trace='none',
          dendrogram = 'col',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = HNSCCCols,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='kendall'))/2),
          main="Expanded TF Network")
dev.off()

# heatmap of enriched TFs network (i.e. including target genes)
# subset genes to those in the TF network I'm working with
expanded_genes <- dplyr::bind_rows(TF_gene_list, .id = "column_label")
expanded_genes <- unique(c(as.character(expanded_genes$column_label),
                           as.character(expanded_genes$TF_genes)))

# heatmap 
pdf("TF_heatmap.HNSCC.TF_network.pdf")
heatmap.2(HNSCC.RNAseq.its[rownames(HNSCC.RNAseq.its) %in% expanded_genes,],
          scale='row',trace='none',
          dendrogram = 'col',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = HNSCCCols,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='kendall'))/2),
          main="Expanded TF Network")
dev.off()

# heatmap of enriched TFs network (i.e. including target genes)
# subset genes to those in the TF network I'm working with
expanded_genes <- dplyr::bind_rows(expanded_TF_gene_list, .id = "column_label")
expanded_genes <- unique(c(as.character(expanded_genes$column_label),
                           as.character(expanded_genes$TF_genes)))

# heatmap 
pdf("TF_heatmap.HNSCC.expanded_TF_network.pdf")
heatmap.2(HNSCC.RNAseq.its[rownames(HNSCC.RNAseq.its) %in% expanded_genes,],
          scale='row',trace='none',
          dendrogram = 'col',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = HNSCCCols,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='kendall'))/2),
          main="Expanded TF Network")
dev.off()

#################################################################
# now plot heatmaps using TCGA samples
#################################################################

# heatmap of enriched TFs
# subset genes to those in the TF network I'm working with
expanded_genes <- names(TF_gene_list)

# heatmap 
pdf("TF_heatmap.TCGA.TFs.pdf")
par(mar=c(1,1,1,1))
heatmap.2(expTCGA.its[rownames(expTCGA.its) %in% expanded_genes.sig,],
          scale='row',
          trace='none', key = F,
          dendrogram = 'none',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = TCGACols,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='pearson'))/2),
          main="Expanded TF Network")

legend('left',
       legend = unique(metaTCGA$project_id),
       col = unique(TCGACols), 
       lty= 1,             
       lwd = 5,          
       cex=.7
)
dev.off()

# heatmap 
# keep <- c('TCGA-ACC','TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC',
#           'TCGA-COAD', 'TCGA-ESCA', 'TCGA-HNSC', 'TCGA-KIRC',
#           'TCGA-KIRP', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC',
#           'TCGA-OV', 'TCGA-PAAD', 'TCGA-PRAD', 'TCGA-READ',
#           'TCGA-STAD', 'TCGA-THCA', 'TCGA-UCEC')

keep <- c('TCGA-BLCA', 'TCGA-CESC',
          'TCGA-ESCA', 'TCGA-HNSC')


keep_samples <- metaTCGA$project_id %in% keep &
  metaTCGA$sample_type == 'Primary Tumor'

# subset genes to those in the TF network I'm working with
expanded_genes <- dplyr::bind_rows(expanded_TF_gene_list, .id = "column_label")
expanded_genes <- unique(c(as.character(expanded_genes$column_label),
                           as.character(expanded_genes$TF_genes)))
pdf("TF_heatmap.TCGAsquamous.TFs.pdf")
par(mar=c(1,1,1,1))
heatmap.2(expTCGA.its[rownames(expTCGA.its) %in% expanded_genes,keep_samples],
          scale='row',trace='none', key = F,
          dendrogram = 'none',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = TCGACols[keep_samples],
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='pearson'))/2),
          main="Expanded TF Network")

legend('left',
       legend = unique(metaTCGA$project_id[keep_samples]),
       col = unique(TCGACols[keep_samples]), 
       lty= 1,             
       lwd = 5,          
       cex=.7
)
dev.off()

# subset genes to those in the TF network I'm working with
expanded_genes <- names(expanded_TF_gene_list)

pdf("TF_heatmap.TCGAsquamous.expanded_TFs.pdf")
par(mar=c(1,1,1,1))
heatmap.2(expTCGA.its[rownames(expTCGA.its) %in% expanded_genes,keep_samples],
          scale='row',trace='none', key = F,
          dendrogram = 'none',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = TCGACols[keep_samples],
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='kendall'))/2),
          main="Expanded TF Network")

legend('left',
       legend = unique(metaTCGA$project_id[keep_samples]),
       col = unique(TCGACols[keep_samples]), 
       lty= 1,             
       lwd = 5,          
       cex=.7
)
dev.off()

# now plot just HNSCC samples from our data set and from TCGA
# first, batch correct using combat

# run combat to remove batch effects
hnscc_samples <- which(metaTCGA$project_id == "TCGA-HNSC")
batchdf <- cbind(HNSCC.RNAseq.its, expTCGA.its[,hnscc_samples])

# first store labels used to denote batches
cohort_labels <- substr(colnames(batchdf),
       1,
       4)

batchCols <- rep("black", ncol(batchdf))
batchCols[which(colnames(batchdf) %in% rownames(HNSCC_labels))] <- HNSCC_labels[colnames(batchdf)[which(colnames(batchdf) %in% rownames(HNSCC_labels))],]
batchCols[which(batchCols == 'T')] <- 'orange'
batchCols[which(batchCols == 'N')] <- 'pink'
hnscc_barcodes <- metaTCGA$barcode[metaTCGA$project_id == "TCGA-HNSC"]
batchCols[which(colnames(batchdf) %in% hnscc_barcodes)] <- 'yellow'
hnscc_normals <- metaTCGA$barcode[metaTCGA$project_id == "TCGA-HNSC" & metaTCGA$sample_type == 'Solid Tissue Normal']
batchCols[which(colnames(batchdf) %in% hnscc_normals)] <- 'grey'

# plot clusters before batch correction
pdf('TCGA_hnscc.combat.clusters.pdf')
plotColoredClusters(standard.pearson(batchdf),
                    labs=cohort_labels,
                    cols=batchCols)
dev.off()

# now run combat batch correction
batchdf.combat <- sva::ComBat(batchdf, cohort_labels)
# now replot clusters post batch correction
pdf('TCGA_hnscc.combat.clusters.pdf')
plotColoredClusters(standard.pearson(plotdf.combat),
                    labs=cohort_labels,
                    cols=batchCols)
dev.off()

# plot heatmap using batch correcte data
# subset genes to those in the TF network I'm working with
expanded_genes <- names(expanded_TF_gene_list)

pdf("TF_heatmap.HNSCC_TCGA.expanded_TFs.pdf")
par(mar=c(1,1,1,1))
heatmap.2(batchdf.combat[rownames(batchdf.combat) %in% expanded_genes,],
          scale='row',trace='none', key = F,
          dendrogram = 'none',
          col=greenred,
          labCol = "",
          labRow = "",
          ColSideColors = batchCols,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1 - cor(t(x),method='kendall'))/2),
          main="Expanded TF Network")

legend('left',
       legend = c('DG-N', 'DG-T', 'DG-O','TCGA-T', 'TCGA-N'),
       col = unique(batchCols), 
       lty= 1,             
       lwd = 5,          
       cex=.7
)
dev.off()
