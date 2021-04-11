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

# source functions
source("TFTargetJaccard.R")

# load HNSCC data
df.tn <- as.data.frame(readRDS("../data/DEG.SED.df.rds"))
# Load JQ1 data
df.047 <- readRDS("../data/JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
df.090 <- readRDS("../data/JQ1_vs_DMSO_DE_results_O90_JQ1.rds")

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
# Now add gene symbols to df
df$gene_symbols <- symbols
# remove any rows with missing gene symbols
df <- subset(df, !is.na(gene_symbols))

####################################################################
# Make custom gene lists representing TF target genes
####################################################################

# read in TF enrichment data
TF.df <- readRDS("../data/TF.enrichment_analysis.rds") # get TFs whose target genes you want to analyze for overlap with genes comprising leading edges of enriched pathways

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
# now expand this list by seeing whether any of the gene in the TF lists are TFs themselves that are not already in the TF list. If they are, add them to the TF list, along with the genes they target.
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
pal <- colorRampPalette(c("blue", "purple",
                          "white", "orange", "red"))
colors <- pal(nColor)
# bin the data but account for extreme
cuts <- apply(nodes[,3:6], 2, cut,
              c(-Inf,seq(-3, 3, 0.5), Inf), labels=1:14)
# set colors using these bins
nodes$log2FC_col <- colors[as.numeric(cuts[,1])]
nodes$log2FC_JQ1047_col <- colors[as.numeric(cuts[,2])]
nodes$log2FC_JQ1090_col <- colors[as.numeric(cuts[,3])]
nodes$log2FC_JQ1_col <- colors[as.numeric(cuts[,4])]
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

# calculate x and y coords for TFs using jaccard dist and mds
# jac_mat <- TFTargetJaccard(TF_gene_list, dist = T)
# jac_mat[is.nan(jac_mat)] <- 1 # deal with Nan
# # now perform mds
# mds <- cmdscale(d = jac_mat, eig = TRUE, k = 2) # k is the number of dim
# x <- mds$points[,1]
# y <- mds$points[,2]
# coords <- cbind(x, y)
# # make sure row order is correct
# coords <- coords[match(nodes$label, rownames(coords)),]
# # rename cols
# colnames(coords) <- c("V1", "V2")

pdf("../figures/TF_network.pdf", height = 15, width = 15)
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

pdf("../figures/TF_network.JQ1.pdf", height = 15, width = 15)
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

# make legend
pdf("../figures/TF_network.legend.pdf")
plot(rep(1, 16), 1:16, col = c("white", colors, "white"),
     pch = 19, cex = 3)
text(x = rep(1, 15), 1:15,
     labels = as.character(c(-Inf,seq(-3, 3, 0.5), Inf)),
     pos  = 4, offset = 2, cex = 1.5)
     #adj = c(0, -1), cex = 0.5)
dev.off()

# plot TF distance matrix
# image(1:ncol(jac_mat), 1:ncol(jac_mat),
#       jac_mat, axes = FALSE, xlab="", ylab="")
# axis(1, 1:ncol(jac_mat), colnames(jac_mat),
#                   cex.axis = 1, las=3)

# now make graphs with fruchterman reingold determined layout
# make edges df using dplyr bind_rows function
# Note: the warnings produced by the following command are fine - they just corce factors into character vectors as thigns are combined
edges <- dplyr::bind_rows(TF_gene_list, .id = "column_label")
# rename cols
names(edges) <- c("label", "pair")

# get rid of self connected nodes
edges <- edges[!edges$label == edges$pair,]

# add edge weights
edges$weights <- table(edges$label)[edges$label] * 0.01
# spread TFs out a bit from one another
edges$weights[edges$label %in% names(TF_gene_list) &
                edges$pair %in% names(TF_gene_list)] <-
  edges$weights[edges$label %in% names(TF_gene_list) &
                  edges$pair %in% names(TF_gene_list)] * 3

# node should be all TFs and Gene
nodes <- as.data.frame(unique(c(names(TF_gene_list),
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

# generate a color ramp over log2FC
# Create a color palette
nColor <- 14
# colors <- paletteer_c(palette = "viridis::viridis", n = nColor)
#colors <- brewer.pal(n = nColor, name = "Blues")
pal <- colorRampPalette(c("blue", "purple",
                          "white", "orange", "red"))
colors <- pal(nColor)
# bin the data but account for extreme
cuts <- apply(nodes[,3:6], 2, cut,
              c(-Inf,seq(-3, 3, 0.5), Inf), labels=1:14)
# set colors using these bins
nodes$log2FC_col <- colors[as.numeric(cuts[,1])]
nodes$log2FC_JQ1047_col <- colors[as.numeric(cuts[,2])]
nodes$log2FC_JQ1090_col <- colors[as.numeric(cuts[,3])]
nodes$log2FC_JQ1_col <- colors[as.numeric(cuts[,4])]
# now assign NA values to a color
nodes$log2FC_col[which(is.na(nodes$log2FC))] <- "darkgrey"
nodes$log2FC_JQ1047_col[which(is.na(nodes$log2FC_JQ1047))] <- "darkgrey"
nodes$log2FC_JQ1090_col[which(is.na(nodes$log2FC_JQ1090))] <- "darkgrey"
nodes$log2FC_JQ1_col[which(is.na(nodes$log2FC_JQ1))] <- "darkgrey"

# redefine edges and nodes

# first make the graph from the data I prepared
ig <- graph_from_data_frame(d = edges,
                            vertices = nodes, directed = FALSE)

#coords <- layout_with_kk(ig,
coords <- layout_with_fr(ig, weights = edges$weights,)# * 0.2#, dim = 3,
#coords <- layout_with_drl(ig) * 3#, options=list(edge.cut=0))

pdf("../figures/TF_network.fr.pdf", height = 18, width = 18)
plot(ig,
     layout = coords,
     rescale = T,
     vertex.label = ifelse(V(ig)$status == "TF", V(ig)$name, NA),
     vertex.label.cex = ifelse(V(ig)$status == "TF", 1.5, NA),
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_col,
     vertex.size = ifelse(V(ig)$status == "TF", 6, 1.2),
     edge.width = 0.5,
     edge.color = "grey",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)

dev.off()

# make another network colored by post JQ1 (mean across both cell lines)

pdf("../figures/TF_network.fr.JQ1.pdf", height = 18, width = 18)
plot(ig,
     layout = coords,
     rescale = T,
     vertex.label = ifelse(V(ig)$status == "TF", V(ig)$name, NA),
     vertex.label.cex = ifelse(V(ig)$status == "TF", 1.5, NA),
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_JQ1_col,
     vertex.size = ifelse(V(ig)$status == "TF", 6, 1.2),
     edge.width = 0.5,
     edge.color = "grey",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)

dev.off()

####################################################################
# Prepare nodes and edges for networks that include the
# target TFs of our enriched TFs
####################################################################

# first redefine only TF list
only_TFs_list <- lapply(expanded_TF_gene_list, function(x){
  hit <- x$TF_genes %in% names(expanded_TF_gene_list)
  return(data.frame(TF_genes = as.character(x$TF_genes[hit])))
}
)

# make edges df using dplyr bind_rows function
# Note: the warnings produced by the following command are fine - they just corce factors into character vectors as thigns are combined
edges <- dplyr::bind_rows(only_TFs_list, .id = "column_label")
# rename cols
names(edges) <- c("label", "pair")
#edges$weights <- log2(table(edges$label)[edges$label] + 3) * 0.01

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
nodes$status[nodes$label %in% names(expanded_TF_gene_list)] <- "TF"

# add transcriptomic data as meta data
mtch <- match(nodes$label, df$gene_symbols)
nodes$log2FC <- df$GENE_log2FC[mtch]
nodes$log2FC_JQ1047 <- df$log2FoldChange_047[mtch]
nodes$log2FC_JQ1090 <- df$log2FoldChange_090[mtch]
# mean of cell line changes post JQ1
nodes$log2FC_JQ1 <- rowMeans(cbind(nodes$log2FC_JQ1047, nodes$log2FC_JQ1090), na.rm=TRUE)

# scale logFCs
# nodes$log2FC <- scale(df$GENE_log2FC[mtch], center = T, scale = T)
# nodes$log2FC_JQ1047 <- scale(df$log2FoldChange_047[mtch], center = T, scale = T)
# nodes$log2FC_JQ1090 <- scale(df$log2FoldChange_090[mtch], center = T, scale = T)

# generate a color ramp over log2FC
# Create a color palette
nColor <- 14
# colors <- paletteer_c(palette = "viridis::viridis", n = nColor)
#colors <- brewer.pal(n = nColor, name = "Blues")
pal <- colorRampPalette(c("blue", "purple",
                          "white", "orange", "red"))
colors <- pal(nColor)
# bin the data but account for extreme
cuts <- apply(nodes[,3:6], 2, cut,
              c(-Inf,seq(-3, 3, 0.5), Inf), labels=1:14)
# set colors using these bins
nodes$log2FC_col <- colors[as.numeric(cuts[,1])]
nodes$log2FC_JQ1047_col <- colors[as.numeric(cuts[,2])]
nodes$log2FC_JQ1090_col <- colors[as.numeric(cuts[,3])]
nodes$log2FC_JQ1_col <- colors[as.numeric(cuts[,4])]
# now assign NA values to a color
nodes$log2FC_col[which(is.na(nodes$log2FC))] <- "darkgrey"
nodes$log2FC_JQ1047_col[which(is.na(nodes$log2FC_JQ1047))] <- "darkgrey"
nodes$log2FC_JQ1090_col[which(is.na(nodes$log2FC_JQ1090))] <- "darkgrey"
nodes$log2FC_JQ1_col[which(is.na(nodes$log2FC_JQ1))] <- "darkgrey"

#############################################################
# Make network plots with only TFs
#############################################################

# plot with iRanges package

# first make the graph from the data I prepared
ig <- graph_from_data_frame(d = edges,
                            vertices = nodes, directed = FALSE)
#coords <- layout_with_kk(ig,
#coords <- layout_with_fr(ig, weights = edges$weights,)# * 0.2#, dim = 3,
#coords <- layout_with_drl(ig) * 3#, options=list(edge.cut=0))

# calculate x and y coords for TFs using jaccard dist and mds
jac_mat <- TFTargetJaccard(expanded_TF_gene_list, dist = T)
jac_mat[is.nan(jac_mat)] <- 1 # deal with Nan
# now perform mds
mds <- cmdscale(d = jac_mat, eig = TRUE, k = 2) # k is the number of dim
x <- mds$points[,1]
y <- mds$points[,2]
coords <- cbind(x, y)
# make sure row order is correct
coords <- coords[match(nodes$label, rownames(coords)),]
# rename cols
colnames(coords) <- c("V1", "V2")

pdf("../figures/TF_network.only_TFs.pdf", height = 18, width = 18)
plot(ig,
     layout = coords,
     rescale = T,
     vertex.label = ifelse(V(ig)$name %in%
                        names(TF_gene_list),
                           V(ig)$name,
                           ifelse(V(ig)$log2FC <= -1.8 |
                                    V(ig)$log2FC >= 1.8,
                                    V(ig)$name,
                                    NA)),
     vertex.label.cex = 1.3,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_col,
     vertex.size = ifelse(V(ig)$name %in% names(TF_gene_list),
                          5, 2),
     edge.width = 0.3,
     edge.color = "grey",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)
dev.off()

# make another network colored by post JQ1 (mean across both cell lines)

pdf("../figures/TF_network.only_TFs.JQ1.pdf", height = 18, width = 18)
plot(ig,
     layout = coords,
     rescale = T,
     vertex.label = ifelse(V(ig)$name %in%
                             names(TF_gene_list),
                           V(ig)$name,
                           ifelse(V(ig)$log2FC_JQ1 <= -1.8 |
                                    V(ig)$log2FC_JQ1 >= 1.8,
                                  V(ig)$name,
                                  NA)),
     vertex.label.cex = 1.3,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_JQ1_col,
     vertex.size = ifelse(V(ig)$name %in% names(TF_gene_list),
                          5, 2),
     edge.width = 0.3,
     edge.color = "grey",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)
dev.off()

# plot TF distance matrix
# image(1:ncol(jac_mat), 1:ncol(jac_mat),
#       jac_mat, axes = FALSE, xlab="", ylab="")
# axis(1, 1:ncol(jac_mat), colnames(jac_mat),
#      cex.axis = 0.1, las=3)

# first redefine only TF list
only_TFs_list <- lapply(expanded_TF_gene_list, function(x){
  hit <- x$TF_genes %in% TF_in_trustt
  return(data.frame(TF_genes = as.character(x$TF_genes[hit])))
}
)


# now make graphs with fruchterman reingold determined layout
# make edges df using dplyr bind_rows function
# Note: the warnings produced by the following command are fine - they just corce factors into character vectors as thigns are combined
edges <- dplyr::bind_rows(only_TFs_list, .id = "column_label")
# rename cols
names(edges) <- c("label", "pair")

# get rid of self connected nodes
edges <- edges[!edges$label == edges$pair,]

# add edge weights
edges$weights <- table(edges$label)[edges$label] * 0.01
# spread TFs out a bit from one another
edges$weights[edges$label %in% names(TF_gene_list) &
                edges$pair %in% names(TF_gene_list)] <-
  edges$weights[edges$label %in% names(TF_gene_list) &
                  edges$pair %in% names(TF_gene_list)] * 3

# node should be all TFs and Gene
nodes <- as.data.frame(unique(c(names(expanded_TF_gene_list),
                                as.character(edges$label),
                                as.character(edges$pair))))
# name col
names(nodes) <- "label"

# now add meta data to nodes df
# first whether the node is a TF or a gene
nodes$status <- "Gene"
nodes$status[nodes$label %in% TF_in_trustt] <- "TF"

# add transcriptomic data as meta data
mtch <- match(nodes$label, df$gene_symbols)
nodes$log2FC <- df$GENE_log2FC[mtch]
nodes$log2FC_JQ1047 <- df$log2FoldChange_047[mtch]
nodes$log2FC_JQ1090 <- df$log2FoldChange_090[mtch]
# mean of cell line changes post JQ1
nodes$log2FC_JQ1 <- rowMeans(cbind(nodes$log2FC_JQ1047, nodes$log2FC_JQ1090), na.rm=TRUE)

# generate a color ramp over log2FC
# Create a color palette
nColor <- 14
# colors <- paletteer_c(palette = "viridis::viridis", n = nColor)
#colors <- brewer.pal(n = nColor, name = "Blues")
pal <- colorRampPalette(c("blue", "purple",
                          "white", "orange", "red"))
colors <- pal(nColor)
# bin the data but account for extreme
cuts <- apply(nodes[,3:6], 2, cut,
              c(-Inf,seq(-3, 3, 0.5), Inf), labels=1:14)
# set colors using these bins
nodes$log2FC_col <- colors[as.numeric(cuts[,1])]
nodes$log2FC_JQ1047_col <- colors[as.numeric(cuts[,2])]
nodes$log2FC_JQ1090_col <- colors[as.numeric(cuts[,3])]
nodes$log2FC_JQ1_col <- colors[as.numeric(cuts[,4])]
# now assign NA values to a color
nodes$log2FC_col[which(is.na(nodes$log2FC))] <- "darkgrey"
nodes$log2FC_JQ1047_col[which(is.na(nodes$log2FC_JQ1047))] <- "darkgrey"
nodes$log2FC_JQ1090_col[which(is.na(nodes$log2FC_JQ1090))] <- "darkgrey"
nodes$log2FC_JQ1_col[which(is.na(nodes$log2FC_JQ1))] <- "darkgrey"

# redefine edges and nodes

# first make the graph from the data I prepared
ig <- graph_from_data_frame(d = edges,
                            vertices = nodes, directed = FALSE)

#coords <- layout_with_kk(ig,
coords <- layout_with_fr(ig, weights = edges$weights,)# * 0.2#, dim = 3,
#coords <- layout_with_drl(ig) * 3#, options=list(edge.cut=0))

pdf("../figures/TF_network.only_TFs.fr.pdf", height = 18, width = 18)
plot(ig,
     layout = coords,
     rescale = T,
     vertex.label = ifelse(V(ig)$name %in%
                             names(TF_gene_list),
                           V(ig)$name,
                           ifelse(V(ig)$log2FC <= -1.8 |
                                    V(ig)$log2FC >= 1.8,
                                  V(ig)$name,
                                  NA)),
     vertex.label.cex = 1.5,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_col,
     vertex.size = ifelse(V(ig)$name %in% names(TF_gene_list),
                        3.5, 1.5),
     edge.width = 0.3,
     edge.color = "grey",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)

dev.off()

# make another network colored by post JQ1 (mean across both cell lines)

pdf("../figures/TF_network.only_TFs.fr.JQ1.pdf", height = 18, width = 18)
plot(ig,
     layout = coords,
     rescale = T,
     vertex.label = ifelse(V(ig)$name %in%
                             names(TF_gene_list),
                           V(ig)$name,
                           ifelse(V(ig)$log2FC_JQ1 <= -1.8 |
                                    V(ig)$log2FC_JQ1 >= 1.8,
                                  V(ig)$name,
                                  NA)),
     vertex.label.cex = 1.5,
     vertex.label.color = "black",
     vertex.label.dist = 0,
     vertex.color = V(ig)$log2FC_JQ1_col,
     vertex.size = ifelse(V(ig)$name %in% names(TF_gene_list),
                          3.5, 1.5),
     edge.width = 0.3,
     edge.color = "grey",
     vertex.frame.color = NA,#plasma(2, alpha = 0.1)[as.factor(V(ig)$sample_type)],
)


dev.off()

####################################################################
##### Run customGSEA analysis (see customGSEA script and url below)
####################################################################

# Here you will run GSEA on two different gene lists:
# 1) The TFs we see dSED enrichment for and their gene targets
# 2) The expanded TF + target + TF-target-targets
# The GSEA will look to see if these genes are enriched for 
# certain pathways.
# To do this, I must first prune the GO Terms pathways so that
# they only include genes from the list I'm analyzing
# (i.e. 1 of 2 above). Otherwise, this would likely bias the
# analysis.

# first make ranked gene list - this will be subset below
# for each analysis
ranked_gene_list = df$GENE_log2FC / df$GENE_lfcSE # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_symbols
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
#myPaths <- fgsea::gmtPathways(c5bp_file)
myPaths <- fgsea::gmtPathways(c6onc_file)

# FIRST, DO ANALYSIS WITH ONLY TF + TARGET GENES
# remove gene from pathways that are not in the master list
# for this GSEA analysis
TF_target_genes <- unique(c(as.character(dplyr::bind_rows(TF_gene_list)$TF_genes), as.character(names(TF_gene_list)))) # should return warnings - they reflect object coercion (factor to character) and are ok to ignore
myPaths.filt <- lapply(myPaths, function(x){
    hits <- x %in% TF_target_genes
    return(x[hits])
  })
# get rid of pathways with no genes left - if there are any
if (sum(sapply(myPaths.filt, function(x){
  length(x) == 0
})) != 0){
  myPaths.filt <- myPaths.filt[-which(sapply(myPaths.filt, function(x){
      length(x) == 0
    }))]
}
# subset ranked gene list to only the genes in analysis
ranked_gene_list.TF_target_genes <-
  ranked_gene_list[TF_target_genes]
# set genes names
names(ranked_gene_list.TF_target_genes) <- TF_target_genes
# sort this list again
ranked_gene_list.TF_target_genes <-
  sort(ranked_gene_list.TF_target_genes, decreasing = T)
# make names lower case and prettier
# names(myPaths.filt) <- sapply(names(myPaths.filt),
#   function(x){
#     gsub( "go_", "", tolower(x),)
#   })

# run fgsea
fgseaRes <- fgsea(pathways = myPaths.filt, 
                  stats = ranked_gene_list.TF_target_genes,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)


# collapse pathways
# collapsedPathways <- collapsePathways(fgseaRes[order(fgseaRes$pval),],
#                                       myPaths.filt,
#                                       ranked_gene_list.TF_target_genes)
# 
# mainPathways <- fgseaRes[fgseaRes$pathway %in% collapsedPathways$mainPathways,]
# upPaths <- mainPathways[mainPathways$NES > 0,]
# upPaths <- upPaths[head(order(upPaths$NES, decreasing = T),10),]
# downPaths <- mainPathways[mainPathways$NES < 0,]
# downPaths <- downPaths[head(order(downPaths$NES, decreasing = F),10),]
#topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

sigPaths <- subset(fgseaRes, padj <= 0.1)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
#downPaths <- sigPaths[sigPaths$NES < 0,]
#downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths)#, downPaths[nrow(downPaths):1,])

pdf("GSEA.TF_target_genes.c6.pdf", height=6, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths.filt[topPaths$pathway],
              ranked_gene_list.TF_target_genes, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# now repeat this GSEA analysis for TF + TF TARGET GENES +
# TARGET TF TARGET GENES
TF_target_genes <- unique(c(as.character(dplyr::bind_rows(expanded_TF_gene_list)$TF_genes), as.character(names(expanded_TF_gene_list)))) # should return warnings - they reflect object coercion (factor to character) and are ok to ignore
myPaths.filt <- lapply(myPaths, function(x){
  hits <- x %in% TF_target_genes
  return(x[hits])
})

# get rid of pathways with no genes left - if there are any
if (sum(sapply(myPaths.filt, function(x){
  length(x) == 0
})) != 0){
  myPaths.filt <- myPaths.filt[-which(sapply(myPaths.filt, function(x){
    length(x) == 0
  }))]
}

# subset ranked gene list to only the genes in analysis
ranked_gene_list.TF_target_genes <-
  ranked_gene_list[TF_target_genes]
# set genes names
names(ranked_gene_list.TF_target_genes) <- TF_target_genes
# sort this list again
ranked_gene_list.TF_target_genes <-
  sort(ranked_gene_list.TF_target_genes, decreasing = T)

# run fgsea
fgseaRes <- fgsea(pathways = myPaths.filt, 
                  stats = ranked_gene_list.TF_target_genes,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)

sigPaths <- subset(fgseaRes, padj <= 0.1)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

pdf("GSEA.TF_target_genes.expanded.c6.pdf", height=6, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths.filt[topPaths$pathway],
              ranked_gene_list.TF_target_genes, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# JQ1 treatment

# calculate means across both cell lines treated with JQ1
df$log2FoldChange_JQ1 <- rowMeans(cbind(
  df$log2FoldChange_047, df$log2FoldChange_090), na.rm=TRUE)
df$lfcSE_JQ1 <- rowMeans(cbind(
  df$lfcSE_047, df$lfcSE_090), na.rm=TRUE)

#set gene against which to filter pathway genes
JQ1_target_genes <- df$gene_symbols
# filter pathway genes
myPaths.filt <- lapply(myPaths, function(x){
  hits <- x %in% JQ1_target_genes
  return(x[hits])
})
# error check - no pathways should have zero genes after filtering 
stopifnot(sum(sapply(myPaths.filt, length) == 0) == 0)

# make new ranked gene list using JQ1 gene expression levels
ranked_gene_list = df$log2FoldChange_JQ1 / df$lfcSE_JQ1 # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_symbols
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# run fgsea
fgseaRes <- fgsea(pathways = myPaths.filt, 
                  stats = ranked_gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)

sigPaths <- subset(fgseaRes, padj <= 0.1)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

# output figure showing positive pathways enrichments
pdf("GSEA.JQ1.c6.pdf", height=13, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths.filt[topPaths$pathway],
              ranked_gene_list, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# now repeat for hallmark pathways
myPaths <- fgsea::gmtPathways(h_all_file)

# FIRST, DO ANALYSIS WITH ONLY TF + TARGET GENES
# remove gene from pathways that are not in the master list
# for this GSEA analysis
TF_target_genes <- unique(c(as.character(dplyr::bind_rows(TF_gene_list)$TF_genes), as.character(names(TF_gene_list)))) # should return warnings - they reflect object coercion (factor to character) and are ok to ignore
myPaths.filt <- lapply(myPaths, function(x){
  hits <- x %in% TF_target_genes
  return(x[hits])
})
# get rid of pathways with no genes left - if there are any
if (sum(sapply(myPaths.filt, function(x){
  length(x) == 0
})) != 0){
  myPaths.filt <- myPaths.filt[-which(sapply(myPaths.filt, function(x){
    length(x) == 0
  }))]
}

# subset ranked gene list to only the genes in analysis
ranked_gene_list.TF_target_genes <-
  ranked_gene_list[TF_target_genes]
# set genes names
names(ranked_gene_list.TF_target_genes) <- TF_target_genes
# sort this list again
ranked_gene_list.TF_target_genes <-
  sort(ranked_gene_list.TF_target_genes, decreasing = T)

# run fgsea
fgseaRes <- fgsea(pathways = myPaths.filt, 
                  stats = ranked_gene_list.TF_target_genes,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)

sigPaths <- subset(fgseaRes, padj <= 0.1)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
#downPaths <- sigPaths[sigPaths$NES < 0,]
#downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths)#, downPaths[nrow(downPaths):1,])

pdf("GSEA.TF_target_genes.h-all.pdf", height=6, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths.filt[topPaths$pathway],
              ranked_gene_list.TF_target_genes, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# now repeat this GSEA analysis for TF + TF TARGET GENES +
# TARGET TF TARGET GENES
TF_target_genes <- unique(c(as.character(dplyr::bind_rows(expanded_TF_gene_list)$TF_genes), as.character(names(expanded_TF_gene_list)))) # should return warnings - they reflect object coercion (factor to character) and are ok to ignore
myPaths.filt <- lapply(myPaths, function(x){
  hits <- x %in% TF_target_genes
  return(x[hits])
})
# get rid of pathways with no genes left - if there are any
if (sum(sapply(myPaths.filt, function(x){
  length(x) == 0
})) != 0){
  myPaths.filt <- myPaths.filt[-which(sapply(myPaths.filt, function(x){
    length(x) == 0
  }))]
}

# subset ranked gene list to only the genes in analysis
ranked_gene_list.TF_target_genes <-
  ranked_gene_list[TF_target_genes]
# set genes names
names(ranked_gene_list.TF_target_genes) <- TF_target_genes
# sort this list again
ranked_gene_list.TF_target_genes <-
  sort(ranked_gene_list.TF_target_genes, decreasing = T)

# run fgsea
fgseaRes <- fgsea(pathways = myPaths.filt, 
                  stats = ranked_gene_list.TF_target_genes,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)

sigPaths <- subset(fgseaRes, padj <= 0.1)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

pdf("GSEA.TF_target_genes.expanded.h-all.pdf", height=6, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths.filt[topPaths$pathway],
              ranked_gene_list.TF_target_genes, fgseaRes, 
              gseaParam = 0.5)
dev.off()

# now repeat this GSEA analysis for all genes after
# JQ1 treatment

# calculate means across both cell lines treated with JQ1
df$log2FoldChange_JQ1 <- rowMeans(cbind(
  df$log2FoldChange_047, df$log2FoldChange_090), na.rm=TRUE)
df$lfcSE_JQ1 <- rowMeans(cbind(
  df$lfcSE_047, df$lfcSE_090), na.rm=TRUE)

#set gene against which to filter pathway genes
JQ1_target_genes <- df$gene_symbols
# filter pathway genes
myPaths.filt <- lapply(myPaths, function(x){
  hits <- x %in% JQ1_target_genes
  return(x[hits])
})
# error check - no pathways should have zero genes after filtering 
stopifnot(sum(sapply(myPaths.filt, length) == 0) == 0)

# make new ranked gene list using JQ1 gene expression levels
ranked_gene_list = df$log2FoldChange_JQ1 / df$lfcSE_JQ1 # use logfoldchange normalized by standard eror, which should recreate the gene level stat
# set names to gene symbols
names(ranked_gene_list) = df$gene_symbols
# sort descending
ranked_gene_list = sort(ranked_gene_list, decreasing = TRUE)
# remove duplicate gene symbols
ranked_gene_list = ranked_gene_list[!duplicated(names(ranked_gene_list))]

# run fgsea
fgseaRes <- fgsea(pathways = myPaths.filt, 
                  stats = ranked_gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)

sigPaths <- subset(fgseaRes, padj <= 0.1)
upPaths <- sigPaths[sigPaths$NES > 0,]
upPaths <- upPaths[order(upPaths$NES, decreasing = T),]
# there are no significant down paths
downPaths <- sigPaths[sigPaths$NES < 0,]
downPaths <- downPaths[order(downPaths$NES, decreasing = F),]
topPaths <- rbind(upPaths, downPaths[nrow(downPaths):1,])

# output figure showing positive pathways enrichments
pdf("GSEA.JQ1.h-all.pdf", height=13, width=18)
# now repeat with just normal specific enhancers
plotGseaTable(myPaths.filt[topPaths$pathway],
              ranked_gene_list, fgseaRes, 
              gseaParam = 0.5)
dev.off()

