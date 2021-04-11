# Using the raw counts from the Cistrome analysis (Ilya/Vanja)
# calculate odds ratios, and pvals, and prioritize TFs (especially ones enriched in tumor DSEs)
# Also, make plots of enrichment odds ratios for each TF
# comparing tumor vs UPP enrichment across promoters,
# typical enhancers (TE), and super enhancers (SEs)
# Note: I will calculate odds ratios two different ways - 
# To calculate odds ratios for tumors and UPPP independently,
# I will use fisher exact test on a 2x2 count matrix (i.e. a separate matrix for tumor and UPPP feature overlap)
# To determine enrichment significance in tumors compared with normals,
# I will use a fisher like framework (FLINT).

######################################################################
# Set up working environment
######################################################################

# packages
library(readxl)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(scales)
library(viridis)
# source my functions
source("ParseCistromeCounts.R")
source("ContingenciesResults.R")

######################################################################
# Read and pre-process data
######################################################################

# I'm using the older dataset to match the TF names Ilya/Vanja used to the typical gene names for TFs
mapping.df <- read.table("../data/sites-vs-nosites_mapping.tsv", sep = "\t", header = T)
# This is the newer dataset with the counts
counts.df <- read.table("../data/sites-vs-nosites.5_24_20.Ilya.tsv", sep = "\t", header = T) # this has correct odds ratios
# set cistrome TF name to character in dataset used for mapping
mapping.df$TF <- as.character(mapping.df$TF)
# add mapping to counts.df
mtch <- match(counts.df$TF, mapping.df$TF)
counts.df$Transcription.factor <- mapping.df$Transcription.factor[mtch]
# rename col names
rename_col_names <- c("TF_cistrome_names",                                                                
                      "UPPP.Promoter.odds.ratio",
                      "UPPP.Promoter.num_sites",
                      "UPPP.Promoter.num_not_sites",
                      "UPPP.Promoter.neg_log10_pval_corr",
                      "UPPP.Promoter.num_sites.10k",
                      "UPPP.Promoter.num_not_sites.10k",
                      "UPPP.Promoter.num_sites.100k",
                      "UPPP.Promoter.num_not_sites.100k",
                      "Tumor.Promoter.odds.ratio",
                      "Tumor.Promoter.num_sites",
                      "Tumor.Promoter.num_not_sites",
                      "Tumor.Promoter.neg_log10_pval_corr",
                      "Tumor.Promoter.num_sites.10k",
                      "Tumor.Promoter.num_not_sites.10k",
                      "Tumor.Promoter.num_sites.100k",
                      "Tumor.Promoter.num_not_sites.100k",
                      "UPPP.Enhancer.odds.ratio",
                      "UPPP.Enhancer.num_sites",
                      "UPPP.Enhancer.num_not_sites",
                      "UPPP.Enhancer.neg_log10_pval_corr",
                      "UPPP.Enhancer.num_sites.10k",
                      "UPPP.Enhancer.num_not_sites.10k",
                      "UPPP.Enhancer.num_sites.100k",
                      "UPPP.Enhancer.num_not_sites.100k",
                      "Tumor.Enhancer.odds.ratio",
                      "Tumor.Enhancer.num_sites",
                      "Tumor.Enhancer.num_not_sites",
                      "Tumor.Enhancer.neg_log10_pval_corr",
                      "Tumor.Enhancer.num_sites.10k",
                      "Tumor.Enhancer.num_not_sites.10k",
                      "Tumor.Enhancer.num_sites.100k",
                      "Tumor.Enhancer.num_not_sites.100k",
                      "UPPP.SE.odds.ratio",
                      "UPPP.SE.num_sites",
                      "UPPP.SE.num_not_sites",
                      "UPPP.SE.neg_log10_pval_corr",
                      "UPPP.SE.num_sites.10k",
                      "UPPP.SE.num_not_sites.10k",
                      "UPPP.SE.num_sites.100k",
                      "UPPP.SE.num_not_sites.100k",
                      "Tumor.SE.odds.ratio",
                      "Tumor.SE.num_sites",
                      "Tumor.SE.num_not_sites",
                      "Tumor.SE.neg_log10_pval_corr",
                      "Tumor.SE.num_sites.10k",
                      "Tumor.SE.num_not_sites.10k",
                      "Tumor.SE.num_sites.100k",
                      "Tumor.SE.num_not_sites.100k",
                      "Transcription.factor"
                      )
names(counts.df) <- rename_col_names

#######################################################
# Parse Counts
#######################################################
# use ParseTumorCounts function I wrote to parse tumor counts
# This will be used to calculate odds ratios of enrichment in tumors and normals (and to compare these odds ratios)

TF.cistrome.p.10k <- ParseCistromeCounts(counts.df, "Promoter", "10k", 5)
TF.cistrome.p.100k <- ParseCistromeCounts(counts.df, "Promoter", "100k", 5)
TF.cistrome.e.10k <- ParseCistromeCounts(counts.df,  "Enhancer", "10k", 5)
TF.cistrome.e.100k <- ParseCistromeCounts(counts.df, "Enhancer", "100k", 5)
TF.cistrome.se.10k <- ParseCistromeCounts(counts.df, "SE", "10k", 5)
TF.cistrome.se.100k <- ParseCistromeCounts(counts.df, "SE", "100k", 5)

#####################################################################
# Try my Z-score 2x2x2 approximation
#####################################################################

# calculate 2x2x2 comparison of 2x2 odds ratios
TF.cistrome.p.10k.res <- ContingenciesResults(TF.cistrome.p.10k)
TF.cistrome.p.100k.res <- ContingenciesResults(TF.cistrome.p.100k)
TF.cistrome.e.10k.res <- ContingenciesResults(TF.cistrome.e.10k)
TF.cistrome.e.100k.res <- ContingenciesResults(TF.cistrome.e.100k)
TF.cistrome.se.10k.res <- ContingenciesResults(TF.cistrome.se.10k)
TF.cistrome.se.100k.res <- ContingenciesResults(TF.cistrome.se.100k)

# merge the results from 10k and 100k approaches
TF.cistrome.p.res <- merge(TF.cistrome.p.10k.res, TF.cistrome.p.100k.res, by = "TF")
TF.cistrome.e.res <- merge(TF.cistrome.e.10k.res, TF.cistrome.e.100k.res, by = "TF")
TF.cistrome.se.res <- merge(TF.cistrome.se.10k.res, TF.cistrome.se.100k.res, by = "TF")

# rename cols
names(TF.cistrome.p.res) <- c("TF", "delta10k", "pval10k","tumorOR10k", "normalOR10k", "delta100k", "pval100k","tumorOR100k", "normalOR100k")
names(TF.cistrome.e.res) <- c("TF", "delta10k", "pval10k","tumorOR10k", "normalOR10k", "delta100k", "pval100k","tumorOR100k", "normalOR100k")
names(TF.cistrome.se.res) <- c("TF", "delta10k", "pval10k","tumorOR10k", "normalOR10k", "delta100k", "pval100k","tumorOR100k", "normalOR100k")

# filter out TFs with -+Inf as there delta10k/delta100k values
TF.cistrome.p.filt <- subset(TF.cistrome.p.res, !is.infinite(delta10k) &
                        !is.infinite(delta100k) &
                        !is.nan(delta10k) &
                        !is.nan(delta100k)
                        )

TF.cistrome.e.filt <- subset(TF.cistrome.e.res, !is.infinite(delta10k) &
                               !is.infinite(delta100k) &
                               !is.nan(delta10k) &
                               !is.nan(delta100k)
                             )

TF.cistrome.se.filt <- subset(TF.cistrome.se.res, !is.infinite(delta10k) &
                               !is.infinite(delta100k) &
                               !is.nan(delta10k) &
                               !is.nan(delta100k)
                              )

# add adjusted pvals
TF.cistrome.p.filt$padj10k <- p.adjust(TF.cistrome.p.filt$pval10k, method = "bonferroni")
TF.cistrome.p.filt$padj100k <- p.adjust(TF.cistrome.p.filt$pval100k, method = "bonferroni")
TF.cistrome.e.filt$padj10k <- p.adjust(TF.cistrome.e.filt$pval10k, method = "bonferroni")
TF.cistrome.e.filt$padj100k <- p.adjust(TF.cistrome.e.filt$pval100k, method = "bonferroni")
TF.cistrome.se.filt$padj10k <- p.adjust(TF.cistrome.se.filt$pval10k, method = "bonferroni")
TF.cistrome.se.filt$padj100k <- p.adjust(TF.cistrome.se.filt$pval100k, method = "bonferroni")

# add annotation
TF.cistrome.p.filt$analysis <- "Promoter"
TF.cistrome.e.filt$analysis <- "Typical Enhancer"
TF.cistrome.se.filt$analysis <- "Super Enhancer"

#####################################################################
# Make Plots
#####################################################################

# Plot tumor OR on yaxis and UPPP OR on xaxis.
# Do this for promoters, enhancers, and SEs.
# Also, on these plots, highlight the
# P/E/SE that are signficant after multiple testing
# correction in both 10k and 100k analyses.

plt.df <- rbind(TF.cistrome.p.filt,
                TF.cistrome.e.filt,
                TF.cistrome.se.filt
                )
# change factor level order for analysis variable
plt.df$analysis <- factor(plt.df$analysis,
                          levels = c("Promoter",
                                     "Typical Enhancer",
                                     "Super Enhancer"))

# average ORs and deltas between 10k and 100k approaches
# for pval - make a new col with the max pval from 10k and 100k approaches
plt.df$tumorOR <- pmin(plt.df$tumorOR10k,
                                 plt.df$tumorOR100k, na.rm=TRUE)
plt.df$normalOR <- pmin(plt.df$normalOR10k,
                                  plt.df$normalOR100k, na.rm=TRUE)
plt.df$delta <- pmin(plt.df$delta10k,
                               plt.df$delta100k, na.rm=TRUE)
plt.df$padj <- pmax(plt.df$padj10k, plt.df$padj100k)

# output this df as an R object in case you need it later
saveRDS(object = plt.df, file = paste0("../data/TF.enrichment_analysis.rds"))

# define signficant rows as padj <= 0.05
# note that this ensures that signficant pvals are
# significant in both 10k and 100k analyses (which is conservative)

# storing the signficant rows is done for convenience, as it
# will allow me too change the signficance threshold
# in one place!

sig.rows <- plt.df$padj <= 0.05 

# define cols
colblind_pal <- colorblind_pal()(8) # colorblind palette - commenting out
pltcols <- rep("black", nrow(plt.df))
pltcols[sig.rows] <- colblind_pal[7]

# plot promoters
p.rows <- plt.df$analysis == "Promoter"
plt.p <- ggplot(plt.df[p.rows,],
            aes(x = normalOR, y = tumorOR)) +
  geom_point(col =pltcols[p.rows]) +
  xlab("UPPP OR") +
  ylab("Tumor OR") +
  ggtitle("Promoters") +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, color="grey", 
              linetype="dashed", size=1) +
  xlim(0, 10) +
  ylim(0, 16)

# plot typical enhancers
e.rows <- plt.df$analysis == "Typical Enhancer"
plt.e <- ggplot(plt.df[e.rows,],
                aes(x = normalOR, y = tumorOR)) +
  geom_point(col = pltcols[e.rows]) +
  xlab("UPPP OR") +
  ylab("Tumor OR") +
  ggtitle("Typical Enhancers") +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, color="grey", 
              linetype="dashed", size=1) +
  xlim(0, 10) +
  ylim(0, 16)


# plot super enhancers
se.rows <- plt.df$analysis == "Super Enhancer"
plt.se <- ggplot(plt.df[se.rows,],
                aes(x = normalOR, y = tumorOR)) +
  geom_point(col = pltcols[se.rows]) +
  xlab("UPPP OR") +
  ylab("Tumor OR") +
  ggtitle("Super Enhancers") +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, color="grey", 
              linetype="dashed", size=1) +
  xlim(0, 10) +
  ylim(0, 16)


# plot the signficant TFs and their delta values (delta represents the log OR tumor - log OR normal)
# This will require some preprocessing
# first find the rows from the df with TFs you want to plot
delta.rows <- (plt.df$analysis == "Super Enhancer" &
  sig.rows & plt.df$delta > 0) | (plt.df$analysis == "Typical Enhancer" &
                 sig.rows & plt.df$delta > 0) |
  (plt.df$analysis == "Promoter" &
     sig.rows & plt.df$delta > 0)

# now subset the df to these rows only
plt.df.delta <- plt.df[which(delta.rows),]
# store unique vector of these TFs
delta.TFs <- unique(as.character(plt.df.delta$TF))

# now assign cols depending on whether these TFs are signficant in TE, SE, or both.
# also parse out the delta and adjp values for plotting
# uses a for loop
delta.cols <- rep("black", length(delta.TFs)) # init color vector as "black"
names(delta.cols) <- delta.TFs # assign TF as rownames
delta.vals <- rep(0, length(delta.TFs)) # init - will update vals in loop below
names(delta.vals) <- delta.TFs # assign TFs as rownames
delta.padj <- rep(0, length(delta.TFs)) # init - will update vals in loop below
names(delta.padj) <- delta.TFs # assign TFs as rownames
# run for loop to do color asssignments and delta value parsing

for (TF in delta.TFs){
  temp.df <- plt.df.delta[plt.df.delta$TF == TF,]
  temp.analysis <- as.character(temp.df$analysis)
  temp.delta <- temp.df$delta
  test1 <- sum(grepl("Promoter", temp.analysis))
  test2 <- sum(grepl("Typical Enhancer", temp.analysis))
  test3 <- sum(grepl("Super Enhancer", temp.analysis))
  if (test1 == 1 & test2 == 1 & test3 == 1){
    delta.cols[TF] <- "P_TE_SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else if (test1 == 1 & test2 == 1){
    delta.cols[TF] <- "P_TE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$padj
  } else if (test1 == 1 & test3 == 1){
    delta.cols[TF] <- "P_SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else if (test2 == 1 & test3 == 1){
    delta.cols[TF] <- "TE_SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else if (test1 == 1){
    delta.cols[TF] <- "P"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Promoter",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Promoter",]$padj
  } else if (test2 == 1){
    delta.cols[TF] <- "TE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$padj
  } else if (test3 == 1){
    delta.cols[TF] <- "SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else{
    stop("This should not be possible - color loop, delta plot")
  }
}

col_map <- c(P_TE_SE = "purple",
             P_TE = "gold",
             P_SE = "lightpink",
             TE_SE = "hotpink",
             P = "blue",
             TE = "lightgreen",
             SE = "cyan"
)

# replace plt.df.delta with df derived from vectors you just made
plt.df.delta <- data.frame(delta = delta.vals,
                           cmap = delta.cols,
                           padj = delta.padj)

plt.df.delta$TF <- rownames(plt.df.delta)
# order TFs be delta value
plt.df.delta$TF <- factor(plt.df.delta$TF,
                          levels = plt.df.delta$TF[order(plt.df.delta$delta, decreasing = T)])

# save this df in a new variable
plt.TFs1 <- cbind(plt.df.delta)

# repeat the processing above, but now for TFs significantly
# enriched in normal-specific features
# first find the rows from the df with TFs you want to plot
delta.rows <- (plt.df$analysis == "Super Enhancer" &
                 sig.rows & plt.df$delta < 0) |
  (plt.df$analysis == "Typical Enhancer" &
     sig.rows & plt.df$delta < 0) |
  (plt.df$analysis == "Promoter" &
     sig.rows & plt.df$delta < 0)

# now subset the df to these rows only
plt.df.delta <- plt.df[which(delta.rows),]
# store unique vector of these TFs
delta.TFs <- unique(as.character(plt.df.delta$TF))

# now assign cols depending on whether these TFs are signficant in TE, SE, or both.
# also parse out the delta and adjp values for plotting
# uses a for loop
delta.cols <- rep("black", length(delta.TFs)) # init color vector as "black"
names(delta.cols) <- delta.TFs # assign TF as rownames
delta.vals <- rep(0, length(delta.TFs)) # init - will update vals in loop below
names(delta.vals) <- delta.TFs # assign TFs as rownames
delta.padj <- rep(0, length(delta.TFs)) # init - will update vals in loop below
names(delta.padj) <- delta.TFs # assign TFs as rownames
# run for loop to do color asssignments and delta value parsing

for (TF in delta.TFs){
  temp.df <- plt.df.delta[plt.df.delta$TF == TF,]
  temp.analysis <- as.character(temp.df$analysis)
  temp.delta <- temp.df$delta
  test1 <- sum(grepl("Promoter", temp.analysis))
  test2 <- sum(grepl("Typical Enhancer", temp.analysis))
  test3 <- sum(grepl("Super Enhancer", temp.analysis))
  if (test1 == 1 & test2 == 1 & test3 == 1){
    delta.cols[TF] <- "P_TE_SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else if (test1 == 1 & test2 == 1){
    delta.cols[TF] <- "P_TE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$padj
  } else if (test1 == 1 & test3 == 1){
    delta.cols[TF] <- "P_SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else if (test2 == 1 & test3 == 1){
    delta.cols[TF] <- "TE_SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else if (test1 == 1){
    delta.cols[TF] <- "P"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Promoter",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Promoter",]$padj
  } else if (test2 == 1){
    delta.cols[TF] <- "TE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Typical Enhancer",]$padj
  } else if (test3 == 1){
    delta.cols[TF] <- "SE"
    delta.vals[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$delta
    delta.padj[TF] <- temp.df[temp.df$analysis == "Super Enhancer",]$padj
  } else{
    stop("This should not be possible - color loop, delta plot")
  }
}

# replace plt.df.delta with df dervied from vectors you just made
plt.df.delta <- data.frame(delta = delta.vals,
                           cmap = delta.cols,
                           padj = delta.padj)

plt.df.delta$TF <- rownames(plt.df.delta)
# order TFs be delta value
plt.df.delta$TF <- factor(plt.df.delta$TF,
                          levels = plt.df.delta$TF[order(plt.df.delta$delta, decreasing = T)])

# combine these TFs with previously determined tumor specific
# TFs

plt.df.delta <- rbind(plt.TFs1, plt.df.delta)
# sort df
plt.df.delta <- plt.df.delta[order(plt.df.delta$delta, decreasing = T),]

# now plot the figure
# Note, this represents plt.d (delta) and plt.d (panel d in the combined figure)
plt.d <- ggplot(plt.df.delta, # plt.dn == delta normal
                aes(x = TF,
                    y = delta,
                    color = cmap)) +
  #geom_point(col = delta.cols, size = -log10(plt.df.delta$padj)) + # note: commenting this out because I don't like how this looks with size scaled to pval - I don't think it adds much either as this is already filtered for signficantly enriched TFs
  geom_point(size = 3) +
  geom_segment(aes(xend = TF), col = "black", yend = 0) +
  labs(x = "TF",
       y = expression(paste(Delta, "log(OR)")),
       color = "Legend") +
  ggtitle("Significantly Enriched TFs") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  scale_color_manual(values = col_map[levels(plt.df.delta$cmap)]) +
  # scale_color_colorblind() + # looks terrible
  # scale_color_viridis(discrete = TRUE, option = "C") + # all options look terrible
  theme(legend.position="bottom") +
  ylim(-3,3) +
  geom_hline(yintercept = 0, linetype="solid", 
             color = "black", size=2)

# combined plots using ggpubr
# using two calls to ggarrange, one inside the other
# this should allow for the alignment I want
pdf("../figures/TFs_enrichments.pdf", width = 9, height = 6,
    onefile = F)
ggarrange(ggarrange(plt.p, plt.e, plt.se,
                    ncol = 3,
                    labels = c("A", "B", "C")),
          # ggarrange(plt.dt,
          #           plt.dn,
          #           ncol = 2, 
          #           labels = c("D", "E")),
          plt.d,
          nrow = 2,
          common.legend = T,
          legend = "bottom",
          labels = c("", "D")
)
dev.off()

#############################################################
# plot 10k vs 100k values for supplemental qc figure
# compare 10k vs 100k
plt.10kvs100k <- ggplot(plt.df,
              aes(x = delta10k,
                  y = delta100k,
                  color = analysis)
                  ) +
              geom_point(shape = 16) +
              theme_classic(base_size = 20) +
              geom_abline(intercept = 0, slope = 1, linetype = 2, col = "darkgrey") +
              scale_color_manual(values=alpha(c("cyan","hotpink", "lightgreen"), 1)) +
              labs(color = "Legend") +
              theme(legend.position = "bottom") +
              xlab(expression(paste(Delta, "10k"))) +
              ylab(expression(paste(Delta, "100k")))

#############################################################
# Evaluate normality assumption
#############################################################

# evaluate to see distributions of delta10k given normality
# assumption inherent to method I used above

# plot hists and qqs for delta10k and delta100k
# Note: save these in case they are useful as supplemental figures

# graph as density histogram
plt.hist.10k <-ggplot(plt.df,
                      aes(x = delta10k)) +
  geom_density() +
  theme_classic() + 
  xlab(expression(paste(Delta, "10k")))

# graph qq with line
plt.qq.10k <- ggplot(plt.df,
                      aes(sample = delta10k)) +
  geom_qq() +
  geom_qq_line() +
  theme_classic() +
  ylab(expression(paste(Delta, "10k")))

# graph as density histogram
plt.hist.100k <-ggplot(plt.df,
                      aes(x = delta100k)) +
  geom_density() +
  theme_classic() +
  xlab(expression(paste(Delta, "100k")))

# graph qq with line
plt.qq.100k <- ggplot(plt.df,
                     aes(sample = delta100k)) +
  geom_qq() +
  geom_qq_line() +
  theme_classic() +
  ylab(expression(paste(Delta, "100k")))

# output figures as a combined figure with 5 panels
pdf("../figures/TFs_enrichments.qc.pdf", onefile = F)
ggarrange(
  ggarrange(
    plt.hist.10k,
    plt.hist.100k,
    plt.qq.10k,
    plt.qq.100k,
    ncol = 2,
    nrow = 2,
    labels = c("A", "B", "C", "D"),
    common.legend = T
  ),
  plt.10kvs100k,
  ncol = 1,
  nrow = 2, 
  labels = c("", "E"),
  common.legend = T,
  legend = "bottom"
  
)
dev.off()

