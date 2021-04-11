# replot Enhancer Figure - Log Fold Change X Distance X Tissue Specificity X Num Close Genes

######################################################################
# Set up working environment
######################################################################

# packages
library(ggplot2)
library(hexbin)
library(ggridges)
library(viridis)
library(ggpubr)
library(reshape)
library(plot3D)
library(org.Hs.eg.db)

#pval_cutoffs <- c(1,0.05,0.01,0.001)

#for (pval_cutoff in pval_cutoffs){
df <- as.data.frame(readRDS("../data/DEG.SED.df.rds"))

# annotate gene symbols
symbols <- mapIds(org.Hs.eg.db, rownames(df), 'SYMBOL', 'ENSEMBL')
# Now add gene symbols to df
df$gene_symbols <- symbols
df$absdistance <- abs(df$distance)
df$SE_type <- "UPPP"
df$SE_type[df$SE_logFC > 0] <- 'Tumor'
df$distbin <- cut(log10(df$absdistance), breaks = c(0, 3, 4 ,5 ,6, 7, 8), right = F, labels = c(3, 4 ,5 ,6, 7, 8))
           
# subset based on gene and se pvals
#subdf <- subset(df, SE_PValue <= 0.05 & GENE_padj <= 0.05 & absdistance < 10000000)
subdf <- subset(df, SE_PValue <= 0.05 & absdistance < 10000000)

subdf$SE_rank <- rank(-subdf$SE_logFC)
subdf$GENE_rank <- rank(-subdf$GENE_log2FC)
cor(subdf$SE_logFC, subdf$GENE_log2FC, use="complete.obs")
cor(subdf$SE_rank, subdf$GENE_rank, use="complete.obs")
ggplot(subdf,
       aes(x = SE_logFC, y = GENE_log2FC)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~distbin) +
  geom_smooth(method='lm', formula = y~x) +
  stat_cor(label.x = -5, label.y = 5, cex = 4)

p1 <- ggplot(subdf, aes(y = distbin)) +   # print figure to out device when in for loop in R
  geom_density_ridges(
    aes(x = GENE_log2FC, fill = SE_type), 
    alpha = .7, color = "white"
  ) +
  #scale_fill_manual(values = c("green", "purple")) +
  scale_fill_manual(values = viridis(2)) +
  labs(
    x = "Fold Change (log2)",
    y = "Distance (MB)",
    title = "Tumor vs Normal Gene Expression Differences By Distance To Nearest SE"
  ) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme(plot.margin = unit(c(1,0,1,1), "lines")) + 
  theme_classic()

p1



# plot pval
# pdf(paste0(outDir,"/pval_dist_SE.pdf"))
# ggplot(df.sub.filt, aes(y = distbin)) +
#   geom_density_ridges(
#     aes(x = padj, fill = SE_type), 
#     alpha = .7, color = "white"
#   ) +
#   scale_fill_manual(values = c("Red", "Blue")) +
#   labs(
#     x = "Fold Change (log2)",
#     y = "Distance (MB)",
#     title = "Tumor vs Normal Gene Expression Differences By Distance To Nearest SE"
#   ) +
#   scale_y_discrete(expand = c(0.01, 0)) +
#   scale_x_continuous(expand = c(0.01, 0)) 
# dev.off()


distbin <- c()
wilc_pval <- c()
shift_pval <- c()
# Use ks test to see the difference between log fold changes for each distance category
for (each in levels(df.sub.filt$distbin)){
  testdf <- subset(df.sub.filt, distbin == each)
  #ks <- ks.test(subset(testdf, SE_type == "Normal")$log2FoldChange,
  #              subset(testdf, SE_type == "Tumor")$log2FoldChange, alternative = "two.sided")
  # save test and pval metadata
  normal.dist <- subset(testdf, SE_type == "Normal")$log2FoldChange
  tumor.dist <- subset(testdf, SE_type == "Tumor")$log2FoldChange
  wilc <- wilcox.test(normal.dist, tumor.dist, alternative = "less")
  distbin <- append(distbin,each)
  #pval <- append(pval,ks$p.value)
  wilc_pval <- append(wilc_pval,wilc$p.value)
  
  # Now, see if the normal distribution is further left than the tumor one is to
  # the right at smaller distance bins
  shift.test <- wilcox.test(abs(normal.dist), abs(tumor.dist), alternative = "greater")
  shift_pval <- append(shift_pval, shift.test$p.value)
}

outdf <- as.data.frame(cbind(distbin, wilc_pval, shift_pval))
write.csv(outdf, file = paste0(outDir,"/lfc_dist_SE.ks_stats.pval",pval_cutoff,".csv"))
cols <- viridis(7)
# attempt to manualy make a line that is colored by a gradient representing wilcoxon sig tests/vals
pdf(paste0("pval_gradient",pval_cutoff,".pdf"), height = 10, width = 10)
outdf$zero = rep(0,nrow(outdf))
outdf$one = rep(1,nrow(outdf))
outdf$wilc_pval <- as.numeric(as.character(outdf$wilc_pval))
outdf$shift_pval <- as.numeric(as.character(outdf$shift_pval))
outdf$distbin <- as.numeric(as.character(outdf$distbin))
p2 <- ggplot() +
  coord_cartesian(xlim = c(-5, 10)) +
  geom_point(data = outdf, aes(x=zero, y = distbin, color = wilc_pval), size = 3) + 
  geom_point(data = outdf, aes(x=one, y = distbin, color = shift_pval), size = 3) +
  theme_void() +
  scale_color_gradient(#colors = rev(magma(7)), #low = low_col, mid = mid_col, high = high_col, n.breaks = 5,
    low = cols[6],
    high = cols[2],
    breaks = c(min(outdf$wilc_pval),0.1,0.3,1.0),
    labels = c(paste0(round(min(outdf$wilc_pval))),expression("0.1"),expression("0.3"),expression("1.0")),
    #midpoint = 0.01,
    limits = c(0,1)) +
  theme(plot.margin = unit(c(1,5,1,8), "lines")) +
  theme(legend.title = element_blank())
print(p2)
dev.off()
