# Identify  current working directory
getwd()
# Download required data into desired working directory, then set the wokring directory to this location
setwd()

# Install and load packages

install.packages('BiocManager')
library(BiocManager)
install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap','tidyverse'))
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
install.packages('RColorBrewer')
library(RColorBrewer)
install.packages('ggrepel')
library(ggrepel)


# Import the count data
sample_table = read.csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')

# Amalgamate data 
files = pull(sample_table, Run)
files = paste0('counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)

# Create and assign a vector for tximport 
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi = tximport(files, 
                 type='salmon',
                 tx2gene=gene_map,
                 ignoreTxVersion=TRUE) # converts transcript level to gene level abundance


# this didnt work originally, click extract zip on files


# Run the DESeq2 command for esimtation of size factors and dispersions and data normalisation
dds = DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)

# Create dispersion estimates plot  
plotDispEsts(dds)

# Transform data to the log2 scale and reduce impact of outliers and create data frame of PCA data 
rld = rlog(dds)
pca_data <- plotPCA(rld, intgroup='Group', returnData = TRUE)
show(pca_plot)

# Create PCA plot, using colour-blind friendly colours, adding titles and changing font sizes and themes of the plot 
pca_1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA plot demonstrating differences between three treatment conditions for alveolar macrophages") + 
  xlab("PC1: 60% variance") +
  ylab("PC2: 18% variance") +
  scale_color_manual(values = c("Naive" = "orange", "Allo24h" = "cornflowerblue", "Allo2h" = "purple")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    text = element_text(size = 12),  # 
    axis.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14) 
  )

# Add ellipses around the clustered samples 
pca_1 <- pca_1 + 
  stat_ellipse(geom = "polygon", aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  scale_fill_manual(values = alpha(c("Naive" = "orange", "Allo24h" = "cornflowerblue", "Allo2h" = "purple"), 0.2

# Generate colour schemes 
sda_colours <- brewer.pal(9,'RdPu')

# Create sample-sample distance matrix based on euclidean distances 
sample_distance = dist(t(assay(rld)), method='euclidian')
sample_distance_matrix = as.matrix(sample_distance)

# Create heatmap using this data, separating the data and borders for visualization, adding a title and matching panel colours to the PCA plot for corresponding samples.
heatmap_annotation = data.frame(group=colData(dds)[,c('Group')], row.names=rownames(colData(dds)))
pheatmap(sample_distance_matrix,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         annotation_col = heatmap_annotation,
         cutree_cols = 3,
         cutree_rows = 3,
         color = sda_colours,
         main = "Heatmap of sample distances for GSE116583", 
         annotation_colors = annotation_colours, 
         border_color = NA)

# Create a scree plot, with threshold line of 1 and sequential colouring
pca_result <- prcomp(sample_distance_matrix, scale. = TRUE)
eigenvalues <- pca_result$sdev^2
install.packages("viridisLite")
library(viridisLite
line_colours <- inferno(length(eigenvalues))
plot(1:length(eigenvalues), eigenvalues, type = "b", pch = 16,
     xlab = "Principal Component", ylab = "Eigenvalue", 
     main = "Scree Plot",
     col = line_colours)
abline(h = 1, col = "red", lty = 2)
total_variance <- sum(eigenvalues)
eigenvalue_percentage <- eigenvalues / total_variance * 100
text(1:length(eigenvalues), eigenvalues, labels = round(eigenvalues, 2), pos = 1, col = "black", cex = 0.8)


# Volcano plot (Allo24h vs naive) 
# DESEq2 differential expression analysis Allo24h vs naive  
results_table = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
summary(results_table)
results_tibble = as_tibble(results_table, rownames='ensembl_gene_id')
filtered_results = filter(results_tibble, complete.cases(results_tibble))
# Mutate the data to include logPVal 
filtered_results = mutate(filtered_results, logPVal = -log10(padj))
# View the data 
view(filtered_results)
# Add the significance column based  on whether data is significantly expressed (TRUEFALSE)
filtered_results = mutate(filtered_results,
                          significant = padj<0.05)

# Extract the top 10 differentially expressed genes and creat a vector required by the volcano plot 
annot_results = left_join(filtered_results, annotation)
annot_results = arrange(annot_results, padj)
View(head(annot_results, 10))
degs = filter(annot_results, abs(log2FoldChange) > 1 & padj < 0.05)
topdegs <- head(degs, 10)

#Create parameters required of the volcano plot  (thresholds) 
filtered_results$combined_significance1 <- ifelse(
  filtered_results$significant == FALSE,
  'non-significant',
  ifelse(filtered_results$log2FoldChange > 1, 'upregulated',
         ifelse(filtered_results$log2FoldChange < -1, 'downregulated', 'p value only'))

#Create volcano plot with colour-blind friendly colours, with  threshold lines, titles and labels of the top 10 differentially expressed genes in boxes 
olcano_1 <- ggplot(filtered_results, aes(x = log2FoldChange, y = logPVal)) +
  geom_point(aes(colour = combined_significance1)) +
  scale_colour_manual(
    values = c(
      'upregulated' = 'lightcoral',
      'non-significant' = 'grey',
      'downregulated' = 'mediumorchid',
      'p value only' = 'blue'),
    name = 'Expression Change') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'black') +
  labs(title = ' Differential Expression Allo24h vs Naive') +
  geom_vline(xintercept = -1, linetype = 2, color = 'black') +
  geom_vline(xintercept = 1, linetype = 2, color = 'black')

volcano_genes <- volcano_1 + geom_label_repel(data = topdegs,
                                              aes(x = log2FoldChange, y = logPVal, label = external_gene_name),
                                              box.padding = 0.5, force = 1, fontface = 'bold', size = 3) +
  geom_segment(data = topdegs,
               aes(x = log2FoldChange, y = logPVal,
                   xend = log2FoldChange, yend = logPVal),
               color = 'black', size = 1) +
  theme_bw() +
  theme(
    text = element_text(size = 12, face = 'bold'),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 10),  # Adjust the size as needed
    legend.title = element_text(size = 12),  # Adjust the size as needed
    legend.text = element_text(size = 10)  # Adjust the si
  )


# Gene annotation 
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                'start_position', 'end_position',
                                'strand', 'gene_biotype', 'external_gene_name',
                                'description'),
                   filters = 'ensembl_gene_id', values = filtered_results$ensembl_gene_id,
                   mart = ensembl108)

# Volcano plot (Allo2h vs naive)
# DESEq2 differential expression analysis Allo2h vs naive
results_table2 = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
summary(results_table2)
results_tibble2 = as_tibble(results_table2, rownames='ensembl_gene_id')
filtered_results2 = filter(results_tibble2, complete.cases(results_tibble2))
# Mutate the data to include logPVal
filtered_results2 = mutate(filtered_results2, logPVal = -log10(padj))
# View the data
view(filtered_results2)
# Add the significance column based  on whether data is significantly expressed (TRUEFALSE)
filtered_results2 = mutate(filtered_results2,
                          significant = padj<0.05)

# Extract the top 10 differentially expressed genes and creat a vector required by the volcano plot
annot_results2 = left_join(filtered_results2, annotation)
annot_results2 = arrange(annot_results, padj)
View(head(annot_results2, 10))
degs2 = filter(annot_results2, abs(log2FoldChange) > 1 & padj < 0.05)
topdegs2 <- head(degs2, 10)

#Create parameters required of the volcano plot  (thresholds)
filtered_results2$combined_significance2 <- ifelse(
  filtered_results2$significant == FALSE,
  'non-significant',
  ifelse(filtered_results2$log2FoldChange > 1, 'upregulated',
         ifelse(filtered_results2$log2FoldChange < -1, 'downregulated', 'p value only'))

#Create volcano plot with colour-blind friendly colours, with  threshold lines, titles and labels of the top 10 differentially exp>olcano_1 <- ggplot(filtered_results, aes(x = log2FoldChange, y = logPVal)) +
volcano_2 <- ggplot(filtered_results2, aes(x = log2FoldChange, y = logPVal)) +
  geom_point(aes(colour = combined_significance2)) +
  scale_colour_manual(
    values = c(
      'upregulated' = 'lightcoral',
      'non-significant' = 'grey',
      'downregulated' = 'mediumorchid',
      'p value only' = 'blue'),
    name = 'Expression Change') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'black') +
  labs(title = 'Differential Expression Allo2h vs Naive') +
  geom_vline(xintercept = -1, linetype = 2, color = 'black') +
  geom_vline(xintercept = 1, linetype = 2, color = 'black')

volcano_2genes <- volcano_2 + geom_label_repel(data = topdegs2,
                                              aes(x = log2FoldChange, y = logPVal, label = external_gene_name),
                                              box.padding = 0.5, force = 1, fontface = 'bold', size = 3) +
  geom_segment(data = topdegs2,
               aes(x = log2FoldChange, y = logPVal,
                   xend = log2FoldChange, yend = logPVal),
               color = 'black', size = 1) +
  theme_bw() +
  theme(
    text = element_text(size = 12, face = 'bold'),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.text = element_text(size = 5),
    element_line(size = 1, colour = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

# Gene annotation 
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                'start_position', 'end_position',
                                'strand', 'gene_biotype', 'external_gene_name',
                                'description'),
                   filters = 'ensembl_gene_id', values = filtered_results2$ensembl_gene_id,
                   mart = ensembl108)
