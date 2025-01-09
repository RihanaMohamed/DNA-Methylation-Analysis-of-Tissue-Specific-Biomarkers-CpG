#=============================== Group 4 R Final Project ===================================
#============== Identifying Significant CpG sites in DNA from Different Body Fluids ===================================
#===== Group 4: John Daniel Esguerra, Shahed Al Asmi, Rihana Mohamed, Rachel Shadoff ===================================

#Set working directory
getwd()

# ========= Load Necessary Libraries =========
# BiocManager::install("Biobase")
# BiocManager::install("GEOquery")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("Rtsne")
# install.packages("umap")
# install.packages("sva")
# BiocManager::install("minfi")

library(Biobase)
library(GEOquery)
library(dplyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(Rtsne)
library(umap)
library(EnhancedVolcano)
library(clusterProfiler)
library(pheatmap)
library(limma)
library(minfi)



#=============================== LOADING DATASET ===================================
gset <- getGEO("GSE55734", GSEMatrix = TRUE, getGPL = TRUE)
# Found 1 file(s)
# GSE55734_series_matrix.txt.gz

gset_forensic = gset[["GSE55734_series_matrix.txt.gz"]]



#==================================== Step 2: Dataset Exploration =====================================================

# Assay Data
assay_data = gset_forensic@assayData[["exprs"]] # Columns are sample IDs, rows are methylation sites 
dim(assay_data) # 485,577 CpG sites and 16 samples
head(assay_data) # Proportion of methylated CpG sites stored in assay_data


# Phenotype Data
pheno_data = pData(gset_forensic) # Columns are tissue and sample type, rows are sample IDs
head(pheno_data) 
dim(pheno_data) # 16 samples and 34 columns
colnames(pheno_data)
head(pheno_data$source_name_ch1) #ch1 contains sample tissue type
sampletypes <- table(pheno_data$source_name_ch1) # 6 blood, 4 saliva, 6 vaginal secretion


# Feature Data 
feature_data = fData(gset_forensic) # Additional information about each sample
dim(feature_data) # 485577   x  37
head(feature_data)


#==================================== Step 3: Data Cleaning =====================================================

# 1. Check for NA values ------------------------------------------------------

summary(assay_data) # Dataset contains NAs
# Adding a column to count NA values
assay_data_NAs <- mutate(as.data.frame(assay_data), 
                         NA_count = rowSums(is.na(assay_data)))

assay_data_NAs_filtered <- assay_data_NAs %>% filter(NA_count != 0)

# Change y axis labels to non-scientific notation - to remove commas (5,000 to 5000)
marks_no_sci <- function(x) format(x, big.mark = ",", scientific = FALSE)


# Barplot showing NA distribution
na_distrib <- ggplot(assay_data_NAs_filtered, aes(x = NA_count)) + 
  geom_bar(fill="purple") +
  labs(title = "Distribution of NA Values", x = "NA Count", y = "Number of CpG Sites") +
  scale_y_continuous(labels = marks_no_sci) +
  scale_x_discrete(name ="NA Count", 
                   breaks = c("1","2","3","4","5","6","7","8","9","10") , 
                   limits = c("1","2","3","4","5","6","7","8","9","10")) +
  geom_text(aes(label = after_stat(count)), stat = "count", 
            vjust = -0.2, 
            colour = "black")


# Removing NA Data
assay_data_NAremoved = na.omit(assay_data)
dim(assay_data_NAremoved) # 9543 incomplete records removed



# 2. Checking for duplicates ------------------------------------------------------------

# Check for duplicates in rows
# nrow(assay_data_NAremoved)  == count(unique(as.data.frame(assay_data_NAremoved)))# returns TRUE
nrow(assay_data_NAremoved)  == nrow(unique(as.data.frame(assay_data_NAremoved)))
unique_rows <- unique(as.data.frame(assay_data_NAremoved))
nrow(assay_data_NAremoved) == nrow(unique_rows)

# Check for duplicates in columns
colnames(assay_data_NAremoved) == (unique(colnames(as.data.frame(assay_data_NAremoved)))) # returns TRUE




# 3. Checking standard deviation --------------------------------------------------------------

# Check distribution of standard deviation at each site
row_stdev <- apply(assay_data_NAremoved, 1, sd)
assay_data_sd <- mutate(as.data.frame(assay_data_NAremoved), 
                        std_dev = apply(as.data.frame(assay_data_NAremoved), 1, sd))


# Barplot showing standard deviation distribution
stdev_plot <- ggplot(assay_data_sd, aes(x = std_dev)) + 
  geom_histogram(fill="purple") +
  labs(title = "Distribution of Standard Deviation", 
       x = "Standard Deviation", 
       y = "Number of CpG Sites") +
  scale_y_continuous(labels = marks_no_sci)


# ISOLATING BY THRESHOLDS 
# Top 25% most variable sites
top25_sd = assay_data_sd %>% filter(std_dev > quantile(std_dev, 0.75))
dim(top25_sd) # 119,009 records

# Top 10% most variable sites
top10_sd = assay_data_sd %>% filter(std_dev > quantile(std_dev, 0.90))
dim(top10_sd) # 47,604 records

# Top 5% most variable sites
top5_sd = assay_data_sd %>% filter(std_dev > quantile(std_dev, 0.95))
dim(top5_sd) # 23,802 records

# Top 1% most variable sites
top1_sd = assay_data_sd %>% filter(std_dev > quantile(std_dev, 0.99))
dim(top1_sd) # 



assay_data_cleaned_25 = top25_sd[,-17]  # removes col 17 (std dev)
# 
assay_data_cleaned_10 = top10_sd[,-17]

assay_data_cleaned = top5_sd[,-17] # removes col 17 (std dev)

# assay_data_cleaned_1 = top1_sd[,-17]




#==================================== Step 4: Dimensionality Reduction =============================================
# Clustering methods will show whether or not tissue type is a significant factor that separates samples. 
# If clusters are visible, we can then move on to looking at which specific sites are associated with each tissue type

# PCA                
pca <- prcomp(t(assay_data_cleaned)) # omits values of NA
pca_df <-data.frame(pca$x) %>%
  mutate(Tissue = pheno_data$`sample type:ch1`)
pca_df <- pca_df[, c(1:2, ncol(pca_df))] %>% 
  dplyr::rename(Dim1 = PC1, Dim2 = PC2) # Standardizing column names for use in create plot function


# pcaPlot<-data.frame(pca$x) 
# names(pca)
# head(pca$x)[,1:5]
# pca_plot <- ggplot(data = data.frame(pca$x),aes(x = PC1, y = PC2, col = pheno_data$source_name_ch1)) +
#   geom_point() +
#   labs(title = "PCA Top 5%") 

# pca_25 <-prcomp(t(assay_data_cleaned_25)) # omits values of NA
# pcaPlot<-data.frame(pca_25$x) 
# names(pca_25)
# head(pca_25$x)[,1:5]
# pca_25_plot <- ggplot(data = data.frame(pca_25$x),aes(x = PC1, y = PC2, col = pheno_data$source_name_ch1)) +
#   geom_point() + 
#   labs(title = "PCA Top 25%")
# 
# pca_10 <-prcomp(t(assay_data_cleaned_10))
# pcaPlot<-data.frame(pca_10$x) 
# names(pca_10)
# head(pca_10$x)[,1:5]
# pca_10_plot <-ggplot(data = data.frame(pca_10$x),aes(x = PC1, y = PC2, col = pheno_data$source_name_ch1)) +
#   geom_point() + 
#   labs(title = "PCA Top 10%") 


# pca_plot+ theme(legend.position = "none")|pca_10_plot+ theme(legend.position = "none")|pca_25_plot #comparison of pca plots per stdev



#tSNE
set.seed(123) # set seed for reproducibility
# Computing perplexity
num_samples <- ncol(assay_data_cleaned)
perplexity_value <- min(5, (num_samples - 1) / 3) # perplexity of 5

tsne_result <- Rtsne(t(assay_data_cleaned), perplexity = perplexity_value, verbose = TRUE)
tsne_df <- data.frame(Dim1 = tsne_result$Y[, 1], Dim2 = tsne_result$Y[, 2], Tissue = pheno_data$source_name_ch1)



#UMAP
umap_result <- umap(t(assay_data_cleaned))
umap_df <- data.frame(Dim1 = umap_result$layout[, 1], Dim2 = umap_result$layout[, 2], Tissue = pheno_data$source_name_ch1)


# Plotting PCA, tSNE, and UMAP using create_plot function

create_plot <- function(data, dim1_col = "Dim1", dim2_col = "Dim2", tissue_col = "Tissue", title) {
  ggplot(data, aes_string(x = dim1_col, y = dim2_col, color = tissue_col)) +
    geom_point() +
    labs(title = paste(title, "Plot"), x = dim1_col, y = dim2_col)
}

PCAplot <- create_plot(pca_df, title = "PCA")
tSNEPlot <- create_plot(tsne_df, title = "t-SNE")
UMAPPlot <- create_plot(umap_df, title = "UMAP")

PCAplot | tSNEPlot | UMAPPlot #

# pca_plot | tsne_plot|umap_plot
# PCAplot | tsne_plot|umap_plot


#==================================== Step 5: MA and MDS plots =====================================================

# MA and MDS plots
densityPlot(as.matrix(assay_data_cleaned), main = "Density of Beta Values")

normalized_data <- normalizeBetweenArrays(assay_data_cleaned, method = "quantile") 
# quantile normalization aims to match the data distribution across all of your samples

par(mfrow = c(1, 1))
densityPlot(as.matrix(assay_data_cleaned), main = "Before Normalization")
densityPlot(as.matrix(normalized_data), main = "After Normalization")

limma::plotMA(assay_data_cleaned, ylim = c(-1, 1))
limma::plotMA(normalized_data, ylim = c(-1, 1)) 

Tissue = pheno_data$source_name_ch1
color_pal <- as.factor(pheno_data$source_name_ch1)
palette_colors <- scales::hue_pal()(length(unique(color_pal)))
tissue_col_map <- setNames(palette_colors, levels(color_pal))

plotMDS(assay_data_cleaned,
        col = tissue_col_map[color_pal],  # Assign colors by tissue type
        main = "MDS Plot Before Normalization")
legend("bottomright", 
       legend = levels(color_pal), 
       fill = palette_colors, 
       title = "Tissue Type")

plotMDS(normalized_data,
        col = tissue_col_map[color_pal],  # Assign colors by tissue type
        main = "MDS Plot After Normalization")
legend("bottomleft", 
       legend = levels(color_pal), 
       fill = palette_colors, 
       title = "Tissue Type")

foren_design <- model.matrix(~ 0 + pheno_data$source_name_ch1)
colnames(foren_design) <- gsub("pheno_data\\$source_name_ch1", "", colnames(foren_design))
colnames(foren_design) <- make.names(colnames(foren_design)) # removes the space in vaginal secretion


# we need m values
beta_to_m <- function(beta) log2(beta / (1 - beta))
m_values <- beta_to_m(normalized_data)

fit <- lmFit(m_values, foren_design) #Fits a linear model to the M-values, accounting for tissue type differences using the matrix created.

contrast_matrix <- makeContrasts(
  Blood_vs_Saliva = Blood - Saliva,
  Blood_vs_Vaginal = Blood - Vaginal.secretion,
  Saliva_vs_Vaginal = Saliva - Vaginal.secretion,
  levels = foren_design
)

fit2 <- contrasts.fit(fit, contrast_matrix) # Contrasts are applied to linear fit model
fit2 <- eBayes(fit2) # Moderates the t-statistics using empirical Bayes to improve stability

limma::plotMA(fit2, main = "MA Plot with Tissue Types")


#==================================== Step 6: Volcano Plots =====================================================

# Volcano plots and contrasts

# Blood vs Saliva
top_blood_vs_saliva <- topTable(fit2, coef = "Blood_vs_Saliva", adjust.method = "BH", number = Inf)
significant_blood_vs_saliva <- top_blood_vs_saliva[
  top_blood_vs_saliva$adj.P.Val < 0.05 & abs(top_blood_vs_saliva$logFC) > 1, ]

# Blood vs Vaginal
top_blood_vs_vaginal <- topTable(fit2, coef = "Blood_vs_Vaginal", adjust.method = "BH", number = Inf)
significant_blood_vs_vaginal <- top_blood_vs_vaginal[
  top_blood_vs_vaginal$adj.P.Val < 0.05 & abs(top_blood_vs_vaginal$logFC) > 1, ]

# Saliva vs Vaginal
top_saliva_vs_vaginal <- topTable(fit2, coef = "Saliva_vs_Vaginal", adjust.method = "BH", number = Inf)
significant_saliva_vs_vaginal <- top_saliva_vs_vaginal[
  top_saliva_vs_vaginal$adj.P.Val < 0.05 & abs(top_saliva_vs_vaginal$logFC) > 1, ]

# Combine Results for All Contrasts - add column for 1 of 3 contrast types
significant_blood_vs_saliva$Contrast <- "Blood_vs_Saliva"
significant_blood_vs_vaginal$Contrast <- "Blood_vs_Vaginal"
significant_saliva_vs_vaginal$Contrast <- "Saliva_vs_Vaginal"

all_significant <- rbind(
  significant_blood_vs_saliva,
  significant_blood_vs_vaginal,
  significant_saliva_vs_vaginal
)

# View combined results
head(all_significant)

# Save Results to CSV
write.csv(all_significant, "significant_features_all_contrasts.csv", row.names = TRUE)

# Add significance labels for all contrasts with Hypermethylated, Hypomethylated, and Not Significant
top_blood_vs_saliva$Significant <- with(
  top_blood_vs_saliva,
  ifelse(adj.P.Val < 0.05 & logFC > 1, "Hypermethylated",
         ifelse(adj.P.Val < 0.05 & logFC < -1, "Hypomethylated", "Not Significant"))
)

top_blood_vs_vaginal$Significant <- with(
  top_blood_vs_vaginal,
  ifelse(adj.P.Val < 0.05 & logFC > 1, "Hypermethylated",
         ifelse(adj.P.Val < 0.05 & logFC < -1, "Hypomethylated", "Not Significant"))
)

top_saliva_vs_vaginal$Significant <- with(
  top_saliva_vs_vaginal,
  ifelse(adj.P.Val < 0.05 & logFC > 1, "Hypermethylated",
         ifelse(adj.P.Val < 0.05 & logFC < -1, "Hypomethylated", "Not Significant"))
)

# Add Contrast column to identify comparisons
top_blood_vs_saliva$Contrast <- "Blood_vs_Saliva"
top_blood_vs_vaginal$Contrast <- "Blood_vs_Vaginal"
top_saliva_vs_vaginal$Contrast <- "Saliva_vs_Vaginal"

# Combine all contrasts into a list for plotting
contrast_list <- list(
  "Blood_vs_Saliva" = top_blood_vs_saliva,
  "Blood_vs_Vaginal" = top_blood_vs_vaginal,
  "Saliva_vs_Vaginal" = top_saliva_vs_vaginal
)

# Loop over contrasts and generate volcano plots
for (contrast_name in names(contrast_list)) {
  contrast_data <- contrast_list[[contrast_name]]
  
  # Debug: Check the count of each category
  print(paste("Counts for", contrast_name, ":"))
  print(table(contrast_data$Significant))
  
  # Create volcano plot
  plot <- ggplot(contrast_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(
      "Not Significant" = "grey",
      "Hypermethylated" = "brown4",
      "Hypomethylated" = "cyan3"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "Log Fold Change (logFC)",
      y = "-log10 Adjusted P-Value",
      color = "Methylation Status"
    ) +
    theme_minimal()
  
  print(plot)  # Display the plot
}

## Unified Volcano Plot Test ##
# Combine full results for all contrasts
combined_results <- rbind(
  top_blood_vs_saliva,
  top_blood_vs_vaginal,
  top_saliva_vs_vaginal
)

# Add a Contrast column for identification
combined_results$Contrast <- c(
  rep("Blood_vs_Saliva", nrow(top_blood_vs_saliva)),
  rep("Blood_vs_Vaginal", nrow(top_blood_vs_vaginal)),
  rep("Saliva_vs_Vaginal", nrow(top_saliva_vs_vaginal))
)

# Add significance label
combined_results$Significant <- with(
  combined_results,
  ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Significant", "Not Significant")
)

# Unified volcano plot
ggplot(combined_results, aes(x = logFC, y = -log10(adj.P.Val), color = Contrast, shape = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  scale_color_manual(values = c("Blood_vs_Saliva" = "red", 
                                "Blood_vs_Vaginal" = "blue", 
                                "Saliva_vs_Vaginal" = "black")) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +
  labs(
    title = "Unified Volcano Plot for All Contrasts",
    x = "Log Fold Change (logFC)",
    y = "-log10 Adjusted P-Value",
    color = "Comparison",
    shape = "Significance"
  ) +
  theme_minimal()

#==================================== Step 7: HEAT MAP and ISOLATION of cpg sites =======================================
# Heatmap

# Blood-specific CpG sites
blood_specific <- significant_blood_vs_saliva[
  rownames(significant_blood_vs_saliva) %in% rownames(significant_blood_vs_vaginal) &
    !rownames(significant_blood_vs_saliva) %in% rownames(significant_saliva_vs_vaginal), ]

# Saliva-specific CpG sites
saliva_specific <- significant_blood_vs_saliva[
  rownames(significant_blood_vs_saliva) %in% rownames(significant_saliva_vs_vaginal) &
    !rownames(significant_blood_vs_saliva) %in% rownames(significant_blood_vs_vaginal), ]

# Vaginal-specific CpG sites
vaginal_specific <- significant_blood_vs_vaginal[
  rownames(significant_blood_vs_vaginal) %in% rownames(significant_saliva_vs_vaginal) &
    !rownames(significant_blood_vs_vaginal) %in% rownames(significant_blood_vs_saliva), ]

# Combine tissue-specific CpG sites
tissue_specific_cpgs <- unique(c(rownames(blood_specific), rownames(saliva_specific), rownames(vaginal_specific)))

# Subset normalized data for significant CpG sites
heatmap_data <- normalized_data[tissue_specific_cpgs, ]
colnames(heatmap_data) <- pheno_data$title

# Create a simplified annotation dataframe
annotated_col <- pheno_data[, c("source_name_ch1"), drop = FALSE]
rownames(annotated_col) <- colnames(heatmap_data)


cpg_heatmap <- pheatmap(
  heatmap_data, 
  annotation_col = annotated_col,  # Add tissue type metadata
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  main = "Heatmap of Tissue-Specific CpG Sites",
)



# ======= Filtering for Tissue Specific Sites =======
# trying with k=3 because there are 3 visible clusters in the heatmap
cpg_clusters <- cbind(heatmap_data, cluster = cutree(cpg_heatmap$tree_row, k = 3))
cpg_df <- as.data.frame(cpg_clusters)
cpg_df_c1 <- cpg_df[cpg_df$cluster == 1, ] # hypomethylated in vaginal secretions + saliva
cpg_df_c2 <- cpg_df[cpg_df$cluster == 2, ] # hypomethylated in blood + saliva, undetermined in vaginal secretions (some samples show hypermethylation, some show no change)
cpg_df_c3 <- cpg_df[cpg_df$cluster == 3, ] # hypermethylated in blood + saliva


blood_specific_pval <- blood_specific[blood_specific$adj.P.Val < 0.05, ] 
blood_specific_ordered <- blood_specific_pval[order(blood_specific_pval$adj.P.Val),]
blood_most_sig <- blood_specific_ordered[1:10, ] %>% 
  mutate(Tissue = "Blood")

blood_specific_FC <- mutate(blood_specific, abs(blood_specific$logFC))
blood_specific_ordered_FC <- blood_specific_FC[order(blood_specific_pval$logFC, decreasing= TRUE),]


saliva_specific_pval <- saliva_specific[saliva_specific$adj.P.Val < 0.05, ]
saliva_specific_ordered <- saliva_specific_pval[order(saliva_specific_pval$adj.P.Val),]
saliva_most_sig <- saliva_specific_ordered[1:10, ] %>%
  mutate(Tissue = "Saliva")

saliva_specific_FC <- mutate(saliva_specific, abs(saliva_specific$logFC))
saliva_specific_ordered_FC <- saliva_specific_FC[order(saliva_specific_pval$logFC, decreasing= TRUE),]


vaginal_specific_pval <- vaginal_specific[vaginal_specific$adj.P.Val < 0.05, ]
vaginal_specific_ordered <- vaginal_specific_pval[order(vaginal_specific_pval$adj.P.Val),]
vaginal_most_sig <- vaginal_specific_ordered[1:10, ] %>% 
  mutate(Tissue = "Vaginal")

vaginal_specific_FC <- mutate(vaginal_specific, abs(vaginal_specific$logFC))
vaginal_specific_ordered_FC <- vaginal_specific_FC[order(vaginal_specific_pval$logFC,decreasing= TRUE),]

most_sig_sites <- rbind(vaginal_most_sig, blood_most_sig, saliva_most_sig) # Combined results with 10 CpG sites per tissue 


# ANSWER to RESEARCH QUESTION: "Which CpG sites are associated with the DNA derived from different body fluid-specific samples?”
cpg_tissue_isolated <- most_sig_sites[8] # isolates the cpg site and the associated tissue

## Further exploration of identified CpG Sites
most_sig_sites <- most_sig_sites %>%
  rownames_to_column(var="CpGSites")

sigsitefeaturedata <- merge(x = most_sig_sites, y = feature_data, by.x="CpGSites", by.y = "ID", all.x=TRUE)

sigsitefeaturedata_filtered <- subset(sigsitefeaturedata, select = c("CpGSites","Tissue", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq", "Infinium_Design_Type", "UCSC_RefGene_Name","UCSC_RefGene_Group","UCSC_CpG_Islands_Name","Relation_to_UCSC_CpG_Island","Regulatory_Feature_Name" , "Regulatory_Feature_Group"))

genetable <- data.frame(
  "CpG Site" = sigsitefeaturedata_filtered$CpGSites,
  "Associated Tissue" = sigsitefeaturedata_filtered$Tissue,
  "Associated Gene"= sigsitefeaturedata_filtered$UCSC_RefGene_Name)

write.csv(genetable, "genetable.csv", row.names = FALSE)


