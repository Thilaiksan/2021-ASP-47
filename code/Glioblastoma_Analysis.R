
############################################################
# PROJECT TITLE:
# Gene Expression-Based Biomarker Discovery and
# Machine Learning Analysis of Glioblastoma
#
# GEO DATASET: GSE4290
# ORGANISM: Homo sapiens
# PLATFORM: GPL570
#
# PROJECT OBJECTIVES:
# 1. Compare glioblastoma and non-tumour brain samples.
# 2. Identify differentially expressed genes.
# 3. Identify potential molecular biomarkers.
# 4. Perform GO and KEGG enrichment analysis.
# 5. Classify samples using SVM and Random Forest.
# 6. Evaluate models using 5-fold cross-validation.
#
# IMPORTANT:
# Run the complete script from beginning to end.
############################################################


############################################################
# 1. LOAD REQUIRED PACKAGES
############################################################

library(GEOquery)
library(Biobase)
library(limma)
library(ggplot2)
library(pheatmap)
library(hgu133plus2.db)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(caret)
library(e1071)
library(randomForest)
library(dplyr)


############################################################
# 2. SET PROJECT DIRECTORY
############################################################

project_path <- "C:/Users/hp/Desktop/Glioblastoma_Bioinformatics_Project"

setwd(project_path)

dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

cat("\nPROJECT DIRECTORY:", getwd(), "\n")


############################################################
# 3. IMPORT GEO DATASET
############################################################

if (!file.exists("data/GSE4290_GEO_Data.rds")) {
  stop("GSE4290_GEO_Data.rds not found in the data folder.")
}

gse_list <- readRDS("data/GSE4290_GEO_Data.rds")

gse <- gse_list[[1]]

expression_data <- exprs(gse)

metadata <- pData(gse)

cat("\nDATASET INFORMATION\n")

cat("Platform:", annotation(gse), "\n")
cat("Expression dimensions:", dim(expression_data), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")

if (annotation(gse) != "GPL570") {
  stop("Unexpected platform. Check the platform annotation.")
}


############################################################
# 4. SELECT GLIOBLASTOMA AND NON-TUMOUR SAMPLES
############################################################

diagnosis <- metadata$characteristics_ch1

gbm_samples <- diagnosis ==
  "Histopathological diagnostic: glioblastoma, grade 4"

control_samples <- diagnosis ==
  "Histopathological diagnostic: non-tumor"

selected_samples <- gbm_samples | control_samples

filtered_expression <- expression_data[
  ,
  selected_samples,
  drop = FALSE
]

filtered_metadata <- metadata[
  selected_samples,
  ,
  drop = FALSE
]

sample_group <- ifelse(
  gbm_samples[selected_samples],
  "Glioblastoma",
  "Non-tumour"
)

sample_group <- factor(
  sample_group,
  levels = c("Non-tumour", "Glioblastoma")
)

cat("\nSAMPLE GROUP COUNTS\n")
print(table(sample_group))

cat("\nFILTERED EXPRESSION DIMENSIONS\n")
print(dim(filtered_expression))

if (ncol(filtered_expression) != length(sample_group)) {
  stop("Sample group and expression matrix dimensions do not match.")
}


############################################################
# 5. QUALITY CONTROL
############################################################

cat("\nQUALITY CONTROL\n")

cat("Missing values:",
    sum(is.na(filtered_expression)), "\n")

cat("Infinite values:",
    sum(is.infinite(filtered_expression)), "\n")

cat("Expression range:\n")

print(range(filtered_expression, na.rm = TRUE))

if (any(!is.finite(filtered_expression))) {
  stop("Expression matrix contains missing or infinite values.")
}


############################################################
# 6. LOG2 TRANSFORMATION
############################################################

# GSE4290 expression values have a wide range.
# Apply log2(x + 1) to retain zero-valued measurements
# without generating missing values.

processed_expression <- log2(filtered_expression + 1)

cat("\nEXPRESSION RANGE AFTER LOG2 TRANSFORMATION\n")

print(range(processed_expression, na.rm = TRUE))


############################################################
# 7. QUANTILE NORMALIZATION
############################################################

normalized_expression <- normalizeBetweenArrays(
  processed_expression,
  method = "quantile"
)

cat("\nNORMALIZED EXPRESSION DIMENSIONS\n")

print(dim(normalized_expression))

cat("\nNORMALIZED EXPRESSION RANGE\n")

print(range(normalized_expression, na.rm = TRUE))

# Save normalized expression matrix

saveRDS(
  normalized_expression,
  "data/GSE4290_Normalized_Expression.rds"
)

# Normalized expression boxplot

png(
  "results/figures/Normalized_Expression_Boxplot.png",
  width = 1200,
  height = 800
)

boxplot(
  normalized_expression,
  outline = FALSE,
  las = 2,
  main = "Normalized Expression Distribution",
  ylab = "Log2 Expression Value"
)

dev.off()


############################################################
# 8. PRINCIPAL COMPONENT ANALYSIS
############################################################

pca_data <- t(normalized_expression)

gene_variance <- apply(
  pca_data,
  2,
  var,
  na.rm = TRUE
)

pca_data <- pca_data[
  ,
  is.finite(gene_variance) & gene_variance > 0,
  drop = FALSE
]

pca_result <- prcomp(
  pca_data,
  center = TRUE,
  scale. = TRUE
)

pca_plot_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = sample_group
)

pca_plot <- ggplot(
  pca_plot_data,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Glioblastoma Samples",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()

print(pca_plot)

ggsave(
  "results/figures/PCA_Glioblastoma.png",
  pca_plot,
  width = 8,
  height = 6
)


############################################################
# 9. DIFFERENTIAL EXPRESSION ANALYSIS USING LIMMA
############################################################

group <- sample_group

design <- model.matrix(~ 0 + group)

colnames(design) <- c(
  "Non_tumour",
  "Glioblastoma"
)

contrast_matrix <- makeContrasts(
  Glioblastoma - Non_tumour,
  levels = design
)

fit <- lmFit(
  normalized_expression,
  design
)

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(fit2)

deg_results <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

write.csv(
  deg_results,
  "results/all_DEG_results.csv"
)

cat("\nTOTAL PROBES ANALYSED:",
    nrow(deg_results), "\n")


############################################################
# 10. IDENTIFY SIGNIFICANT DEGs
############################################################

significant_degs <- deg_results[
  !is.na(deg_results$adj.P.Val) &
    deg_results$adj.P.Val < 0.05 &
    abs(deg_results$logFC) > 1,
  ,
  drop = FALSE
]

upregulated_genes <- significant_degs[
  significant_degs$logFC > 1,
  ,
  drop = FALSE
]

downregulated_genes <- significant_degs[
  significant_degs$logFC < -1,
  ,
  drop = FALSE
]

cat("\nDIFFERENTIALLY EXPRESSED GENES\n")

cat("Significant probes:",
    nrow(significant_degs), "\n")

cat("Upregulated probes:",
    nrow(upregulated_genes), "\n")

cat("Downregulated probes:",
    nrow(downregulated_genes), "\n")

write.csv(
  significant_degs,
  "results/significant_DEGs.csv"
)

write.csv(
  upregulated_genes,
  "results/Upregulated_Genes.csv"
)

write.csv(
  downregulated_genes,
  "results/Downregulated_Genes.csv"
)


############################################################
# 11. GENE ANNOTATION
############################################################

probe_ids <- rownames(deg_results)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

entrez_ids <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

deg_results$GeneSymbol <- as.character(
  gene_symbols[rownames(deg_results)]
)

deg_results$EntrezID <- as.character(
  entrez_ids[rownames(deg_results)]
)

annotated_degs <- deg_results[
  !is.na(deg_results$GeneSymbol) &
    deg_results$GeneSymbol != "",
  ,
  drop = FALSE
]

write.csv(
  annotated_degs,
  "results/annotated_DEG_results.csv"
)


############################################################
# 11A. FINAL SIGNIFICANT DEG SUMMARY AND BIOMARKER TABLES
############################################################

# Re-apply the project significance criteria after annotation.
# This keeps the annotated and unannotated probe counts explicit.

significant_check <- deg_results[
  !is.na(deg_results$adj.P.Val) &
    deg_results$adj.P.Val < 0.05 &
    abs(deg_results$logFC) > 1,
  ,
  drop = FALSE
]

upregulated_significant <- significant_check[
  significant_check$logFC > 1,
  ,
  drop = FALSE
]

downregulated_significant <- significant_check[
  significant_check$logFC < -1,
  ,
  drop = FALSE
]

annotated_significant <- significant_check[
  !is.na(significant_check$GeneSymbol) &
    significant_check$GeneSymbol != "",
  ,
  drop = FALSE
]

unannotated_significant_count <-
  nrow(significant_check) - nrow(annotated_significant)

unique_annotated_gene_count <- length(
  unique(annotated_significant$GeneSymbol)
)

# Save significant annotated probe rows separately.
write.csv(
  annotated_significant,
  "results/Significant_Annotated_DEGs.csv"
)

# Save the top 10 upregulated and downregulated annotated probes.
top10_upregulated <- annotated_significant[
  annotated_significant$logFC > 1,
  ,
  drop = FALSE
]
top10_upregulated <- top10_upregulated[
  order(-top10_upregulated$logFC),
  ,
  drop = FALSE
]
top10_upregulated <- head(top10_upregulated, 10)

write.csv(
  top10_upregulated,
  "results/Top10_Upregulated_Genes.csv"
)

top10_downregulated <- annotated_significant[
  annotated_significant$logFC < -1,
  ,
  drop = FALSE
]
top10_downregulated <- top10_downregulated[
  order(top10_downregulated$logFC),
  ,
  drop = FALSE
]
top10_downregulated <- head(top10_downregulated, 10)

write.csv(
  top10_downregulated,
  "results/Top10_Downregulated_Genes.csv"
)

# Create a reproducible summary using calculated values, not hard-coded counts.
deg_summary <- data.frame(
  Category = c(
    "Total significant probe sets",
    "Upregulated significant probe sets",
    "Downregulated significant probe sets",
    "Significant probes with gene symbols",
    "Significant probes without gene symbols",
    "Unique annotated gene symbols"
  ),
  Count = c(
    nrow(significant_check),
    nrow(upregulated_significant),
    nrow(downregulated_significant),
    nrow(annotated_significant),
    unannotated_significant_count,
    unique_annotated_gene_count
  )
)

cat("\nFINAL DEG SUMMARY\n")
print(deg_summary)

if (nrow(upregulated_significant) +
    nrow(downregulated_significant) != nrow(significant_check)) {
  stop("Upregulated and downregulated probe counts do not sum to all significant probes.")
}

write.csv(
  deg_summary,
  "results/DEG_Summary.csv",
  row.names = FALSE
)

cat("\nTOP 10 UPREGULATED ANNOTATED PROBES\n")
print(
  top10_upregulated[, c("GeneSymbol", "logFC", "adj.P.Val"), drop = FALSE]
)

cat("\nTOP 10 DOWNREGULATED ANNOTATED PROBES\n")
print(
  top10_downregulated[, c("GeneSymbol", "logFC", "adj.P.Val"), drop = FALSE]
)

cat("\nTOP 10 ANNOTATED GENES\n")

print(
  head(
    annotated_degs[
      order(annotated_degs$adj.P.Val),
      c("GeneSymbol", "logFC", "P.Value", "adj.P.Val")
    ],
    10
  )
)


############################################################
# 12. VOLCANO PLOT
############################################################

deg_results$Significance <- "Not significant"

deg_results$Significance[
  !is.na(deg_results$adj.P.Val) &
    deg_results$adj.P.Val < 0.05 &
    deg_results$logFC > 1
] <- "Upregulated"

deg_results$Significance[
  !is.na(deg_results$adj.P.Val) &
    deg_results$adj.P.Val < 0.05 &
    deg_results$logFC < -1
] <- "Downregulated"

volcano_plot <- ggplot(
  deg_results,
  aes(
    x = logFC,
    y = -log10(pmax(adj.P.Val, 1e-300)),
    color = Significance
  )
) +
  geom_point(alpha = 0.6) +
  labs(
    title = "Volcano Plot: Glioblastoma vs Non-tumour",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

print(volcano_plot)

ggsave(
  "results/figures/Volcano_Plot_GBM.png",
  volcano_plot,
  width = 8,
  height = 6
)


############################################################
# 13. TOP 20 DEG HEATMAP
############################################################

top20_degs <- head(
  significant_degs[
    order(significant_degs$adj.P.Val),
    ,
    drop = FALSE
  ],
  20
)

if (nrow(top20_degs) > 0) {
  
  top20_expression <- normalized_expression[
    rownames(top20_degs),
    ,
    drop = FALSE
  ]
  
  heatmap_data <- t(scale(t(top20_expression)))
  
  heatmap_data <- heatmap_data[
    complete.cases(heatmap_data),
    ,
    drop = FALSE
  ]
  
  sample_annotation <- data.frame(
    Group = sample_group
  )
  
  rownames(sample_annotation) <- colnames(heatmap_data)
  
  png(
    "results/figures/Top20_DEG_Heatmap.png",
    width = 1200,
    height = 900
  )
  
  pheatmap(
    heatmap_data,
    annotation_col = sample_annotation,
    show_colnames = FALSE,
    main = "Top 20 Differentially Expressed Probes"
  )
  
  dev.off()
}


############################################################
# 14. PREPARE GENE LISTS FOR ENRICHMENT
############################################################

# Convert all measured probes into Entrez Gene IDs.
# This provides the background universe for enrichment.

all_entrez <- mapIds(
  hgu133plus2.db,
  keys = rownames(normalized_expression),
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

gene_universe <- unique(
  na.omit(as.character(all_entrez))
)

# Upregulated Entrez IDs

up_entrez <- mapIds(
  hgu133plus2.db,
  keys = rownames(upregulated_genes),
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

up_gene_list <- unique(
  na.omit(as.character(up_entrez))
)

# Downregulated Entrez IDs

down_entrez <- mapIds(
  hgu133plus2.db,
  keys = rownames(downregulated_genes),
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

down_gene_list <- unique(
  na.omit(as.character(down_entrez))
)

cat("\nENRICHMENT GENE COUNTS\n")

cat("Background genes:",
    length(gene_universe), "\n")

cat("Upregulated genes:",
    length(up_gene_list), "\n")

cat("Downregulated genes:",
    length(down_gene_list), "\n")


############################################################
# 15. GO ENRICHMENT ANALYSIS
############################################################

go_up <- NULL
go_down <- NULL

if (length(up_gene_list) >= 5) {
  
  go_up <- enrichGO(
    gene = up_gene_list,
    universe = gene_universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
}

if (length(down_gene_list) >= 5) {
  
  go_down <- enrichGO(
    gene = down_gene_list,
    universe = gene_universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
}

go_up_df <- if (!is.null(go_up)) {
  as.data.frame(go_up)
} else {
  data.frame()
}

go_down_df <- if (!is.null(go_down)) {
  as.data.frame(go_down)
} else {
  data.frame()
}

cat("\nGO ENRICHMENT RESULTS\n")

cat("Upregulated GO terms:",
    nrow(go_up_df), "\n")

cat("Downregulated GO terms:",
    nrow(go_down_df), "\n")

write.csv(
  go_up_df,
  "results/GO_upregulated.csv",
  row.names = FALSE
)

write.csv(
  go_down_df,
  "results/GO_downregulated.csv",
  row.names = FALSE
)


############################################################
# 16. KEGG ENRICHMENT ANALYSIS
############################################################

kegg_up <- NULL
kegg_down <- NULL

if (length(up_gene_list) >= 5) {
  
  kegg_up <- tryCatch(
    enrichKEGG(
      gene = up_gene_list,
      universe = gene_universe,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    ),
    error = function(e) {
      message("KEGG upregulated error: ", e$message)
      NULL
    }
  )
  
}

if (length(down_gene_list) >= 5) {
  
  kegg_down <- tryCatch(
    enrichKEGG(
      gene = down_gene_list,
      universe = gene_universe,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    ),
    error = function(e) {
      message("KEGG downregulated error: ", e$message)
      NULL
    }
  )
  
}

kegg_up_df <- if (!is.null(kegg_up)) {
  as.data.frame(kegg_up)
} else {
  data.frame()
}

kegg_down_df <- if (!is.null(kegg_down)) {
  as.data.frame(kegg_down)
} else {
  data.frame()
}

cat("\nKEGG ENRICHMENT RESULTS\n")

cat("Upregulated KEGG pathways:",
    nrow(kegg_up_df), "\n")

cat("Downregulated KEGG pathways:",
    nrow(kegg_down_df), "\n")

write.csv(
  kegg_up_df,
  "results/KEGG_upregulated.csv",
  row.names = FALSE
)

write.csv(
  kegg_down_df,
  "results/KEGG_downregulated.csv",
  row.names = FALSE
)


############################################################
# 17. ENRICHMENT VISUALIZATION
############################################################

# clusterProfiler returns enrichResult objects.
# Create plots directly from the result data frames to keep
# the figure output consistent and reproducible.

create_enrichment_plot <- function(data, plot_title, file_name) {
  
  required_columns <- c(
    "Description",
    "RichFactor",
    "Count",
    "p.adjust"
  )
  
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required enrichment columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  if (nrow(data) == 0) {
    message("No enrichment terms available for: ", plot_title)
    return(NULL)
  }
  
  # Ensure plotting columns are numeric.
  data$RichFactor <- suppressWarnings(
    as.numeric(as.character(data$RichFactor))
  )
  
  data$Count <- suppressWarnings(
    as.numeric(as.character(data$Count))
  )
  
  data$p.adjust <- suppressWarnings(
    as.numeric(as.character(data$p.adjust))
  )
  
  # Remove rows with unusable plotting values.
  plot_data <- data[
    !is.na(data$Description) &
      is.finite(data$RichFactor) &
      is.finite(data$Count) &
      !is.na(data$p.adjust),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_data) == 0) {
    message("No valid rows available for: ", plot_title)
    return(NULL)
  }
  
  # Select the 15 most significant terms.
  plot_data <- plot_data[
    order(plot_data$p.adjust),
    ,
    drop = FALSE
  ]
  
  plot_data <- head(plot_data, 15)
  
  # Avoid infinite values when adjusted p-values are zero.
  plot_data$LogP <- -log10(
    pmax(plot_data$p.adjust, .Machine$double.xmin)
  )
  
  # Preserve the selected order in the plot.
  plot_data$Description <- factor(
    plot_data$Description,
    levels = rev(unique(plot_data$Description))
  )
  
  enrichment_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = RichFactor,
      y = Description,
      size = Count,
      color = LogP
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = plot_title,
      x = "Rich Factor",
      y = "Enriched Term",
      size = "Gene Count",
      color = "-log10 Adjusted P-value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(size = 8)
    )
  
  ggplot2::ggsave(
    filename = file_name,
    plot = enrichment_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
  
  return(enrichment_plot)
}

# Generate the four enrichment figures.

go_up_plot <- create_enrichment_plot(
  go_up_df,
  "GO Biological Processes - Upregulated",
  "results/figures/GO_Upregulated.png"
)

go_down_plot <- create_enrichment_plot(
  go_down_df,
  "GO Biological Processes - Downregulated",
  "results/figures/GO_Downregulated.png"
)

kegg_up_plot <- create_enrichment_plot(
  kegg_up_df,
  "KEGG Pathways - Upregulated",
  "results/figures/KEGG_Upregulated.png"
)

kegg_down_plot <- create_enrichment_plot(
  kegg_down_df,
  "KEGG Pathways - Downregulated",
  "results/figures/KEGG_Downregulated.png"
)

cat("\nENRICHMENT FIGURES CHECK\n")
print(
  list.files(
    "results/figures",
    pattern = "^(GO|KEGG)_.*\\.png$",
    full.names = FALSE
  )
)


############################################################
# 18. MACHINE LEARNING DATA PREPARATION
############################################################

# Samples are rows and probes are columns.

ml_data <- as.data.frame(t(normalized_expression))

ml_data$Group <- sample_group

# Remove probes with zero variance

feature_variance <- sapply(
  ml_data[, setdiff(names(ml_data), "Group"), drop = FALSE],
  var,
  na.rm = TRUE
)

valid_features <- names(feature_variance)[
  is.finite(feature_variance) & feature_variance > 0
]

ml_data <- ml_data[
  ,
  c(valid_features, "Group"),
  drop = FALSE
]

# Stratified 80:20 train-test split

set.seed(42)

train_index <- createDataPartition(
  ml_data$Group,
  p = 0.80,
  list = FALSE
)

train_data <- ml_data[train_index, , drop = FALSE]

test_data <- ml_data[-train_index, , drop = FALSE]

cat("\nTRAINING GROUP COUNTS\n")
print(table(train_data$Group))

cat("\nTEST GROUP COUNTS\n")
print(table(test_data$Group))


############################################################
# 19. FEATURE SELECTION FUNCTION
############################################################

select_top_features <- function(
    train_expression,
    train_group,
    number_of_features = 100
) {
  
  train_group <- factor(
    train_group,
    levels = c("Non-tumour", "Glioblastoma")
  )
  
  feature_design <- model.matrix(
    ~ 0 + train_group
  )
  
  colnames(feature_design) <- c(
    "Non_tumour",
    "Glioblastoma"
  )
  
  feature_contrast <- makeContrasts(
    Glioblastoma - Non_tumour,
    levels = feature_design
  )
  
  feature_fit <- lmFit(
    train_expression,
    feature_design
  )
  
  feature_fit <- contrasts.fit(
    feature_fit,
    feature_contrast
  )
  
  feature_fit <- eBayes(feature_fit)
  
  feature_results <- topTable(
    feature_fit,
    coef = 1,
    number = Inf,
    sort.by = "P"
  )
  
  selected_features <- rownames(
    feature_results
  )[seq_len(min(
    number_of_features,
    nrow(feature_results)
  ))]
  
  selected_features
}


############################################################
# 20. PREPARE TRAINING AND TEST FEATURES
############################################################

train_expression <- normalized_expression[
  ,
  train_index,
  drop = FALSE
]

selected_features <- select_top_features(
  train_expression,
  train_data$Group,
  number_of_features = 100
)

train_matrix <- t(
  normalized_expression[
    selected_features,
    train_index,
    drop = FALSE
  ]
)

test_matrix <- t(
  normalized_expression[
    selected_features,
    -train_index,
    drop = FALSE
  ]
)

# Rename features safely for modelling

colnames(train_matrix) <- paste0(
  "Probe_",
  seq_len(ncol(train_matrix))
)

colnames(test_matrix) <- colnames(train_matrix)


############################################################
# 21. SCALE TRAINING AND TEST DATA
############################################################

# Scaling parameters are calculated only from training data.

training_means <- colMeans(train_matrix)

training_sds <- apply(
  train_matrix,
  2,
  sd
)

training_sds[!is.finite(training_sds) |
               training_sds == 0] <- 1

train_scaled <- scale(
  train_matrix,
  center = training_means,
  scale = training_sds
)

test_scaled <- scale(
  test_matrix,
  center = training_means,
  scale = training_sds
)


############################################################
# 22. SUPPORT VECTOR MACHINE CLASSIFICATION
############################################################

svm_model <- svm(
  x = train_scaled,
  y = train_data$Group,
  kernel = "radial",
  cost = 1,
  gamma = 1 / ncol(train_scaled)
)

svm_predictions <- predict(
  svm_model,
  newdata = test_scaled
)

svm_predictions <- factor(
  svm_predictions,
  levels = levels(sample_group)
)

svm_confusion <- confusionMatrix(
  svm_predictions,
  test_data$Group,
  positive = "Glioblastoma"
)

cat("\nSVM TEST RESULTS\n")

print(svm_confusion)

write.csv(
  data.frame(
    Actual = test_data$Group,
    Predicted = svm_predictions
  ),
  "results/SVM_Test_Predictions.csv",
  row.names = FALSE
)


############################################################
# 23. RANDOM FOREST CLASSIFICATION
############################################################

rf_train_data <- data.frame(
  train_matrix,
  Group = train_data$Group,
  check.names = FALSE
)

rf_test_data <- data.frame(
  test_matrix,
  check.names = FALSE
)

set.seed(42)

rf_model <- randomForest(
  Group ~ .,
  data = rf_train_data,
  ntree = 500,
  importance = TRUE
)

rf_predictions <- predict(
  rf_model,
  newdata = rf_test_data
)

rf_predictions <- factor(
  rf_predictions,
  levels = levels(sample_group)
)

rf_confusion <- confusionMatrix(
  rf_predictions,
  test_data$Group,
  positive = "Glioblastoma"
)

cat("\nRANDOM FOREST TEST RESULTS\n")

print(rf_confusion)

write.csv(
  data.frame(
    Actual = test_data$Group,
    Predicted = rf_predictions
  ),
  "results/RandomForest_Test_Predictions.csv",
  row.names = FALSE
)


############################################################
# 24. RANDOM FOREST FEATURE IMPORTANCE
############################################################

rf_importance <- importance(rf_model)

rf_importance_table <- data.frame(
  Probe = rownames(rf_importance),
  rf_importance,
  row.names = NULL,
  check.names = FALSE
)

rf_importance_table <- rf_importance_table[
  order(
    rf_importance_table$MeanDecreaseGini,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

write.csv(
  rf_importance_table,
  "results/RandomForest_Feature_Importance.csv",
  row.names = FALSE
)

# Map selected probes back to their original IDs.

selected_probe_symbols <- mapIds(
  hgu133plus2.db,
  keys = selected_features,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

rf_importance_table$OriginalProbe <- selected_features[
  match(
    rf_importance_table$Probe,
    colnames(train_matrix)
  )
]

rf_importance_table$GeneSymbol <- as.character(
  selected_probe_symbols[
    rf_importance_table$OriginalProbe
  ]
)

write.csv(
  rf_importance_table,
  "results/RandomForest_Annotated_Genes.csv",
  row.names = FALSE
)

top20_rf <- head(
  rf_importance_table,
  20
)

write.csv(
  top20_rf,
  "results/Top20_Important_Genes.csv",
  row.names = FALSE
)

top20_rf$GeneLabel <- ifelse(
  is.na(top20_rf$GeneSymbol),
  top20_rf$OriginalProbe,
  top20_rf$GeneSymbol
)

importance_plot <- ggplot(
  top20_rf,
  aes(
    x = reorder(GeneLabel, MeanDecreaseGini),
    y = MeanDecreaseGini
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 Random Forest Important Genes",
    x = "Gene",
    y = "Mean Decrease Gini"
  ) +
  theme_minimal()

ggsave(
  "results/figures/Top20_Important_Genes.png",
  importance_plot,
  width = 9,
  height = 7
)


############################################################
# 25. MODEL EVALUATION FUNCTION
############################################################

calculate_metrics <- function(
    actual,
    predicted
) {
  
  actual <- factor(
    actual,
    levels = c("Non-tumour", "Glioblastoma")
  )
  
  predicted <- factor(
    predicted,
    levels = c("Non-tumour", "Glioblastoma")
  )
  
  cm <- table(
    Actual = actual,
    Predicted = predicted
  )
  
  TN <- cm["Non-tumour", "Non-tumour"]
  
  TP <- cm["Glioblastoma", "Glioblastoma"]
  
  FP <- cm["Non-tumour", "Glioblastoma"]
  
  FN <- cm["Glioblastoma", "Non-tumour"]
  
  accuracy <- (TP + TN) / sum(cm)
  
  sensitivity <- ifelse(
    TP + FN == 0,
    NA,
    TP / (TP + FN)
  )
  
  specificity <- ifelse(
    TN + FP == 0,
    NA,
    TN / (TN + FP)
  )
  
  precision <- ifelse(
    TP + FP == 0,
    NA,
    TP / (TP + FP)
  )
  
  f1_score <- ifelse(
    is.na(precision) ||
      is.na(sensitivity) ||
      precision + sensitivity == 0,
    NA,
    2 * precision * sensitivity /
      (precision + sensitivity)
  )
  
  data.frame(
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1_Score = f1_score
  )
}


############################################################
# 26. HOLDOUT MODEL COMPARISON
############################################################

svm_metrics <- calculate_metrics(
  test_data$Group,
  svm_predictions
)

rf_metrics <- calculate_metrics(
  test_data$Group,
  rf_predictions
)

holdout_comparison <- rbind(
  cbind(Model = "SVM", svm_metrics),
  cbind(Model = "Random Forest", rf_metrics)
)

print(holdout_comparison)

write.csv(
  holdout_comparison,
  "results/Holdout_Model_Comparison.csv",
  row.names = FALSE
)


############################################################
# 27. FIVE-FOLD STRATIFIED CROSS-VALIDATION
############################################################

# Feature selection and scaling are repeated separately
# inside every training fold to prevent data leakage.

set.seed(42)

folds <- createFolds(
  ml_data$Group,
  k = 5,
  list = TRUE,
  returnTrain = FALSE
)

svm_cv_predictions <- rep(
  NA_character_,
  nrow(ml_data)
)

rf_cv_predictions <- rep(
  NA_character_,
  nrow(ml_data)
)

for (i in seq_along(folds)) {
  
  cat("\nRunning CV Fold:", i, "\n")
  
  test_indices <- folds[[i]]
  
  train_indices <- setdiff(
    seq_len(nrow(ml_data)),
    test_indices
  )
  
  # Split fold data
  
  fold_train_group <- droplevels(
    ml_data$Group[train_indices]
  )
  
  fold_test_group <- ml_data$Group[test_indices]
  
  fold_train_expression <- t(
    normalized_expression[
      ,
      train_indices,
      drop = FALSE
    ]
  )
  
  fold_test_expression <- t(
    normalized_expression[
      ,
      test_indices,
      drop = FALSE
    ]
  )
  
  # Feature selection using training fold only
  
  fold_features <- select_top_features(
    t(fold_train_expression),
    fold_train_group,
    number_of_features = 100
  )
  
  fold_train_matrix <- fold_train_expression[
    ,
    fold_features,
    drop = FALSE
  ]
  
  fold_test_matrix <- fold_test_expression[
    ,
    fold_features,
    drop = FALSE
  ]
  
  # Safe feature names
  
  colnames(fold_train_matrix) <- paste0(
    "Probe_",
    seq_len(ncol(fold_train_matrix))
  )
  
  colnames(fold_test_matrix) <- colnames(
    fold_train_matrix
  )
  
  # Scaling from training fold only
  
  fold_means <- colMeans(fold_train_matrix)
  
  fold_sds <- apply(
    fold_train_matrix,
    2,
    sd
  )
  
  fold_sds[
    !is.finite(fold_sds) |
      fold_sds == 0
  ] <- 1
  
  fold_train_scaled <- scale(
    fold_train_matrix,
    center = fold_means,
    scale = fold_sds
  )
  
  fold_test_scaled <- scale(
    fold_test_matrix,
    center = fold_means,
    scale = fold_sds
  )
  
  # SVM model
  
  fold_svm <- svm(
    x = fold_train_scaled,
    y = fold_train_group,
    kernel = "radial",
    cost = 1,
    gamma = 1 / ncol(fold_train_scaled)
  )
  
  fold_svm_pred <- predict(
    fold_svm,
    newdata = fold_test_scaled
  )
  
  svm_cv_predictions[test_indices] <- as.character(
    fold_svm_pred
  )
  
  # Random Forest model
  
  fold_rf_train <- data.frame(
    fold_train_matrix,
    Group = fold_train_group,
    check.names = FALSE
  )
  
  fold_rf_test <- data.frame(
    fold_test_matrix,
    check.names = FALSE
  )
  
  set.seed(42 + i)
  
  fold_rf <- randomForest(
    Group ~ .,
    data = fold_rf_train,
    ntree = 500
  )
  
  fold_rf_pred <- predict(
    fold_rf,
    newdata = fold_rf_test
  )
  
  rf_cv_predictions[test_indices] <- as.character(
    fold_rf_pred
  )
  
}


############################################################
# 28. FIVE-FOLD CROSS-VALIDATION EVALUATION
############################################################

svm_cv_predictions <- factor(
  svm_cv_predictions,
  levels = levels(sample_group)
)

rf_cv_predictions <- factor(
  rf_cv_predictions,
  levels = levels(sample_group)
)

svm_cv_metrics <- calculate_metrics(
  ml_data$Group,
  svm_cv_predictions
)

rf_cv_metrics <- calculate_metrics(
  ml_data$Group,
  rf_cv_predictions
)

cv_comparison <- rbind(
  cbind(
    Model = "SVM",
    svm_cv_metrics
  ),
  cbind(
    Model = "Random Forest",
    rf_cv_metrics
  )
)

cat("\nFIVE-FOLD CROSS-VALIDATION RESULTS\n")

print(cv_comparison)

write.csv(
  cv_comparison,
  "results/Cross_Validation_Comparison.csv",
  row.names = FALSE
)

# Save all out-of-fold predictions

cv_predictions_table <- data.frame(
  Actual = ml_data$Group,
  SVM_Predicted = svm_cv_predictions,
  RF_Predicted = rf_cv_predictions
)

write.csv(
  cv_predictions_table,
  "results/Cross_Validation_Predictions.csv",
  row.names = FALSE
)


############################################################
# 29. CROSS-VALIDATION VISUALIZATION
############################################################

cv_plot_data <- cv_comparison %>%
  select(
    Model,
    Accuracy,
    Sensitivity,
    Specificity
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      Accuracy,
      Sensitivity,
      Specificity
    ),
    names_to = "Metric",
    values_to = "Score"
  )

cv_plot <- ggplot(
  cv_plot_data,
  aes(
    x = Model,
    y = Score,
    fill = Metric
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "SVM and Random Forest: Five-Fold CV",
    x = "Machine Learning Model",
    y = "Performance Score"
  ) +
  theme_minimal()

print(cv_plot)

ggsave(
  "results/figures/Cross_Validation_Comparison.png",
  cv_plot,
  width = 9,
  height = 6
)


############################################################
# 30. SAVE FINAL PROJECT OUTPUTS
############################################################

# Save important R objects for future use

saveRDS(
  list(
    normalized_expression = normalized_expression,
    deg_results = deg_results,
    significant_degs = significant_degs,
    annotated_degs = annotated_degs,
    significant_check = significant_check,
    annotated_significant = annotated_significant,
    deg_summary = deg_summary,
    top10_upregulated = top10_upregulated,
    top10_downregulated = top10_downregulated,
    go_up = go_up,
    go_down = go_down,
    kegg_up = kegg_up,
    kegg_down = kegg_down,
    svm_model = svm_model,
    rf_model = rf_model,
    cv_comparison = cv_comparison
  ),
  "results/Final_Analysis_Objects.rds"
)


############################################################
# 31. SAVE SESSION INFORMATION
############################################################

writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo.txt"
)


############################################################
# 32. FINAL PROJECT SUMMARY
############################################################

cat("\n")
cat("============================================\n")
cat(" GLIOBLASTOMA ANALYSIS COMPLETED\n")
cat("============================================\n")

cat("Total selected samples:",
    length(sample_group), "\n")

cat("Glioblastoma samples:",
    sum(sample_group == "Glioblastoma"), "\n")

cat("Non-tumour samples:",
    sum(sample_group == "Non-tumour"), "\n")

cat("Significant DEGs:",
    nrow(significant_degs), "\n")

cat("Upregulated probes:",
    nrow(upregulated_genes), "\n")

cat("Downregulated probes:",
    nrow(downregulated_genes), "\n")

cat("Annotated significant probes:",
    nrow(annotated_significant), "\n")

cat("Unannotated significant probes:",
    unannotated_significant_count, "\n")

cat("Unique annotated genes:",
    unique_annotated_gene_count, "\n")

cat("GO upregulated terms:",
    nrow(go_up_df), "\n")

cat("GO downregulated terms:",
    nrow(go_down_df), "\n")

cat("KEGG upregulated pathways:",
    nrow(kegg_up_df), "\n")

cat("KEGG downregulated pathways:",
    nrow(kegg_down_df), "\n")

cat("\nFIVE-FOLD CROSS-VALIDATION\n")

print(cv_comparison)

cat("\nAll available results saved in the results folder.\n")

############################################################
# 33. FINAL AUTOMATED VALIDATION AND CONFUSION MATRIX EXPORT
############################################################

cat("\n============================================\n")
cat(" FINAL AUTOMATED VALIDATION\n")
cat("============================================\n")

# A. Check DEG arithmetic and annotation counts.
stopifnot(
  nrow(upregulated_significant) +
    nrow(downregulated_significant) ==
    nrow(significant_check)
)

stopifnot(
  nrow(annotated_significant) +
    unannotated_significant_count ==
    nrow(significant_check)
)

stopifnot(
  unique_annotated_gene_count == length(
    unique(annotated_significant$GeneSymbol)
  )
)

# B. Check the top-10 gene tables.
stopifnot(
  nrow(top10_upregulated) == min(10, nrow(annotated_significant[annotated_significant$logFC > 1, , drop = FALSE])),
  nrow(top10_downregulated) == min(10, nrow(annotated_significant[annotated_significant$logFC < -1, , drop = FALSE]))
)

if (nrow(top10_upregulated) > 0) {
  stopifnot(all(top10_upregulated$logFC > 1))
}

if (nrow(top10_downregulated) > 0) {
  stopifnot(all(top10_downregulated$logFC < -1))
}

# C. Export out-of-fold cross-validation confusion matrices.
# These are based on the five-fold CV predictions, not the separate
# holdout test predictions.
if (exists("cv_predictions_table")) {
  svm_cm_table <- as.data.frame(
    table(
      Actual = cv_predictions_table$Actual,
      Predicted = cv_predictions_table$SVM_Predicted
    )
  )
  svm_cm_table$Model <- "SVM"
  
  rf_cm_table <- as.data.frame(
    table(
      Actual = cv_predictions_table$Actual,
      Predicted = cv_predictions_table$RF_Predicted
    )
  )
  rf_cm_table$Model <- "Random Forest"
  
  confusion_results <- rbind(svm_cm_table, rf_cm_table)
  
  write.csv(
    confusion_results,
    "results/Confusion_Matrices.csv",
    row.names = FALSE
  )
  
  cat("Five-fold CV confusion matrices saved.\n")
} else {
  warning("Cross-validation prediction table was not found; confusion matrix CSV was not created.")
}

# D. Check required output files after all exports.
required_output_files <- c(
  "results/DEG_Summary.csv",
  "results/annotated_DEG_results.csv",
  "results/Significant_Annotated_DEGs.csv",
  "results/Top10_Upregulated_Genes.csv",
  "results/Top10_Downregulated_Genes.csv",
  "results/Cross_Validation_Comparison.csv",
  "results/Cross_Validation_Predictions.csv",
  "results/Confusion_Matrices.csv",
  "results/GO_upregulated.csv",
  "results/GO_downregulated.csv",
  "results/KEGG_upregulated.csv",
  "results/KEGG_downregulated.csv",
  "results/sessionInfo.txt",
  "results/figures/PCA_Glioblastoma.png",
  "results/figures/Volcano_Plot_GBM.png",
  "results/figures/Top20_DEG_Heatmap.png",
  "results/figures/Cross_Validation_Comparison.png"
)

missing_output_files <- required_output_files[
  !file.exists(required_output_files)
]

if (length(missing_output_files) == 0) {
  cat("All required output files are present.\n")
} else {
  cat("The following required output files are missing:\n")
  print(missing_output_files)
}

# E. Save final validation status.
validation_summary <- data.frame(
  Check = c(
    "Significant probe count arithmetic",
    "Annotated/unannotated count arithmetic",
    "Unique annotated gene count",
    "Required output files"
  ),
  Status = c(
    "PASS",
    "PASS",
    "PASS",
    if (length(missing_output_files) == 0) "PASS" else "CHECK"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  validation_summary,
  "results/Final_Validation_Summary.csv",
  row.names = FALSE
)

cat("\nFINAL VALIDATION SUMMARY\n")
print(validation_summary)

cat("\n============================================\n")
cat(" END OF GLIOBLASTOMA ANALYSIS SCRIPT\n")
cat("============================================\n")

