# Gene Expression-Based Biomarker Discovery and Machine Learning Analysis of Glioblastoma
[View Glioblastoma Bioinformatics Project Report](https://github.com/Thilaiksan/2021-ASP-47/blob/main/Glioblastoma_Project_Report.html)

## Abstract

Glioblastoma (GBM) is an aggressive brain tumour. This project analyses
the public GEO dataset **GSE4290** using gene expression analysis,
functional enrichment and machine learning.

A total of **100 samples** were selected: 77 glioblastoma samples and 23
non-tumour samples. Following quality control and normalization,
differential expression analysis identified **7,210 significant probe
sets**, including 3,355 upregulated and 3,855 downregulated probe sets.
Gene Ontology (GO) and KEGG enrichment analyses were performed to
explore associated biological functions and pathways. Support Vector
Machine (SVM) and Random Forest (RF) classifiers were evaluated using
five-fold cross-validation, achieving accuracies of 0.92 and 0.91,
respectively.

These findings are exploratory and require independent validation before
clinical interpretation.

## 1. Introduction

Gene expression analysis can help identify genes with altered expression
and explore biological processes associated with glioblastoma.

### Aim

To identify differentially expressed genes associated with glioblastoma
and evaluate machine learning models for classifying glioblastoma and
non-tumour samples.

### Objectives

-   Obtain and inspect the GSE4290 dataset.
-   Select glioblastoma and non-tumour samples.
-   Perform quality control and normalization.
-   Identify differentially expressed genes.
-   Perform GO and KEGG enrichment analyses.
-   Develop and evaluate machine learning models.

## 2. Dataset Description

The dataset used in this project is **GSE4290**, obtained from the
National Center for Biotechnology Information (NCBI) Gene Expression
Omnibus (GEO).

-   **GEO dataset:**
    [GSE4290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4290)
-   **GEO platform:**
    [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
-   **Organism:** *Homo sapiens*

### Dataset summary

  Description                       Value
  ---------------------- ----------------
  GEO accession                   GSE4290
  Platform                         GPL570
  Organism                 *Homo sapiens*
  Initial samples                     180
  Selected samples                    100
  Glioblastoma samples                 77
  Non-tumour samples                   23
  Initial probe sets               54,613

**Note:** The comparison group is described as *non-tumour* and is not
automatically considered a group of healthy controls.

## 3. Materials and Methods

### 3.1 Data preprocessing

The expression matrix was filtered to include the selected 100 samples.
Quality control checks were performed for missing and infinite values.

Expression values were log2-transformed and quantile-normalized using
the project preprocessing workflow.

### 3.2 Exploratory data analysis

The project generated the following visualizations:

-   Normalized expression boxplot
-   Principal Component Analysis (PCA)
-   Differential expression volcano plot
-   Heatmap of the top 20 differentially expressed probes

### 3.3 Differential expression analysis

The `limma` package was used to identify differentially expressed probe
sets.

The significance criteria were:

-   Adjusted p-value \< 0.05
-   Absolute log fold change \> 1

### 3.4 Functional enrichment

Gene Ontology (GO) and KEGG pathway enrichment analyses were conducted
separately for upregulated and downregulated genes.

### 3.5 Machine learning

Two classification algorithms were evaluated:

-   Support Vector Machine (SVM)
-   Random Forest (RF)

Five-fold cross-validation was used as the primary model evaluation
approach.

## 4. Results

### 4.1 Quality control

  Measure                            Result
  ------------------------------- ---------
  Selected samples                      100
  Missing values                          0
  Infinite values                         0
  Normalized expression minimum      0.2637
  Normalized expression maximum     15.1306

### 4.2 Differential gene expression

A total of **7,210 significant probe sets** were identified.

  Category                         Number
  ------------------------------ --------
  Upregulated probe sets            3,355
  Downregulated probe sets          3,855
  Total significant probe sets      7,210
  Annotated probe sets              6,490
  Unannotated probe sets              720
  Unique annotated genes            4,534

Examples of the top upregulated genes include **TOP2A, ASPM, NDC80,
IGFBP2, RRM2, LTF, NUSAP1, COL4A1** and **EGFR**.

Examples of the top downregulated genes include **GABRA1, SLC60A1,
CACNG3, SV2B, GJB6, CCK, SLITRK4, VSNL1, DDN** and **GABRA5**.

### 4.3 GO enrichment analysis

The analysis identified:

-   1,001 GO terms for upregulated genes.
-   631 GO terms for downregulated genes.

### 4.4 KEGG pathway enrichment

The analysis identified:

-   58 enriched KEGG pathways for upregulated genes.
-   71 enriched KEGG pathways for downregulated genes.

### 4.5 Machine learning results

Five-fold cross-validation produced the following results:

  Metric             SVM   Random Forest
  ------------- -------- ---------------
  Accuracy          0.92            0.91
  Sensitivity     0.9351          0.9351
  Specificity     0.8696          0.8261
  Precision       0.9600          0.9474
  F1-score        0.9474          0.9412

SVM achieved a slightly higher cross-validation accuracy than Random
Forest in this analysis.

## 5. Discussion

The analysis identified differentially expressed genes between
glioblastoma and non-tumour samples. TOP2A, ASPM, NDC80 and RRM2 were
among the top upregulated results.

GO and KEGG enrichment analyses provided additional information about
biological functions and pathways associated with the selected gene
lists.

The SVM model achieved a slightly higher cross-validation accuracy than
Random Forest. However, these results should be interpreted cautiously
because the dataset contains only 100 samples and the sample groups are
imbalanced.

## 6. Limitations

-   The dataset contains a relatively small number of samples.
-   The sample groups are imbalanced.
-   The non-tumour group should not automatically be considered healthy
    controls.
-   No independent external dataset was used for validation.
-   Candidate genes are not clinically validated biomarkers.
-   Machine learning performance may differ on new datasets.

## 7. Conclusion and Future Work

This project developed a bioinformatics workflow to analyse glioblastoma
gene expression data. Differential expression analysis identified 7,210
significant probe sets. GO and KEGG enrichment analyses provided
functional information, while machine learning models achieved five-fold
cross-validation accuracies of 0.92 for SVM and 0.91 for Random Forest.

Future work should validate candidate genes using independent datasets,
investigate additional machine learning methods and assess model
performance on external samples.

## 8. Project Structure

The main R analysis script is located at `code/Glioblastoma_Analysis.R`.
Generated tables, figures and analysis outputs are stored in the
`results/` directory.

``` text
Glioblastoma_Bioinformatics_Project/
├── code/
│   └── Glioblastoma_Analysis.R
├── data/
├── results/
│   └── figures/
├── Glioblastoma_Project_Report.Rmd
├── Glioblastoma_Project_Report.html
└── README.md
```

The directory listing above summarizes the main project components; the
`data/` and `results/` folders contain the project files and generated
outputs.

## 9. Reproducibility

The main analysis script is:

``` text
code/Glioblastoma_Analysis.R
```

From the project root directory in R or RStudio, run:

``` r
source("code/Glioblastoma_Analysis.R")
```

The analysis depends on the project data files and required R packages
being available. The generated tables, figures and analysis outputs are
stored in the `results/` directory.

## 10. References

1.  Edgar, R., Domrachev, M., and Lash, A. E. (2002). Gene Expression
    Omnibus: NCBI gene expression and hybridization array data
    repository. *Nucleic Acids Research*, 30(1), 207--210.
2.  GEO Series GSE4290.
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4290
3.  GEO Platform GPL570.
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
4.  Ritchie, M. E., et al. (2015). limma powers differential expression
    analyses for RNA-sequencing and microarray studies. *Nucleic Acids
    Research*, 43(7), e47.
5.  The Gene Ontology Consortium. (2021). The Gene Ontology resource.
    *Nucleic Acids Research*, 49(D1), D325--D334.
6.  Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5--32.
7.  Cortes, C., and Vapnik, V. (1995). Support-vector networks. *Machine
    Learning*, 20, 273--297.
