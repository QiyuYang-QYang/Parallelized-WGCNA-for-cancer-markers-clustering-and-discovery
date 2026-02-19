# WGCNA-PLUMBER  
**A High-Performance Weighted Gene Co-Expression Network Analysis Pipeline**

**Authors:** Lily Li, Charlie Yang, You Li, River Zhu  
**Origin:** Programming for Scientists (02-601), Carnegie Mellon University  

---

## Project Overview

Weighted Gene Co-Expression Network Analysis (WGCNA) is a widely adopted framework for identifying gene modules and hub genes from high-dimensional transcriptomic datasets. Standard R-based implementations, while statistically robust, often encounter runtime and memory constraints when scaling to tens of thousands of genes.

This project presents **WGCNA-PLUMBER**, a hybrid computational pipeline designed to improve scalability by:

- Re-engineering computationally intensive matrix operations in **Go**
- Preserving clustering, visualization, and biological interpretation in **R**
- Enabling efficient co-expression network analysis for large RNA-seq datasets

The pipeline was applied to thyroid transcriptomic data for module discovery and hub gene identification.

---

## Data Sources

### RNA-seq Expression Data

Initial exploratory analysis:

- **GTEx Bulk Tissue Expression Data**  
  https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression  
  Subset: *thyroid gene read counts*

> Note: Early experiments were restricted to thyroid samples.

---

### Gene Annotation (TPM Normalization)

- **GENCODE Human Release 44**  
  https://www.gencodegenes.org/human/release_44.html  

Used for:

- Accurate gene length calculation  
- TPM normalization

---

### Updated Dataset (Revision 11.26)

Following refinement inspired by thyroid cancer literature:

- TCGA-THCA dataset  
- Manual selection of tumor samples

---

## Pipeline Architecture

The pipeline consists of five computational phases:

---

### Phase 1 — Data Cleaning & Normalization

**Inputs:**
- `gencode.v44.annotation.gtf.gz`
- `gene_reads_v10_thyroid.gct.gz`

**Steps:**

1. **Gene Length Extraction**  
   `gtf_parser.go` parses GTF annotations to compute gene lengths.

2. **Two-Pass TPM Normalization** (`gct_processor.go`)  
   - Pass 1: Compute TPM denominator (`perSampleRPKSum`)  
   - Pass 2:  
     - Convert read counts → TPM  
     - Apply `log2(TPM + 1)` transformation

3. **Gene Filtering**
   - Low-expression filtering  
   - Low-variance filtering  

**Output:**  
`clean_thyroid_matrix.csv`

---

### Phase 2 — Correlation Matrix Computation

**Steps:**

1. Precompute gene-wise mean and standard deviation  
2. Parallelized Pearson correlation computation using Go worker pools  
3. ~1.4 × 10⁸ gene-pair correlations evaluated  

**Output:**  
`correlation_matrix.csv`

---

### Phase 3 — Adjacency Matrix Construction

Soft thresholding:

\[
adjacency_{ij} = |correlation_{ij}|^{\beta}, \quad \beta = 6
\]

**Purpose:**
- Amplify strong correlations  
- Approximate scale-free topology  

**Output:**  
`adjacency_matrix.csv`

---

### Phase 4 — Topological Overlap Matrix (TOM)

TOM quantifies shared network neighborhood structure:

\[
TOM_{ij} =
\frac{\sum_{u \neq i,j} a_{iu} a_{ju} + a_{ij}}
{\min(k_i, k_j) + 1 - a_{ij}}
\]

**Purpose:**
- Stabilize network similarity  
- Capture indirect gene relationships  

**Output:**  
`tom_matrix.csv`

---

### Phase 5 — Dissimilarity Matrix

\[
dissimilarity_{ij} = 1 - TOM_{ij}
\]

**Purpose:**
- Convert similarity → distance  
- Prepare for clustering  

**Output:**  
`dissimilarity_matrix.csv`

---

## Downstream Analysis (R)

Performed in R:

- Hierarchical clustering  
- Dynamic Tree Cut  
- Module merging via eigengene correlation  
- Hub gene identification (kME)  
- GO / KEGG enrichment  

Visualization:

- Gene dendrogram  
- Module eigengene heatmap  
- Module-trait relationships  
- Cytoscape hub networks  
- RShiny interactive interface  

---

## Key Results

- ~50 initial modules detected  
- Merged into **15 consensus modules**

### Biological Insights

- Significant module–trait associations  
- ECM-related modules identified  

### Hub Gene Discovery

- **PODN** identified as a central hub gene  
- Consistent with literature on ECM remodeling and tumor progression

---

## Computational Contributions

- Parallelized Pearson correlation computation in Go  
- Scalable adjacency and TOM matrix calculations  
- Reduced runtime compared to sequential implementations  
- Hybrid Go + R analytical architecture  

---

## Biological Contributions

- System-level module detection  
- Hub gene identification  
- Functional enrichment analysis  
- Hypothesis generation for thyroid cancer biology  

---

## How to Reproduce

```bash
go run main.go
