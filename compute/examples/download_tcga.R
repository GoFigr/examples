## This script downloads TCGA-LUAD and TCGA-LUSC RNA-Seq data from GDC and
## subsets to curated cancer genes + top 1000 expressed genes. If you'd rather
## skip this step, a pre-built RDS file is available at:
##   https://cdn.gofigr.io/data/TCGA-LUNG_TPM_SE_subset.rds

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(arrow)

## --- Download / cache each cohort -------------------------------------------
download_tcga_cohort <- function(project, cache_file) {
  if (file.exists(cache_file)) {
    cat(sprintf("Loading cached %s from %s\n", project, cache_file))
    return(readRDS(cache_file))
  }

  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    experimental.strategy = "RNA-Seq"
  )

  GDCdownload(query)
  se <- GDCprepare(query, save = FALSE)
  saveRDS(se, cache_file, compress = FALSE)
  return(se)
}

luad_se <- download_tcga_cohort("TCGA-LUAD", "TCGA-LUAD_TPM_SE.rds")
lusc_se <- download_tcga_cohort("TCGA-LUSC", "TCGA-LUSC_TPM_SE.rds")

## --- Merge LUAD + LUSC into a single SummarizedExperiment -------------------
## Keep only genes present in both cohorts
shared_genes <- intersect(rownames(luad_se), rownames(lusc_se))
cat(sprintf("Shared genes between LUAD and LUSC: %d\n", length(shared_genes)))

luad_se <- luad_se[shared_genes, ]
lusc_se <- lusc_se[shared_genes, ]

## Add a cohort column before merging
luad_se$cohort <- "LUAD"
lusc_se$cohort <- "LUSC"

## Harmonize colData columns (LUAD and LUSC may differ) before combining
shared_cols <- intersect(colnames(colData(luad_se)), colnames(colData(lusc_se)))
colData(luad_se) <- colData(luad_se)[, shared_cols]
colData(lusc_se) <- colData(lusc_se)[, shared_cols]

## Combine — cbind works on SummarizedExperiment when rows match
lung_se <- cbind(luad_se, lusc_se)
cat(sprintf("Combined: %d samples (%d LUAD, %d LUSC), %d genes\n",
            ncol(lung_se), ncol(luad_se), ncol(lusc_se), nrow(lung_se)))

## ---------------------------------------------------------------------------
## Expanded gene list for LUAD / pan-cancer demos (~250 genes)
## Sources:
##   - TCGA LUAD marker paper (Cancer Genome Atlas Research Network, Nature 2014)
##   - COSMIC Cancer Gene Census (Sondka et al., Nat Rev Cancer 2018)
##   - Hallmarks of Cancer gene sets (Hanahan & Weinberg 2011)
##   - MSigDB Hallmark / C6 oncogenic signatures
##   - vanHelden et al. (2023) immune checkpoint landscape
##   - Housekeeping genes: Eisenberg & Levanon, Trends Genet 2013
## ---------------------------------------------------------------------------

# --- LUAD key drivers & frequently mutated genes ---
luad_drivers <- c(
  "EGFR", "KRAS", "TP53", "ALK", "ROS1", "BRAF", "MET", "ERBB2",
  "STK11", "KEAP1", "NF1", "RIT1", "MAP2K1", "PIK3CA", "SMARCA4",
  "RB1", "CDKN2A", "NKX2-1", "SETD2", "ARID1A", "U2AF1", "RBM10",
  "ATM", "CTNNB1", "APC", "SMAD4", "FGFR1", "FGFR2", "FGFR3",
  "DDR2", "NTRK1", "NTRK2", "NTRK3", "RET", "ERBB3", "ERBB4",
  "HRAS", "NRAS", "MAP2K2", "ARAF", "RAF1"
)

# --- Oncogenes & tumor suppressors (pan-cancer, COSMIC CGC) ---
oncogenes_tsg <- c(
  "MYC", "MYCN", "MYCL", "JUN", "FOS", "CCND1", "CCNE1", "CDK4",
  "CDK6", "MDM2", "MDM4", "BCL2", "BCL2L1", "MCL1", "BIRC5",
  "TERT", "IDH1", "IDH2", "NOTCH1", "NOTCH2", "NOTCH3",
  "PTEN", "VHL", "BRCA1", "BRCA2", "PALB2", "RAD51",
  "MLH1", "MSH2", "MSH6", "PMS2",
  "AKT1", "AKT2", "MTOR", "TSC1", "TSC2",
  "FBXW7", "KMT2A", "KMT2D", "KDM6A", "DNMT3A", "TET2",
  "EZH2", "SUZ12", "ARID2", "PBRM1", "BAP1",
  "SRC", "ABL1", "JAK2", "STAT3", "STAT5B"
)

# --- Signaling pathways (RTK/RAS/MAPK, PI3K/AKT, Wnt, Hippo, Notch) ---
signaling <- c(
  "SOS1", "GRB2", "SHC1", "GAB1", "PTPN11",
  "MAPK1", "MAPK3", "MAP3K1", "DUSP6",
  "PIK3CB", "PIK3R1", "RICTOR", "RPTOR",
  "GSK3B", "AXIN1", "LEF1", "TCF7L2",
  "YAP1", "WWTR1", "LATS1", "LATS2", "MST1", "STK3",
  "DLL1", "DLL4", "JAG1", "HES1", "HEY1",
  "SHH", "SMO", "PTCH1", "GLI1", "GLI2",
  "WNT3A", "WNT5A", "WNT7A", "FZD7", "LRP5", "LRP6"
)

# --- DNA damage response & cell cycle ---
ddr_cellcycle <- c(
  "CHEK1", "CHEK2", "ATR", "PARP1", "PARP2",
  "CDK1", "CDK2", "CDKN1A", "CDKN1B", "CDKN2B",
  "CCNA2", "CCNB1", "PLK1", "AURKA", "AURKB",
  "BUB1", "MAD2L1", "TTK", "CDC20", "CENPE",
  "E2F1", "E2F3", "FOXM1", "MKI67", "TOP2A"
)

# --- Apoptosis & autophagy ---
apoptosis <- c(
  "BAX", "BAK1", "BID", "BIM", "CASP3", "CASP8", "CASP9",
  "FASLG", "FAS", "TNFRSF10A", "TNFRSF10B",
  "XIAP", "CFLAR", "BECN1", "ATG5", "ATG7", "SQSTM1"
)

# --- Tumor microenvironment & immune markers ---
immune_tme <- c(
  # Immune checkpoints
  "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
  "TIGIT", "ICOS", "CD28", "IDO1", "VTCN1", "VISTA",
  # T-cell / NK markers
  "CD3D", "CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "GZMA", "GZMB",
  "PRF1", "IFNG", "TNF", "IL2", "NCAM1", "KLRK1", "NKG7",
  # Myeloid / macrophage
  "CD68", "CD163", "CSF1R", "ITGAM", "ARG1", "NOS2",
  # B-cell
  "CD19", "MS4A1", "CD79A",
  # Cytokines & chemokines
  "CXCL9", "CXCL10", "CXCL13", "CCL2", "CCL5",
  "IL6", "IL10", "TGFB1", "TGFB2", "VEGFA", "VEGFB"
)

# --- Lung-specific / differentiation markers ---
lung_markers <- c(
  "SFTPA1", "SFTPB", "SFTPC", "SFTPD",
  "SCGB1A1", "MUC1", "MUC5AC", "MUC5B",
  "KRT5", "KRT7", "KRT14", "KRT18", "KRT19",
  "NAPSA", "TTF1", "TP63", "SOX2", "SOX9"
)

# --- EMT markers ---
emt <- c(
  "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2",
  "TWIST1", "ZEB1", "ZEB2", "FN1", "MMP2", "MMP9"
)

# --- Metabolism (Warburg, lipid, amino acid) ---
metabolism <- c(
  "HK2", "PKM", "LDHA", "SLC2A1", "GLS",
  "FASN", "ACLY", "SLC7A11", "GPX4", "NFE2L2"
)

# --- Housekeeping / reference genes ---
housekeeping <- c(
  "ACTB", "GAPDH", "B2M", "RPL13A", "RPLP0",
  "HPRT1", "HMBS", "SDHA", "TBP", "UBC",
  "YWHAZ", "RPS18", "RPL27", "EEF1A1", "PPIA"
)

target_genes <- unique(c(
  luad_drivers, oncogenes_tsg, signaling, ddr_cellcycle,
  apoptosis, immune_tme, lung_markers, emt, metabolism, housekeeping
))

cat(sprintf("Curated cancer gene list: %d genes\n", length(target_genes)))

## Also include the top 1000 most highly expressed genes (by mean TPM across
## all samples) to provide background variance for unsupervised analyses
## (PCA, UMAP, clustering). Without this, PCA on cancer-only genes just
## captures variance among hand-picked features.
tpm_mat <- assay(lung_se, "tpm_unstrand")
gene_mean_expr <- rowMeans(tpm_mat, na.rm = TRUE)
top_1000 <- names(sort(gene_mean_expr, decreasing = TRUE))[1:1000]
top_1000_symbols <- rowData(lung_se[top_1000, ])$gene_name

all_genes <- unique(c(target_genes, top_1000_symbols))
cat(sprintf("Combined gene list (curated + top 1000 expressed): %d unique genes\n",
            length(all_genes)))

rows_to_keep <- rowData(lung_se)$gene_name %in% all_genes
lung_se_subset <- lung_se[rows_to_keep, ]

cat(sprintf("Matched %d / %d genes in combined LUAD/LUSC data\n",
            length(unique(rowData(lung_se_subset)$gene_name)),
            length(all_genes)))

dir.create("data", showWarnings = FALSE)
saveRDS(lung_se_subset, "data/TCGA-LUNG_TPM_SE_subset.rds")

## --- Export tumor expression + clinical as parquet for the Python notebook --
## The downstream Python classifier reads these directly, so we precompute them
## here rather than in the Qmd.
expr_raw <- assay(lung_se_subset, "tpm_unstrand")
expr <- log2(expr_raw + 1)
rownames(expr) <- make.unique(rowData(lung_se_subset)$gene_name)

clinical <- as.data.frame(colData(lung_se_subset))
tumor_idx <- clinical$sample_type == "Primary Tumor"
expr_tumor <- expr[, tumor_idx]
clin_tumor <- clinical[tumor_idx, ]

cat(sprintf("Tumor samples: %d (%d LUAD, %d LUSC) | Genes: %d\n",
            ncol(expr_tumor),
            sum(clin_tumor$cohort == "LUAD"),
            sum(clin_tumor$cohort == "LUSC"),
            nrow(expr_tumor)))


# Export to Parquet for use in Python
write_parquet(
  as.data.frame(t(expr_tumor)),
  "data/expr_tumor.parquet"
)
write_parquet(
  data.frame(
    sample_id = colnames(expr_tumor),
    cohort = clin_tumor$cohort,
    gender = clin_tumor$gender,
    stage = clin_tumor$ajcc_pathologic_stage,
    stringsAsFactors = FALSE
  ),
  "data/clinical_tumor.parquet"
)
