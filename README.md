# GoFigr examples

Example analyses demonstrating how to use [GoFigr](https://gofigr.io) — a
zero-effort figure repository for Jupyter and R that automatically captures,
versions, and indexes every figure you produce.

## Webinars

### `webinar_2026_04_08/` — TCGA Lung Cancer (LUAD vs LUSC)

A multi-language analysis of TCGA-LUAD and TCGA-LUSC RNA-Seq data, used to
demonstrate GoFigr across both R and Python in the same project.

- **`download_tcga.R`** — Downloads TCGA-LUAD and TCGA-LUSC cohorts from GDC
  via `TCGAbiolinks`, merges them, subsets to a curated cancer gene list plus
  the top 1000 highly-expressed genes, and exports `data/expr_tumor.parquet`
  and `data/clinical_tumor.parquet` for downstream Python use.
- **`tcga_lung_analysis.qmd`** — Quarto/R analysis: PCA, top PC1 loadings,
  LUSC-vs-LUAD volcano plot, UMAP, interactive k-means clustering, and
  Kaplan-Meier survival analysis by proliferation score and NKX2-1 (TTF-1).
  Demonstrates `gofigR::enable`, `sync_file`, `publish`, and `reproducible()`
  for interactive re-stratification in the GoFigr portal.
- **`tcga_lung_classifier.ipynb`** — Python follow-up: trains an L1-regularized
  logistic regression classifier on the parquet exports, plots top predictive
  genes, an interactive ROC curve (re-stratifiable by gender/stage via
  `@reproducible`), and compares ML feature importance against the R volcano
  plot's differential expression.
- **`render.sh`** — Renders the Quarto document to HTML/PDF.

Pre-built parquet inputs are also hosted at
`https://cdn.gofigr.io/data/expr_tumor.parquet` and
`https://cdn.gofigr.io/data/clinical_tumor.parquet` if you'd rather skip the
R download step.
