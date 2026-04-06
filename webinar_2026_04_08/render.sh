#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

QMD="tcga_lung_analysis.qmd"

echo "Rendering HTML..."
quarto render "$QMD" --to html

echo "Done."
