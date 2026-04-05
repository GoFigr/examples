#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

QMD="tcga_luad_analysis.qmd"

echo "Rendering HTML..."
quarto render "$QMD" --to html

echo "Rendering PDF..."
quarto render "$QMD" --to pdf

echo "Done."
