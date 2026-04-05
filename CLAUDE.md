# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains example Jupyter notebooks demonstrating how to use GoFigr — a zero-effort figure repository for Jupyter and R that automatically syncs and indexes figures.

## Repository Contents

- **GoFigr_demo.ipynb** — Introductory demo showing GoFigr basics: setup, automatic figure publishing, dataframe attachment, and Plotly support. Designed for Google Colab.
- **hlca_lung_cpu_analysis.ipynb** — Single-cell RNA-seq analysis (Scanpy) with GoFigr integration, demonstrating usage in a real scientific workflow.

## GoFigr Integration Pattern

Notebooks follow this pattern:
1. `%load_ext gofigr` to load the extension
2. `from gofigr.jupyter import configure, publish, FindByName`
3. `configure(analysis=FindByName("Name", create=True))` to set the target analysis
4. Figures are auto-published; use `publish(fig, dataframes={...})` for explicit control with attached data

## Running Notebooks

Requires `pip install gofigr`. Account setup via `gfconfig` CLI tool. The demo notebook is designed for Google Colab (includes Drive mounting); the RAPIDS notebook needs `scanpy` and related bioinformatics packages.
