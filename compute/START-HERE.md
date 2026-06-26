# Welcome to your GoFigr workspace 👋

This is your own cloud machine, running JupyterLab, VS Code, and R — with the
GoFigr client **already installed and signed in**. Anything you save here stays
here: your files persist across restarts.

## Publish your first figure

You don't need to set up credentials — this box is already connected to your
GoFigr account. In a notebook cell:

```python
%load_ext gofigr

configure(analysis=FindByName("My First Analysis", create=True))

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([0, 1, 2, 3], [0, 1, 4, 9])
ax.set_title("Hello, GoFigr")

publish(fig, target="My first figure")
```

Run it, then open [gofigr.io](https://gofigr.io) — your figure is there, with the
code, data, and environment that produced it captured alongside it. That's the
whole idea: every figure you publish is reproducible and traceable back to its
source.

**Working in R?** `library(gofigR)` is installed too — see the
`tcga_lung_analysis.qmd` example.

## Explore the examples

The `examples/` folder next to this file has a complete, runnable analysis — a
TCGA lung-cancer classifier (Python) and report (R) — showing GoFigr used in a
real workflow. Open `examples/tcga_lung_classifier.ipynb` and run it top to
bottom.

## Good to know

- **Your work persists.** Files in your home directory survive stop/start, so
  you can pick up exactly where you left off.
- **The machine auto-stops when idle** to save you money. Nothing is lost — just
  start it again from the GoFigr app and your files are waiting.
- **It's pre-loaded.** Common Python and R packages are installed; add your own
  with `pip install` / `install.packages()` and they'll persist too.

## Learn more

- Getting started: <https://gofigr.io/getting-started-with-gofigr/>
- Documentation: <https://gofigr.io/docs>
- Python client guide: <https://gofigr.io/docs/gofigr-python/latest/start.html>

Questions? Email us at [support@gofigr.io](mailto:support@gofigr.io).
