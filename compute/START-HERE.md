# Welcome to your GoFigr workspace 👋

This is your own cloud machine, running JupyterLab, a VS Code–compatible
editor (code-server), and R Server — with the GoFigr client **already
installed and signed in**. Anything you save here stays here: your files
persist across restarts.

## Try it — publish your first figure

You don't need to set up credentials: this machine is already connected to your
GoFigr account. The quickest way to see it work is the bundled example — open
**[`examples/tcga_lung_classifier.ipynb`](examples/tcga_lung_classifier.ipynb)**
and run it top to bottom.

As it runs, each figure is published to your GoFigr account together with the
code, data, and environment that produced it — so every figure is reproducible
and traceable back to its source. Open
[app.gofigr.io](https://app.gofigr.io) to watch them appear.

**Working in R?** `library(gofigR)` is installed too — see
**[`examples/tcga_lung_analysis.qmd`](examples/tcga_lung_analysis.qmd)**, the
same analysis as a Quarto report.

## Ask AI as you work

Jupyter AI is preconfigured and connected to our managed LLM gateway — no API
keys needed, and free AI credits are included with your account. Open the chat
panel in JupyterLab and ask about your code, your data, or that error you just
hit. Prefer working in cells? Run `%load_ext jupyter_ai_magics` once and use
`%%ai` directly in your notebook.

## Good to know

- **Your work persists.** Files in your home directory survive stop/start, so
  you can pick up exactly where you left off.
- **The machine auto-stops when idle** to save you money. Nothing is lost — just
  start it again from the GoFigr app and your files are waiting.
- **It's pre-loaded.** Common Python and R packages are installed; add your own
  with `pip install` / `install.packages()` and they'll persist too.

## Learn more

- Working in your instance: <https://docs.gofigr.io/compute/working-in-the-workspace>
- Instance lifecycle (stop/start, idle auto-stop): <https://docs.gofigr.io/compute/lifecycle>
- Documentation: <https://docs.gofigr.io>

Questions? Email us at [support@gofigr.io](mailto:support@gofigr.io).
