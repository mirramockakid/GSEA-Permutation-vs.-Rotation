# GSEA-Permutation-vs.-Rotation

The approach used to determine significance of a GSEA test:
  - (i) create a test statistic for the observed data
  - (ii) randomise some aspect of the original dataset that breaks the association
  - (iii) a new tests statistic

## Quarto

This repo now uses Quarto documents (`.qmd`) instead of R Markdown (`.Rmd`).

Render the walkthrough report:

```bash
quarto render gsea_rotation_vs_permutation.qmd
```

Local renders write output to `report/` (see `_quarto.yml`).

Render all project documents (outputs are written to `public/` via `_quarto.yml`):

```bash
quarto render
```

GitHub Actions uses the `github` Quarto profile (`_quarto-github.yml`) so CI output is written to `public/` for Pages deployment.
