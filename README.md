# Covariate-adjusted statistical dependence representation through partial copulas

This repository contains the code and manuscript for the paper *"Covariate-adjusted statistical dependence representation through partial copulas: bounds and new insights"* by Vinícius Litvinoff Justus and Felipe Fontana Vieira.

## Repository structure

```
.
├── manuscript.Rmd          # Main manuscript (text, figures, and simulation code)
├── simulation_and_plots.R  # Standalone simulation and plotting script
├── references.bib          # BibTeX references
├── gen_env.R               # Script to generate the Nix environment (uses {rix})
├── default.nix             # Nix expression defining the reproducible environment
├── .github/workflows/
│   └── render-paper.yml    # CI workflow to render the manuscript on push
└── copula_project.Rproj    # RStudio project file
```

## Reproducing the manuscript

The project uses [Nix](https://nixos.org/) to guarantee a fully reproducible environment (R version, packages, and LaTeX distribution). No manual installation of R packages or LaTeX is required.

### Prerequisites

Install Nix by following the instructions at https://nixos.org/download or using the [Determinate Systems installer](https://github.com/DeterminateSystems/nix-installer):

```bash
curl --proto '=https' --tlsv1.2 -sSf -L https://install.determinate.systems/nix | sh -s -- install
```

### Steps

1. **Clone the repository:**

   ```bash
   git clone https://github.com/felipelfv/copula_project.git
   cd copula_project
   ```

2. **Enter the Nix shell** (this downloads and builds all dependencies on first run):

   ```bash
   nix-shell
   ```

3. **Render the manuscript:**

   ```bash
   Rscript -e "rmarkdown::render('manuscript.Rmd')"
   ```

   This produces `manuscript.pdf` in the project root.

### Regenerating the Nix environment

If you need to modify R or LaTeX dependencies, edit `gen_env.R` and run:

```bash
Rscript gen_env.R
```

This regenerates `default.nix` using the [{rix}](https://docs.ropensci.org/rix/) package. Then exit and re-enter `nix-shell` for the changes to take effect.

## Continuous integration

The GitHub Actions workflow (`.github/workflows/render-paper.yml`) automatically renders the manuscript on every push to `main` and commits the updated PDF.
