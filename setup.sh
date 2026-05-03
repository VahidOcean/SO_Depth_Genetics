#!/usr/bin/env bash
##############################################################################
##  SETUP & PRE-FLIGHT CHECK
##  Run this before the R pipeline to verify all external dependencies
##############################################################################

echo "========================================"
echo "  Southern Ocean COI Pipeline Setup"
echo "========================================"

# ── Check R packages
echo ""
echo "[1/4] Checking R..."
Rscript -e "cat('R version:', R.version\$major, '.', R.version\$minor, '\n')"

# ── Check MAFFT
echo ""
echo "[2/4] Checking MAFFT..."
if command -v mafft &> /dev/null; then
    mafft --version 2>&1 | head -1
    echo "  OK: MAFFT found"
else
    echo "  MISSING: MAFFT not found"
    echo "  Install:"
    echo "    macOS:  brew install mafft"
    echo "    Ubuntu: sudo apt install mafft"
    echo "    conda:  conda install -c bioconda mafft"
fi

# ── Check geosphere (needed for geographic distances in Mantel)
echo ""
echo "[3/4] Checking optional R packages..."
Rscript -e "
pkgs <- c('geosphere','gridExtra')
missing <- pkgs[!pkgs %in% installed.packages()[,'Package']]
if (length(missing) > 0) {
    cat('Installing:', paste(missing, collapse=', '), '\n')
    install.packages(missing, repos='https://cloud.r-project.org')
} else {
    cat('  OK: geosphere, gridExtra present\n')
}
"

# ── Create directory structure
echo ""
echo "[4/4] Creating directory structure..."
mkdir -p data/raw data/aligned data/filtered \
         results/tables results/figures results/networks
echo "  OK: Directories created"
echo ""
echo "Setup complete. Run the pipeline with:"
echo "  Rscript SO_depth_genetics_pipeline.R"
echo ""
