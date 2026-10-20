#!/bin/bash
set -eo pipefail     # no -u: the module system uses unset variables internally

PROJECT=$SCRATCH/thesis        # adapt
cd "$PROJECT"

# ---- internet (proxy for shell tools and for Java) ----
module purge
module load eth_proxy
if [ -n "${https_proxy:-}" ]; then
  p=${https_proxy#*://}; p=${p%%/*}
  export JAVA_TOOL_OPTIONS="-Dhttps.proxyHost=${p%%:*} -Dhttps.proxyPort=${p##*:} -Dhttp.proxyHost=${p%%:*} -Dhttp.proxyPort=${p##*:}"
fi

# ---- BEAST packages ----
module load stack/2024-06 gcc/12.2.0 openjdk beast2 libbeagle
packagemanager -add BDMM-Prime
packagemanager -add feast
packagemanager -add SA -version 2.1.1
unset JAVA_TOOL_OPTIONS             # don't carry the proxy flags into later runs

# ---- R environment ----
module load stack/2024-06 gcc/12.2.0 r/4.4.0 r-ggplot2/3.4.4-dgi26ie # R + ggplot2
module load cmake udunits gdal geos proj 
module load curl/8.4.0-s6dtj75 # curl
module load openssl/3.1.3-zhfub4o # openssl
module load libpng/1.6.39-fz4tvmr # libpng
module load icu4c/67.1-ke2gg6x # icu4c
module load sqlite # sqlite
module load libtiff/4.5.1-cbrwrtj # libtiff

Rscript -e 'renv::install("sf@1.1-0", type="source")'
Rscript -e 'renv::restore(prompt=FALSE)'

# mkdir -p "$SCRATCH/Rlibs" "$SCRATCH/thesis/renv"
# export R_LIBS_USER=$SCRATCH/thesis/Rlibs
# export RENV_PATHS_ROOT=$SCRATCH/thesis/renv
# Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org")'
# Rscript -e 'renv::restore(prompt = FALSE)'

# ---- Python environment ----
# export PATH="$HOME/.local/bin:$PATH"
# export UV_CACHE_DIR=$SCRATCH/uv-cache
# export UV_PYTHON_INSTALL_DIR=$SCRATCH/uv-python
export UV_NO_CONFIG=1
uv sync --frozen

echo "Setup finished."