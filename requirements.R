# R packages used by the CNS_identification pipeline.
# Run from inside the conda environment created by environment.yml:
#
#   conda activate cns-nfc
#   Rscript requirements.R
#
# Bioconductor packages (ComplexHeatmap) are installed through BiocManager.

cran_packages <- c(
    "docopt",
    "vroom",
    "tidyverse",
    "circlize",
    "splitstackshape",
    "plyr",
    "BiocManager"
)

to_install <- cran_packages[!cran_packages %in% rownames(installed.packages())]
if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
}

# Bioconductor packages
bioc_packages <- c("ComplexHeatmap")
to_install_bioc <- bioc_packages[!bioc_packages %in% rownames(installed.packages())]
if (length(to_install_bioc) > 0) {
    BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)
}

# Sankey rendering uses webshot + PhantomJS. PhantomJS must be installed once,
# outside of a batch context, with:
#
#   Rscript -e "webshot::install_phantomjs()"
#
# It is intentionally not invoked here so that downstream Snakemake rules do
# not depend on a network call at runtime.
if (!"webshot" %in% rownames(installed.packages())) {
    install.packages("webshot", repos = "https://cloud.r-project.org")
}

cat("All required R packages are installed.\n")
