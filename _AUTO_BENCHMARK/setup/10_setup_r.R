# Notes:
#   This file is used for installing the required R packages that this repo
#   will be using. We can use it to install both CRAN and GitHub packages.

CRAN_PKGS <- c(
    "devtools",
    "data.table",
    "haven",
    "here",
    "fixest",
    "glmmML",
    "alpaca",
    "plm"
)

GH_PKGS <- c("sgaure/lfe@v2.8-4")


#' Return whether the package with the given name is installed
#'
#' @param pkg The name of the package to check if it is installed
is_installed <- function(pkg) {
    return(pkg %in% installed.packages()[, 1])
}

#' This is a quick function to perform package installs for both the
#' CRAN package as well as the GitHub packages via devtools
#'
#' @param pkg_vec A vector of packages to install
#' @param cran If TRUE, install from CRAN, else install from github
#' @param force If TRUE, install even if the package is already installed
install_pkgs <- function(pkg_vec, cran = TRUE, force = FALSE) {
    if (cran) {
        for (pkg in pkg_vec) {
            if (force | !is_installed(pkg)) {
                install.packages(
                    pkg,
                    repos = "https://cran.r-project.org",
                    Ncpus = 3
                )
            }
        }
    } else {
        for (pkg in pkg_vec) {
            if (force | !is_installed(pkg)) {
                devtools::install_github(pkg)
            }
        }
    }
}


install <- function() {
    install_pkgs(CRAN_PKGS)
    install_pkgs(GH_PKGS, cran = FALSE)
}

install()