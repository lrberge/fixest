#!/bin/bash

if [ "$EUID" -ne 0 ]; then
  echo "Please run as root"
  exit
fi

JULIA_VERSION="1.5.3"

# Setup apt for HTTPS traffic
apt-get update -y
apt-get install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common

# Setup R
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
apt-get update -y
apt install -y r-base
chmod 777 /usr/local/lib/R/site-library

# Setup Julia
wget "https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
tar xzvf "julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
sudo mv "julia-${JULIA_VERSION}" /opt
sudo ln -s "/opt/julia-${JULIA_VERSION}/bin/julia" /usr/local/bin
rm "julia-${JULIA_VERSION}-linux-x86_64.tar.gz"