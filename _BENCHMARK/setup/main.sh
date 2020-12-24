#!/bin/sh

JULIA=julia
R=Rscript
STATA=/Applications/Stata/StataSE.app/Contents/MacOS/StataSE

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo ${machine}

if [[ "${machine}" == "Mac" ]]; then
    sh ./00_macos.sh
fi;

${R} --vanilla ./10_setup_r.R
${JULIA} ./20_setup_julia.jl
${STATA} -b do ./30_setup_stata.do