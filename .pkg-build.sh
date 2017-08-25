#!/bin/bash
set -ev
R CMD build .
R CMD check --no-build-vignettes ImmuneSpaceR*tar.gz
