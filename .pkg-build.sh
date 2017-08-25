#!/bin/bash
set -ev
travis_wait 30 R CMD build .
travis_wait 30 R CMD check --no-build-vignettes ImmuneSpaceR*tar.gz
