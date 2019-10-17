#!/bin/bash
set -ev

if [[ "${TRAVIS_BRANCH}" = "master" ]] ; then
    export ISR_machine=www.immunespace.org
else
    export ISR_machine=test.immunespace.org
fi

R CMD build .
R CMD check --no-build-vignettes ImmuneSpaceR*tar.gz
