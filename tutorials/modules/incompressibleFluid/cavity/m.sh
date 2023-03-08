#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

wmake -j /home/zhzhq/OpenFOAM/OpenFOAM-dev/applications/solvers/modules/incompressibleFluid
