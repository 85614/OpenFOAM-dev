#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

./m.sh

./Allrerun

git status