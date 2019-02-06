#!/bin/bash

echo --- updating repo ---
cd ~/repos/gnm
git pull

echo --- updating scripts ---
find ~/results -maxdepth 1 -type f -delete
cp ~/repos/abmianalytics/birds/cc/* ~/results

echo -- done --
