#!/bin/bash

echo --- updating repo ---
cd ~/repos/gnm
git pull

echo --- updating scripts ---
find ~/bam -maxdepth 1 -type f -delete
cp ~/repos/GNM/graham/* ~/bam

echo -- done --
