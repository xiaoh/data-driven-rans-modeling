#!/usr/bin/env bash

# Use: image_montage.sh ncols outfile
#   ncols: number of columns
#   outfile: name of output image. Optional. Default 'out.png'.
# >>image_montage.sh 4 montage.png

# density x: resolution of output in pixels/inch.
# geometry +x+y-: space between images (in pixels) in the x and y directions.
# border x: put a border of width x (pixels) around each image.
# tile ax0: create 'a' columns and as many rows as needed.
# >>montage -density 300 -tile 2x0 -geometry +5+50 -border 10 *.png out.png

ncol=${1?Error: Number of columns not specified.}
outfile=${2:-montage.png}

counter=0
for a in $(ls -1v *.png); do counter=$((counter+1)); convert "$a" -set label "Mode $counter" "$a"; done
montage -density 100 -tile ${ncol}x0 -geometry +5+5- $(ls -1v *.png) $outfile
