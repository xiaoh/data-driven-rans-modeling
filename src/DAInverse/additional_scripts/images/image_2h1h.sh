#!/bin/bash

for a in *.png
do
 convert -crop 740x888 "$a" "$a" 
 convert -trim "$a" "$a"
 convert -border 2x2 -bordercolor "#FFFFFF" "$a" "$a"
done

