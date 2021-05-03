#!/bin/bash

for a in *.png; do convert -trim "$a" "$a"; convert -border 2x2 -bordercolor "#FFFFFF" "$a" "$a"; done

