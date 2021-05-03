#!/bin/bash

convert -trim "$1" "$1"
convert -border 2x2 -bordercolor "#FFFFFF" "$1" "$1"

