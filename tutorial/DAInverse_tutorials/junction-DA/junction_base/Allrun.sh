#!/bin/sh
rm -r processor*
decomposePar > log.decomposePar
mpirun -np 8 simpleFoam -parallel > log.simple &
# ----------------------------------------------------------------- end-of-file
