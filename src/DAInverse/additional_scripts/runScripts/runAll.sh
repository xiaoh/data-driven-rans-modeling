#!/usr/bin/env bash

. ~/.openfoam
. ./init.sh


################################################################################
cd SE5
echo 'Running SE5 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 1/8 done'


################################################################################
cd GC5
echo 'Running GC5 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 2/8 done'


################################################################################
cd SE10
echo 'Running SE10 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 3/8 done'


################################################################################
cd GC10
echo 'Running GC10 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 4/8 done'

################################################################################
cd GC20
echo 'Running GC20 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 5/8 done'

################################################################################
cd SE20
echo 'Running SE20 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 6/8 done'

################################################################################
cd GC3
echo 'Running GC3 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 7/8 done'

################################################################################
cd SE3
echo 'Running SE3 ...'
./run.sh
echo 'Plotting ...'
UPlotWithMeanDNS.py
echo 'Deleting intermediate values ...'
./delint.sh
cd ../
echo 'case 8/8 done'

echo 'DONE.'
