#!/bin/bash

for dir in pehill_base-*/
do
    cd $dir
    rm -rf 0/uniform
    rm -rf 0/dz.dat
    cp 0/U 0/Uorg
    cp ../sampleDict-org ./system/sampleDict
    sample -time 0
    sample -time 3000.000000
    sample -latestTime
    cp ../sampleDict-nonorg ./system/sampleDict
    sample -time 0
    sample -latestTime
    cd ../
done
