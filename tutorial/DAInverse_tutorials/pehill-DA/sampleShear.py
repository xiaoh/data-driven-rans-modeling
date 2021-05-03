#!/bin/bash

for dir in pehill-*/
do
    cd $dir
    wallShearStress
	cp ../sampleDict.surface system/sampleDict
    sample
    cd ../
done
