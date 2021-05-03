
for d in ./pehill_base-* ; do
    cp pehill_base/system/sampleDict $d/system/sampleDict
    cd $d
    cp system/controlDict.org system/controlDict
    sample > log.sample3
    cd ../
done
