#!/usr/bin/env bash

START=$(date +%s.%N)
mfuMain.py MainInput.in &> log.enkf
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF  > log.time

# to monitor use one of these in command line:
    # watch -n 5 "ps -ef | grep mfuMain.py"
    # tail log.enkf
