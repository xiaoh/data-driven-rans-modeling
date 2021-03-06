#!/bin/bash

# First, you need go to MLearn/RANSTemp to generate RANS case, this only need to be generate once. So, I didn't put it
# into this run.sh script. Please follow the guide:
# (1) Go to the src/MLearn/RANSTemp folder, change the paraInput.in. Specifying which mesh, which turbulence model, 
#     which openFoam solver, final time step, and a group of Reynolds numbers that you want to simulate. 
# (2) Run generateRANSCases.py
# (3) Run generateFeatureFields.py

# Only after done all above once (got converged RANS cases), then you can run the machine learning.

mLearnMain.py mainInput.in
