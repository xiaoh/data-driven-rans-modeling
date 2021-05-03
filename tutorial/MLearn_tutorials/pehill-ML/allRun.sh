#!/bin/bash

# First, you need go to MLearn/RANSTemp to generate RANS case, this only need to be generate once.
# (1) Go to the src/MLearn/RANSTemp folder, change the paraInput.in. Specifying which mesh, which turbulence model, 
#     which openFoam solver, final time step, and a group of Reynolds numbers that you want to simulate. 
cd ../../../src/MLearn/RANSTemp/peHill/
# (2) Run generateRANSCases.py
python generateRANSCases.py
# (3) Run generateFeatureFields.py
python  generateFeatureFields.py
# Only after done all above once (got converged RANS cases), then you can run the machine learning.
cd -
mLearnMain.py mainInput.in
