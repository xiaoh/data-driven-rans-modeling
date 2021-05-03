# Run this script for each new terminal

export PATH=$PATH:$(pwd)/src/DAInverse
export PATH=$PATH:$(pwd)/src/DAInverse/tools
export PATH=$PATH:$(pwd)/src/DAInverse/postprocessing/pehill
export PATH=$PATH:$(pwd)/src/DAInverse/postprocessing/wingBody
export PATH=$PATH:$(pwd)/src/library/randomField
export PATH=$PATH:$(pwd)/src/randomMatrix/solver
export PATH=$PATH:$(pwd)/src/randomMatrix/preprocess
export PATH=$PATH:$(pwd)/src/randomMatrix/postprocess/common
export PATH=$PATH:$(pwd)/src/randomMatrix/postprocess/forSpecialPhy/pehill
export PATH=$PATH:$(pwd)/src/utilities/postprocessing/pehill
export PATH=$PATH:$(pwd)/src/utilities/postprocessing/channel
export PATH=$PATH:$(pwd)/src/utilities/postprocessing/duct
export PATH=$PATH:$(pwd)/src/utilities/tools
export PATH=$PATH:$(pwd)/src/utilities/postprocessing/diagnosis
export PATH=$PATH:$(pwd)/sandbox/BaycentricMapping-python
export PATH=$PATH:$(pwd)/sandbox/pseduObservationGeneration/pythonScript
export PATH=$PATH:$(pwd)/sandbox/generateCasesFromOmega/pehill


export PYTHONPATH=${PYTHONPATH}:$(pwd)/src
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/DAInverse
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/DAInverse/postprocessing/pehill
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/library/randomField



export PYTHONPATH=${PYTHONPATH}:$(pwd)/sandbox/MLearn/utility
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/randomMatrix/postprocess/common
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/randomMatrix/postprocess/forSpecialPhy
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/utilities/postprocessing/pehill


# initial path for MLearn
export PATH=${PATH}:$(pwd)/src/MLearn/src
export PATH=${PATH}:$(pwd)/src/MLearn/utility/peHill
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/MLearn/src
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/MLearn/src/postprocess
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/MLearn/src/preprocess
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/MLearn/src/postprocess/peHill
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src/MLearn/utility

# initial path for mcmc utility
export PYTHONPATH=${PYTHONPATH}:$(pwd)/sandbox/mcmc/mcmc_python


# 

# echo "export PATH='$PATH':$(pwd)/src/python"  >> ~/.bashrc
