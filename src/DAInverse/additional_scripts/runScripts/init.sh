#!/usr/bin/env bash

# PYTHONPATH
export PYTHONPATH="$GIT/src:$PYTHONPATH"
export PYTHONPATH="$GIT/src/library/randomField/:$PYTHONPATH"

export PYTHONPATH="$GIT/src/DAInverse:$PYTHONPATH"
export PYTHONPATH="$GIT/src/DAInverse/postprocessing/pehill:$PYTHONPATH"

# PATH
export PATH="$GIT/src/DAInverse:$PATH"
export PATH="$GIT/src/DAInverse/postprocessing/pehill:$PATH"
