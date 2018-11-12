#!/bin/bash
# Create a fooof conda environment to run ipython notebooks

# Download MNE conda envirnoment yaml:
#curl -O https://raw.githubusercontent.com/mne-tools/mne-python/master/environment.yml
#mv environment.yml ${root_dir}/PRJ_Stroop/scripts/mne_env.yml

# Create conda environment with all the MNE dependencies
conda env create -f mne_env.yml jupyter
#conda create -n fooofenv python=3.6 jupyter
source activate mne
# If you are on macOS, you need to manually update PyQt5. This step is not needed on Linux, and even breaks things on Windows.
pip install --upgrade "pyqt5>=5.10"
#NOT needed with MNE yaml: pip install numpy scipy matplotlib pytest
pip install fooof

## To make sure ipython notebooks can access the py3.6 kernel from this environment:
#NOT needed with MNE yaml: conda install ipython ipykernel 
python -m ipykernel install --user --name mne --display-name "Py36 (mne)"

# mne version 0.16.2 didn't have the pymatloader functions to access fieldtrip, so upgrading:
pip install --upgrade --no-deps git+https://github.com/mne-tools/mne-python.git
# installed v0.17.dev0
