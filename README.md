# xgems-jupyter

[xGEMS](https://gemshub.github.io/site/start/gemstandalone/download/#xgems) examples using Jupyter notebooks

## Try xGEMS (demo) in your browser

Using Jupyter Lab Notebooks and Binder. This is just a demo and the work will be lost once you close your browser.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gemshub/xgems-jupyter/main?urlpath=lab/tree/how-to-use-xgems-examples.ipynb)  

Wait until the Jupyter Lab Notebook server starts (~1 min) then double click on any `how-to-...` tutorial notebook.

More information on Jupyter Notebooks: [Jupyter Documentation](https://jupyter.readthedocs.io/en/latest/index.html)

## Run xGEMS locally on your computer

Run Jupyter Notebooks examples/tutorials on jupyter notebook/lab running locally on your computer. First checkout this repository or simply click on code download zip. 

### Install xgems

For this you first need to have [Miniforge](https://github.com/conda-forge/miniforge/releases) or [Miniconda](https://conda.io/miniconda.html) installed (follow the instructions on their websites).

In the Linux/MacOS terminal or Windows miniconda/miniforge command prompt or terminal change to the directory jupyter (that you cloned from this repository) and execute the following commands: 

```
conda env create -f environment.yml
conda activate xgems
``` 
First command will create the `xgems` conda environment congaing all the necessary libraries to run xGEMS and associated codes. The second command will activate the environment and your terminal line should show `xgems`. The third command installs jupyter lab in this environment.

### Start jupyter lab

To start jupyter lab in your terminal/command prompt navigate to the desired folder (preferably where you have the notebook files) and execute:

```
jupyter lab
``` 

## xGEMS examples

This repository contains examples of uses cases and functionality of xGEMS. Additional information about the methods and data used in xGEMS you can find in the [documentation](https://xgems.readthedocs.io/en/latest/) or [repository](https://bitbucket.org/gems4/xgems).

Start with how-to-use-xgems-examples.ipynb or look in the folder how-to-use-xgems-examples. 
