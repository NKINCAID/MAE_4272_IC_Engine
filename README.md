Guide for Installing Cantera with a Conda environment
=========


Conda
---------
1. Install [Anaconda](https://www.anaconda.com/download/) on your machine (make sure to select Python 3!)

2. Make sure the installation went smoothly:
```bash
conda --version
```
If you get an error or don't see the version output, it's possible that your current bash doesn't know where conda is! This can be fixed by adding the anaconda or miniconda bin location to your systems PATH variable to your *.bash_profile* or *.bashrc* 
```bash
# add path to conda
export PATH="Your_path_to_anaconda/anaconda3/bin:$PATH"
```
You may need to `source .bash_profile` or just restart your bash.

Creating a Virtual Envrionment and Installing Cantera
---------

1. Use Conda to create a virtual environment. I'll name the environment `ct-env`, although you can name it whatever you'd like!
```bash
conda create --name ct-env --channel cantera cantera ipython matplotlib jupyter scipy
```
2. If `ct-env` is activated, you'll see `(ct-env)` to the left of the terminal prompt. If not, you can activate it by
```bash
conda activate ct-env
```

3. Onece `ct-env` is activated, test out Cantera in ipython
```bash
ipython
```
```python
import cantera as ct
print(ct.avogadro)
```

For more information on installing Cantera checkout the [Cantera website](https://cantera.org/install/conda-install.html#sec-install-conda) 