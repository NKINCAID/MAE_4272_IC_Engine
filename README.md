Guide for Installing Cantera with a Conda Environment
=========

Getting this repository
---------
To start with, clone this repository to your local machine. To do this, in your terminal enter:
```bash
git clone https://github.com/NKINCAID/MAE_4272_IC_Engine.git
```
or download the folder from [Github](https://github.com/NKINCAID/MAE_4272_IC_Engine).

Conda
---------
1. Install [Anaconda](https://www.anaconda.com/download/) on your machine (make sure to select Python 3!)

2. Make sure the installation went smoothly:
```bash
conda --version
```

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

3. Onece `ct-env` is activated, test out Cantera in ipython. In your terminal enter:
```bash
ipython
```
Once ipython is running, you can use it similar to the command prompt in MATLAB. Test out Cantera with:
```python
import cantera as ct
print(ct.avogadro)
```

For more information on installing Cantera checkout the [Cantera website](https://cantera.org/install/conda-install.html#sec-install-conda) 

Running IC Engine
---------
Run the IC engine code in ipython using:
```python
run ic_engine
```
or in your terminal with:
```bash
python ic_engine.py
```