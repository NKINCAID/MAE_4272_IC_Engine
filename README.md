Guide for Installing Cantera with a Conda Environment
=========

The lab notebook can be run either locally on your machine, or on a cloud service such as [Google Colab](https://colab.research.google.com).

* To run on on a cloud server, simply download the notebook from here or Canvas and upload it to Google Colab (or similar).

* If you're interested in running the code locally, I'd recommend using Conda. To do so, follow the instructions below. 


Instructions for running locally
---------
1. To start with, clone this repository to your local machine. To do this, in your terminal enter:
```bash
git clone https://github.com/NKINCAID/MAE_4272_IC_Engine.git
```
or download the folder from [Github](https://github.com/NKINCAID/MAE_4272_IC_Engine) (click on green <>code button and then Download ZIP).

2. Install [Anaconda](https://www.anaconda.com/download/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) on your machine (would recommend Miniconda as it's a faster download/install)

3. Make sure the installation went smoothly by typing the following into your terminal:
```bash
conda --version
```

4. Use Conda to create a virtual environment. I'll name the environment `ct-env`, although you can name it whatever you'd like!
```bash
conda create --name ct-env --channel cantera cantera ipython matplotlib jupyter scipy
```
5. If `ct-env` is activated, you'll see `(ct-env)` to the left of the terminal prompt. If not, you can activate it by
```bash
conda activate ct-env
```

(Alternate route) You can also install via `pip install` to your local python version without installing a virtual environment
```bash
pip install cantera
```


For more information on installing Cantera checkout the [Cantera website](https://cantera.org/install/conda-install.html#sec-install-conda) 
