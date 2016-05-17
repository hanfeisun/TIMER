

[![Build Status](https://travis-ci.org/hanfeisun/TIMER.svg?branch=master)](https://travis-ci.org/hanfeisun/TIMER)



STEP 1: Installing Miniconda3 if you don't have conda installed
---------------------------------------------------------------
If you are using Linux or sshed into a Linux system:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
If you are using Mac locally:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

Follow the commands listed on screen, press the enter key to scroll down.
**Make sure to answer yes when asked if you want to prepend Miniconda3 to PATH**.

Close your terminal, open a new one and you should now have Conda working! Test by entering:
```
conda update conda
```

you may need to add channels for convenience
```
conda config --add channels r
```


STEP 2: get TIMER source code from github, and create TIMER env
---------------------------------------------------------------
```
git clone https://github.com/hanfeisun/TIMER
conda env create -f TIMER/environment.yml
```

activate the TIMER conda environment
```
source activate timer
```

If you want to deactivate the TIMER enviroment, type
```
source deactivate
```

