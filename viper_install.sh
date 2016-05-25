
#install miniconda
if [ "$(uname)" == "Darwin" ]; then
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
fi
conda update conda -y

#reload bash settings
source ~/.bashrc

#clone viper source code
git clone git@bitbucket.org:cfce/viper.git

#create conda viper env
conda env create -f viper/environment.yml

#to use the env
source activate viper

#to leave the env
source deactivate

