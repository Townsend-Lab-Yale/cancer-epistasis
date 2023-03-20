#!/bin/bash

echo "Checking for python version..."
python_version=`python --version`
if [ "$python_version" != "Python 3.8.1" ]
then
    echo -e "\r\nWARNING: We are working with python version 3.8.1, "`
    `"it has to be installed before running this setup.\r\n\r\n"`
    `"A good way to have multiple python versions on your system is "`
    `"by installing pyenv (https://github.com/pyenv/pyenv).\r\n\r\n"`
    `"If pyenv is installed, run:\r\n\r\n"`
    `"env PYTHON_CONFIGURE_OPTS=\"--enable-shared\" pyenv install 3.9.5\r\n"`
    `"pyenv local 3.8.1\r\n"
    read -p "If you still want to continue with a different python version, "`
    `"press enter to continue or Ctrl-C to abort setup."
fi
echo $python_version
echo ""


echo "Creating a python virtual environment..."
python -m venv .venv
echo "done."
echo ""


echo "Installing required python packages locally..."
echo ""
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
echo "done installing python packages."
echo ""
echo ""


echo "Creating directories for results"
cd ..
mkdir -p results
echo ""
echo "...done creating directories."
echo ""


echo "Downloading TCGA data..."
echo ""
mkdir -p 'data'
cd 'data'
wget 'https://misc.cidma.us/data/tcga.luad.maf.txt'
echo "...done."
cd ..

echo "Setup ready."
