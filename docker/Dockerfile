# Base image
FROM continuumio/anaconda3:2021.05

# Install Vim
RUN apt-get update -y && apt-get install -y vim

# Add channel
RUN conda config --add channels conda-forge

# package install
RUN conda install -y wheel==0.36.2 twine==3.4.1
RUN conda install -y numpy==1.20.1 pandas==1.2.4