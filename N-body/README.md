# N-body

This is an open-source code to compute the non-linear evolution of the cosmological matter density contrast. It is based on the Particle Mesh (PM) technique. It is written in C. The code has been parallelized (Mondal et al. 2015) for shared-memory machines using Open Multi-Processing (OpenMP).

Read the user's guide 'nbody_doc.pdf' for a detailed description.

_____________________________________

Requirments for Running the program on debian based system.
1. FFTW -3 installation command
sudo apt-get install -y fftw3

2. Make and gcc
sudo apt-get install -y make gcc

3. Git to downlnload the code
sudo apt-get git

Single line of code for compleat setup:
sudo apt update && sudo apt install git make gcc fftw3 


Requirments for Running the program on Arch based system.
1. FFTW -3 installation command

2. Make and gcc
sudo pacman -S make gcc

3. Git to downlnload the code
sudo pacman -S git

Single line of code for compleat setup:
sudo pacman -Syu git make gcc

_____________________________________


installing fftw3
curl https://www.fftw.org/fftw-3.3.10.tar.gz --output fftw-3.3.10.tar.gz
tar -zxvf fftw-3.3.10.tar.gz
cd fftw-3.3.10

./configure --enable-float --enable-threads --enable-openmp
make
make install

_____________________________________


Download the code by cloning the git repository using

$ git clone https://github.com/rajeshmondal18/N-body.git

You need to install FFTW-3.x.x with the following flags: '--enable-float',  '--enable-threads' and '--enable-openmp' to compile this set of codes.
Look at the installation instruction http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix
_____________________________________


 

---
Use the makefile for compilation in the following manner:

$ make nbody_comp

It will create the executable 'nbody_comp'.

To Run:

$ ./nbody_comp
_____________________________________
Please acknowledge these papers 
1. Bharadwaj and Srikant 2004 (http://adsabs.harvard.edu/abs/2004JApA...25...67B), if you are using the serial version of the code
2. Mondal et al. 2015 (http://adsabs.harvard.edu/abs/2015MNRAS.449L..41M), if you are using the current parallelized version of the code.