# Quick Start
 
# Disclaimer and Licensing
 
OpenDca is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
OpenDca is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with OpenDca. If not, see <http://www.gnu.org/licenses/>.
The full software license for OpenDca version 1.0.0 
can be found in
file COPYING. 

# Please cite this work

OpenDca is a free and open source implementation of the 
DCA algorithm for models of strongly correlated electrons. 
The full software license for OpenDca version 1.0.0 
can be found in
file COPYING. 
You are welcomed to use it and publish data 
obtained with OpenDca. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.


# Hash of the latest commit 

Hash of the latest commit is also posted at
https://web.ornl.gov/~gz1/hashes.html

# Building and Running OpenDca

## Required Software

* GNU C++
* The LAPACK and BLAS libraries
* The GSL library
* PsimagLite (see below)
* DMRG++ (see below)
* Lanczos++ (see below)

## Optional Software

* make or gmake (only needed to use the Makefile)
* perl (may be needed to run some auxiliary script)
* doxygen (to build the manual) 

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
the LAPACK and BLAS libraries, the gsl library, make, and git 
if you don't have them.

2. Issue

    cd someDirectory/

    git clone https://github.com/g1257/PsimagLite.git

    git clone https://github.com/g1257/dmrgpp.git

    git clone https://github.com/g1257/LanczosPlusPlus.git

    git clone https://github.com/g1257/OpenDca.git

3. Compile PsimagLite

    cd PsimagLite/lib/

    make -f Makefile.sample

    cd ../../

4. Now issue

    cd OpenDca/src/

    make

5. You can run it with

    ./openDca -f TestSuite/inputs/input1.inp &> output

6. Post-process with

    perl ../scripts/multiExtract.pl output 0 3 100 200 -10 0.01

Run 

    perl ../scripts/multiExtract.pl

to see the meaning of the arguments.

# Manual

To compile the manual please follow these steps

    cd OpenDca

    doxygen Doxygen

    cd latex/

    make

and you'll get the file `refman.pdf`
