 
# Licensing
The full software license for OpenDca version 1.0.0 
can be found in
file LICENSE. 
OpenDca is a free and open source implementation of the 
DCA algorithm for models of strongly correlated electrons. 
You are welcomed to use it and publish data 
obtained with OpenDca. If you do, please cite this
work .

-------------------------------------------------------------------------------

# Hash of the latest commit is also posted at

https://web.ornl.gov/~gz1/hashes.html

-------------------------------------------------------------------------------

# How To Cite This Work

TBW.

-------------------------------------------------------------------------------

# DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

-------------------------------------------------------------------------------

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

## Quick Start

1. Use your distribution repository tool to install gcc with support for C++,
the LAPACK and BLAS libraries, the gsl library, make, and git 
if you don't have them.

2. Issue

```
cd someDirectory/
git clone https://github.com/g1257/PsimagLite.git
git clone https://github.com/g1257/dmrgpp.git
git clone https://github.com/g1257/LanczosPlusPlus.git
git clone https://github.com/g1257/OpenDca.git
```

3. Compile PsimagLite

```
cd PsimagLite/lib/
make -f Makefile.sample
cd ../../
```

4. Now issue

```
cd OpenDca/src/
make
```

5. You can run it with

```
./openDca -f TestSuite/inputs/input1.inp &> output
```

6. Post-process with

```
perl ../scripts/multiExtract.pl output 0 3 100 200 -10 0.01
```

Run 

```
perl ../scripts/multiExtract.pl
```

to see the meaning of the arguments.


