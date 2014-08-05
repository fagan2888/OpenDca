/*
Copyright (c) 2014, UT-Battelle, LLC

OpenDca, Version 1.0

This file is part of OpenDca.
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
*/
#ifndef DCA_PARAMETERS_H
#define DCA_PARAMETERS_H
#include "String.h"

namespace OpenDca {

template<typename RealType_>
struct DcaParameters {

	typedef RealType_ RealType;

	/** @class hide_DcaParameters
	  - Orbitals=integer
	  - BathSitesPerSite=integer
	  - DcaOmegaBegin=real Start of real freq. domain
	  - DcaOmegaTotal=integer Total number of real frequencies
	  - DcaOmegaStep=real Step in real freq. domain
	  - DcaLargeKs=integer Number of ``cluster'' points (exclude the bath here)
	  - NumberOfMatsubaras=integer Total number of matsubara freq.
	  - DcaFinePoints=integer Number of points for the coarse-graining. This is a
	    maximum, and the geometry mesh will refine it if needed. Must be a perfect square.
	  - DcaBeta=real Inverse temperature in hopping units
	  - DcaMu=real Chemical potential (ignored used potentialV for now)
	  - DcaDelta=real
	  - DcaSolver=string Either Lanczos or Dmrg
	  - DcaMatsubaraIterations=integer Number of self-consistent iterations done
	    with Matsubara freq.
	  - DcaRealFreqIterations=integer Number of self-consistent iterations done
	    with real freq.
	  - DcaOptions=string Either none or nomufeature
	  - DcaAndersonFitCutoff=real Maximum Matsubara freq. in absolute value to be
	                              considered for the Anderson fitting.
	  - Threads=integer [Optional]
	*/
	template<typename SomeInputType>
	DcaParameters(SomeInputType& io)
	{
		io.readline(orbitals,"Orbitals=");
		io.readline(nofPointsInBathPerClusterPoint,"BathSitesPerSite=");
		io.readline(electronsUp,"TargetElectronsUp");
		io.readline(electronsDown,"TargetElectronsDown");
		io.readline(omegaBegin,"DcaOmegaBegin=");
		io.readline(omegas,"DcaOmegaTotal=");
		io.readline(omegaStep,"DcaOmegaStep=");
		io.readline(largeKs,"DcaLargeKs=");
		io.readline(numberOfMatsubaras,"NumberOfMatsubaras=");
		io.readline(smallKs,"DcaFinePoints=");
		io.readline(beta,"DcaBeta=");
		io.readline(mu,"DcaMu=");
		io.readline(targetDensity,"DcaTargetDensity=");
		io.readline(delta,"DcaDelta=");
		io.readline(dcaSolver,"DcaSolver=");
		io.readline(imagIterations,"DcaMatsubaraIterations=");
		io.readline(realIterations,"DcaRealFreqIterations=");
		io.readline(dcaOptions,"DcaOptions=");
		io.readline(andersonFitCutoff,"DcaAndersonFitCutoff=");

		nthreads = 1;
		try {
			io.readline(nthreads,"Threads=");
		} catch (std::exception& e) {}

		try {
			io.read(potentialV,"potentialV");
		} catch (std::exception& e) {
			std::cerr<<"DcaParameters: No potentialV found\n";
		}
	}

	SizeType imagIterations;
	SizeType realIterations;
	SizeType electronsUp;
	SizeType electronsDown;
	SizeType omegas;
	SizeType largeKs;
	SizeType smallKs;
	SizeType numberOfMatsubaras;
	SizeType nofPointsInBathPerClusterPoint;
	SizeType orbitals;
	SizeType nthreads;
	SizeType andersonFitCutoff;
	RealType omegaBegin;
	RealType omegaStep;
	RealType beta;
	mutable RealType mu;
	RealType targetDensity;
	RealType delta;
	PsimagLite::String dcaSolver;
	PsimagLite::String dcaOptions;
	typename PsimagLite::Vector<RealType>::Type potentialV;
}; // struct DcaParameters

}

#endif

