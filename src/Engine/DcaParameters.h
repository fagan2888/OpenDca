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
	  - DcaNumberOfMatsubaras=integer Total number of matsubara freq.
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
		io.readline(numberOfMatsubaras,"DcaNumberOfMatsubaras=");
		io.readline(smallKs,"DcaFinePoints=");
		io.readline(beta,"DcaBeta=");
		io.readline(mu,"DcaMu=");
		io.readline(delta,"DcaDelta=");
		io.readline(dcaSolver,"DcaSolver=");
		io.readline(imagIterations,"DcaMatsubaraIterations=");
		io.readline(realIterations,"DcaRealFreqIterations=");
		io.readline(dcaOptions,"DcaOptions=");

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
	RealType omegaBegin;
	RealType omegaStep;
	RealType beta;
	RealType mu;
	RealType delta;
	PsimagLite::String dcaSolver;
	PsimagLite::String dcaOptions;
	typename PsimagLite::Vector<RealType>::Type potentialV;
}; // struct DcaParameters

}

#endif

