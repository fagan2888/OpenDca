#ifndef DCA_PARAMETERS_H
#define DCA_PARAMETERS_H
#include "String.h"

namespace OpenDca {

template<typename RealType_>
struct DcaParameters {

	typedef RealType_ RealType;

	template<typename SomeInputType>
	DcaParameters(SomeInputType& io)
	{
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
		io.readline(iterations,"DcaIterations=");
	}

	SizeType iterations;
	SizeType electronsUp;
	SizeType electronsDown;
	SizeType omegas;
	SizeType largeKs;
	SizeType smallKs;
	SizeType numberOfMatsubaras;
	SizeType nofPointsInBathPerClusterPoint;
	RealType omegaBegin;
	RealType omegaStep;
	RealType beta;
	RealType mu;
	RealType delta;
	PsimagLite::String dcaSolver;
}; // struct DcaParameters

}

#endif

