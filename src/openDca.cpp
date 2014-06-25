#include "DcaLoop.h"
#include "Concurrency.h"
#include "InputCheck.h"
#include "Provenance.h"
#include "DcaLoopGlobals.h"
#include "FreqEnum.h"
#include "DispersionSimple.h"
#include "Geometry/Geometry.h"

int main(int argc, char *argv[])
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef PsimagLite::Geometry<ComplexType,
	                             InputNgType::Readable,
	                             Dmrg::ProgramGlobals> GeometryType;
	typedef OpenDca::DispersionSimple<RealType, GeometryType> DispersionType;
	typedef OpenDca::DcaLoop<DispersionType,InputNgType> DcaLoopType;
	typedef DispersionType::ParametersType ParametersType;
	typedef PsimagLite::Concurrency ConcurrencyType;

	ConcurrencyType concurrency(&argc,&argv,1);
	Dmrg::InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	PsimagLite::String insitu("");
	while ((opt = getopt(argc, argv,"f:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="") {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<OpenDca::DcaLoopGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersType params(io);
#ifndef USE_PTHREADS
        inputCheck.checkForThreads(params.nthreads);
#endif

        ConcurrencyType::npthreads = params.nthreads;

	DcaLoopType dcaLoop(params,io);

	dcaLoop.main(OpenDca::FREQ_MATSUBARA,params.imagIterations);

	dcaLoop.main(OpenDca::FREQ_REAL,params.realIterations);
}

