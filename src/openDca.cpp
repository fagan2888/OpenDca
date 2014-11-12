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

	dcaLoop.main(PsimagLite::FREQ_MATSUBARA,params.imagIterations);

	dcaLoop.main(PsimagLite::FREQ_REAL,params.realIterations);
}

