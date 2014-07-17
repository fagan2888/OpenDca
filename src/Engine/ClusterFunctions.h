#ifndef CLUSTERFUNCTIONS_H
#define CLUSTERFUNCTIONS_H
#include "DcaSolverLanczos.h"
#include "DcaSolverDmrg.h"
#include "../../../dmrgpp/src/Engine/ProgramGlobals.h"
#include "Geometry/Geometry.h"
#include "Complex.h"

namespace OpenDca {

template<typename DcaToDmrgType>
class ClusterFunctions {

	typedef typename DcaToDmrgType::ParametersType ParamtersType;
	typedef typename DcaToDmrgType::InputNgType InputNgType;
	typedef typename DcaToDmrgType::RealType RealType_;
	typedef PsimagLite::Geometry<RealType_,
	                             DcaToDmrgType,
	                             Dmrg::ProgramGlobals> VaryingGeometryType;
	typedef DcaSolverLanczos<DcaToDmrgType, VaryingGeometryType> DcaSolverLanczosType;
	typedef DcaSolverDmrg<DcaToDmrgType, VaryingGeometryType>  DcaSolverDmrgType;
	typedef typename DcaSolverDmrgType::DcaSolverBaseType DcaSolverBaseType;
	typedef typename DcaSolverBaseType::VectorSizeType VectorSizeType;
	typedef typename DcaSolverBaseType::PlotParamsType PlotParamsType;
	typedef std::complex<RealType_> ComplexType;

public:

	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef RealType_ RealType;

	ClusterFunctions(DcaToDmrgType* myInputPtr,
	                 const ParamtersType& params,
	                 typename InputNgType::Readable& io)
	    : myInputPtr_(myInputPtr), params_(params), io_(io)
	{}

	RealType operator()(RealType mu) const
	{

		return 0.0;
	}

	void findGf(MatrixType& gfCluster)
	{
		DcaToDmrgType& myInput = *myInputPtr_;

		VaryingGeometryType geometry2(myInput,false,params_.smallKs);

		DcaSolverBaseType* solver = allocateSolverPtr(myInput,geometry2);

		RealType omegaEnd = params_.omegas * params_.omegaStep + params_.omegaBegin;

		SizeType matsubaras = params_.numberOfMatsubaras;
		if (params_.dcaOptions.find("lanczosreal") != PsimagLite::String::npos)
			matsubaras = 0;

		PlotParamsType plotParams(params_.omegaBegin,
		                          omegaEnd,
		                          params_.omegaStep,
		                          params_.delta,
		                          params_.beta,
		                          matsubaras);

		SizeType Nc = params_.largeKs;
		VectorSizeType sites(2);
		for (SizeType i = 0; i < Nc; ++i) {
			for (SizeType j = i; j < Nc; ++j) {
				sites[0] = i;
				sites[1] = j; // dca indexing
				solver->solve(gfCluster,sites,plotParams);
			}
		}

		for (SizeType orb = 0; orb< params_.orbitals; ++orb)
			for (SizeType i = 0; i < Nc; ++i)
				for (SizeType j = 0; j < i; ++j)
					for (SizeType x = 0; x < gfCluster.n_row(); ++x)
						gfCluster(x,i+j*Nc+orb*Nc*Nc) =
						          gfCluster(x,j+i*Nc+orb*Nc*Nc);

		deAllocateSolverPtr(&solver);
	}

	void sweepParticleSectors()
	{
		DcaToDmrgType& myInput = *myInputPtr_;

		RealType Eg = 1e6;
		SizeType iMin = 0;
		for (SizeType i = 0; i < myInput.muFeatureSize(); ++i) {
			myInput.muFeatureSet(i);
			SizeType total = myInput.electrons(DcaToDmrgType::SPIN_UP) +
			        myInput.electrons(DcaToDmrgType::SPIN_DOWN);
			if (total == 0) continue;
			std::cout<<"Trying with "<<myInput.electrons(DcaToDmrgType::SPIN_UP);
			std::cout<<" electrons up and ";
			std::cout<<myInput.electrons(DcaToDmrgType::SPIN_DOWN)<<" electrons down\n";
			VaryingGeometryType geometry2(myInput,false,params_.smallKs);
			DcaSolverBaseType* solver = allocateSolverPtr(myInput,geometry2);
			RealType tmp = solver->findLowestEnergy();
			deAllocateSolverPtr(&solver);
			if (i > 0 && tmp > Eg) continue;
			iMin = i;
			Eg = tmp;
		}

		if (myInput.muFeatureSize() > 0) {
			myInput.muFeatureSet(iMin);
			std::cout<<"EffectiveHamiltonian: found lowest energy "<<Eg;
			std::cout<<" in mu sector "<<iMin;
			std::cout<<" electrons up "<<myInput.electrons(DcaToDmrgType::SPIN_UP);
			std::cout<<" electrons down "<<myInput.electrons(DcaToDmrgType::SPIN_DOWN);
			std::cout<<"\n";
		}
	}

private:

	DcaSolverBaseType* allocateSolverPtr(DcaToDmrgType&myInput,
	                                     const VaryingGeometryType& geometry2)
	{
		PsimagLite::String dcaSolver = params_.dcaSolver;
		DcaSolverBaseType* solver = 0;
		if (dcaSolver == "Dmrg") {
			solver = new DcaSolverDmrgType(myInput,geometry2,io_);
		} else if (dcaSolver == "Lanczos") {
			solver = new DcaSolverLanczosType(myInput,geometry2,io_);
		} else {
			PsimagLite::String str("makeHubbardParams(): Unknown solver ");
			str += dcaSolver;
			throw PsimagLite::RuntimeError(str);
		}

		return solver;
	}

	void deAllocateSolverPtr(DcaSolverBaseType** solver) const
	{
		DcaSolverBaseType* solver2 = *solver;
		delete solver2;
		solver2 = 0;
	}

	DcaToDmrgType* myInputPtr_;
	const ParamtersType& params_;
	typename InputNgType::Readable& io_;
}; // ClusterFunctions

} // namespace OpenDca

#endif // CLUSTERFUNCTIONS_H

