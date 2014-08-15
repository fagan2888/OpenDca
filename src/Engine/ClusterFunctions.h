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

	typedef typename DcaToDmrgType::ParametersType ParametersType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef RealType_ RealType;

	ClusterFunctions(const ParametersType& params,
	                 typename InputNgType::Readable& io)
	    : params_(params), io_(io), myInputPtr_(0)
	{}

	void build(DcaToDmrgType* myInputPtr)
	{
		myInputPtr_ = myInputPtr;
	}

	// <gs|n|gs>
	RealType operator()(RealType mu) const
	{
		params_.mu = mu;

		updatePotentialV();

		RealType n = sweepParticleSectors();

		return n - params_.targetDensity;
	}

	void findGf(MatrixType& gfCluster)
	{
		if (!myInputPtr_) {
			PsimagLite::String str("findGf: must call build first\n");
			throw PsimagLite::RuntimeError("ClusterFunctions::" + str);
		}

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

	RealType sweepParticleSectors(bool last = false) const
	{
		if (!myInputPtr_) {
			PsimagLite::String str("sweepParticleSectors: must call build first\n");
			throw PsimagLite::RuntimeError("ClusterFunctions::" + str);
		}

		DcaToDmrgType& myInput = *myInputPtr_;

		RealType Eg = 1e6;
		SizeType iMin = 0;
		RealType densityGs = 0;
		for (SizeType i = 0; i < myInput.particleSectors(); ++i) {
			myInput.particleSectorsSet(i);
			SizeType total = myInput.electrons(SPIN_UP) + myInput.electrons(SPIN_DOWN);
			if (total == 0) continue;
			std::cout<<"Trying with "<<myInput.electrons(SPIN_UP);
			std::cout<<" electrons up and ";
			std::cout<<myInput.electrons(SPIN_DOWN)<<" electrons down\n";
			VaryingGeometryType geometry2(myInput,false,params_.smallKs);
			DcaSolverBaseType* solver = allocateSolverPtr(myInput,geometry2);
			RealType density = solver->density(geometry2.numberOfSites(),params_.largeKs);
			RealType tmp = solver->findLowestEnergy();
			deAllocateSolverPtr(&solver);
			if (i > 0 && tmp > Eg) continue;
			iMin = i;
			Eg = tmp;
			densityGs = density;
		}

		if (myInput.particleSectors() > 0) {
			myInput.particleSectorsSet(iMin);
			std::cout<<"EffectiveHamiltonian: found lowest energy "<<Eg;
			if (last) std::cout<<" BEFORE GF ";
			std::cout<<" in mu sector "<<iMin;
			std::cout<<" electrons up "<<myInput.electrons(SPIN_UP);
			std::cout<<" electrons down "<<myInput.electrons(SPIN_DOWN);
			std::cout<<"mu = "<<params_.mu<<"\n";
		}

		return densityGs;
	}

	void updatePotentialV() const
	{
		assert(myInputPtr_);
		myInputPtr_->updatePotentialV();
	}

	SizeType electrons(SpinEnum upOrDown) const
	{
		assert(myInputPtr_);
		return myInputPtr_->electrons(upOrDown);
	}

private:

	DcaSolverBaseType* allocateSolverPtr(DcaToDmrgType&myInput,
	                                     const VaryingGeometryType& geometry2) const
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

	const ParametersType& params_;
	typename InputNgType::Readable& io_;
	DcaToDmrgType* myInputPtr_;
}; // ClusterFunctions

} // namespace OpenDca

#endif // CLUSTERFUNCTIONS_H

