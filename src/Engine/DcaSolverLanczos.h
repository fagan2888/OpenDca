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
#ifndef LANCZOS_DCA_SOLVER_H
#define LANCZOS_DCA_SOLVER_H
#include "ContinuedFractionCollection.h"
#include "../../LanczosPlusPlus/src/Engine/ProgramGlobals.h"
#include "../../LanczosPlusPlus/src/Engine/DefaultSymmetry.h"
#include "../../LanczosPlusPlus/src/Engine/InternalProductStored.h"
#include "../../LanczosPlusPlus/src/Engine/Engine.h"
#include "../../LanczosPlusPlus/src/Engine/ModelBase.h"
#include "../../LanczosPlusPlus/src/Engine/BasisBase.h"
#include "../../LanczosPlusPlus/src/Engine/ModelSelector.h"
#include "DcaSolverBase.h"
#include "../../PsimagLite/src/FreqEnum.h"

namespace OpenDca {

template<typename DcaToDmrgType,typename VaryingGeometryType>
class DcaSolverLanczos : public DcaSolverBase<DcaToDmrgType, VaryingGeometryType> {

	typedef typename DcaToDmrgType::RealType RealType;
	typedef typename DcaToDmrgType::InputNgType InputNgType;
	typedef LanczosPlusPlus::ModelBase<RealType,
	                                   VaryingGeometryType,
	                                   DcaToDmrgType> ModelType;
	typedef LanczosPlusPlus::BasisBase<VaryingGeometryType> BasisType;
	typedef LanczosPlusPlus::DefaultSymmetry<VaryingGeometryType,
	                                         BasisType> SpecialSymmetryType;
	typedef LanczosPlusPlus::Engine<ModelType,
	                                LanczosPlusPlus::InternalProductStored,
	                                SpecialSymmetryType> EngineType;
	typedef std::pair<SizeType,SizeType> PairType;

public:

	typedef DcaSolverBase<DcaToDmrgType, VaryingGeometryType> DcaSolverBaseType;
	typedef typename DcaSolverBaseType::VectorSizeType VectorSizeType;
	typedef typename DcaSolverBaseType::ContinuedFractionCollectionType
	                                    ContinuedFractionCollectionType;
	typedef typename DcaSolverBaseType::PlotParamsType PlotParamsType;
	typedef typename DcaSolverBaseType::MatrixType MatrixType;
	typedef typename MatrixType::value_type ComplexType;

	DcaSolverLanczos(DcaToDmrgType& myInput,
	                 const VaryingGeometryType& geometry2,
	                 typename InputNgType::Readable& io)
	 : myInput_(myInput),
	   modelSelector_(myInput,geometry2),
	   model_(modelSelector_()),
	   engine_(model_,geometry2.numberOfSites(),myInput),
	   cicjMatrix_(geometry2.numberOfSites(),geometry2.numberOfSites())
	{
		std::cout<<"DcaSolverLanczos: geometry\n";
		std::cout<<geometry2;
		RealType Eg = engine_.gsEnergy();
		std::cout.precision(8);
		std::cout<<"DcaSolverLanczos::Energy= "<<Eg<<"\n";
	}

	void solve(MatrixType& gf,
	           const VectorSizeType& sites,
	           const PlotParamsType& plotParams)
	{
		SizeType gfOp = LanczosPlusPlus::ProgramGlobals::operator2id("c");
		typename PsimagLite::Vector<PairType>::Type spins(1);
		spins[0] = PairType(0,0);

		if (sites.size() != 2)
			throw PsimagLite::RuntimeError("No sites in input file!\n");

		std::cout<<"#dca indexing: gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";
		VectorSizeType sitesLanczos(2);
		sitesLanczos[0] = myInput_.dcaIndexToDmrgIndex(sites[0]);
		sitesLanczos[1] = myInput_.dcaIndexToDmrgIndex(sites[1]);
		std::cout<<"#lanczos indexing (i="<<sitesLanczos[0];
		std::cout<<",j="<<sitesLanczos[1]<<")\n";

		SizeType matsubaras = plotParams.numberOfMatsubaras;
		PsimagLite::FreqEnum freqEnum = (matsubaras == 0) ? PsimagLite::FREQ_REAL : PsimagLite::FREQ_MATSUBARA;

		SizeType norbitals = maxOrbitals(model_);
		for (SizeType orb1 = 0;orb1 < norbitals; ++orb1) {
			SizeType orb2 = orb1;
			ContinuedFractionCollectionType cfCollection(freqEnum);
			engine_.spectralFunction(cfCollection,
			                         gfOp,
			                         sitesLanczos[0],
			                         sitesLanczos[1],
			                         spins,
			                         std::pair<SizeType,SizeType>(orb1,orb2));

			plotAll(gf,sites,orb1,cfCollection,plotParams);
		}
	}

	RealType findLowestEnergy() const
	{
		return engine_.gsEnergy();
	}

	RealType density(SizeType total1, SizeType total2) const
	{
		PsimagLite::Vector<PairType>::Type spins(1,PairType(0,0));
		SizeType cicj = LanczosPlusPlus::ProgramGlobals::operator2id("c");
		SizeType norbitals = maxOrbitals(model_);

		RealType sum = 0;
		for (SizeType orb1=0;orb1<norbitals;orb1++) {
			cicjMatrix_.setTo(0.0);
			engine_.twoPoint(cicjMatrix_,cicj,spins,PairType(orb1,orb1));
			for (SizeType i = 0; i < total2; ++i)
				sum += cicjMatrix_(i,i);
		}

		std::cout<<"WARNING: Density assumes spin sectors are equal\n";
		return 2.0 * sum;
	}

private:

	void plotAll(MatrixType& gf,
	             const VectorSizeType& sites,
	             SizeType orb,
	             const ContinuedFractionCollectionType& cfCollection,
	             const PlotParamsType& plotParams) const
	{
		SizeType Nc = gf.n_col();
		Nc = static_cast<SizeType>(sqrt(Nc));
		typename ContinuedFractionCollectionType::PlotDataType v;
		cfCollection.plot(v,plotParams);
		std::cout.precision(12);
		assert(sites.size() == 2);
		assert(gf.n_row() <= v.size());
		assert(sites[0] + sites[1]*Nc < gf.n_col());
		for (SizeType x = 0; x < gf.n_row(); ++x) {
			ComplexType tmp = lanczosFilterY(v[x].second);
			SizeType x2 = lanczosFilterX(x,gf.n_row());
			gf(x2,sites[0] + sites[1]*Nc + orb*Nc*Nc) = tmp;
		}
	}

	ComplexType lanczosFilterY(ComplexType z) const
	{
		z.real() *= (-0.25);
		z.imag() *= (-0.25);
		return z;
	}

	SizeType lanczosFilterX(int x, int total) const
	{
		assert(x < total);
		return total - x - 1;
	}

	template<typename SomeLanczosModelType>
	SizeType maxOrbitals(const SomeLanczosModelType& model) const
	{
		SizeType res=0;
		for (SizeType i=0;i<model.geometry().numberOfSites();i++) {
			if (res<model.orbitals(i)) res=model.orbitals(i);
		}

		return res;
	}

	DcaToDmrgType& myInput_;
	LanczosPlusPlus::ModelSelector<RealType,
	                               VaryingGeometryType,
	                               DcaToDmrgType> modelSelector_;
	const ModelType& model_;
	EngineType engine_;
	mutable PsimagLite::Matrix<RealType> cicjMatrix_;
};
}

#endif

