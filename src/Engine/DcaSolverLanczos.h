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

namespace OpenDca {

template<typename DcaToDmrgType,typename VaryingGeometryType>
class DcaSolverLanczos : public DcaSolverBase<DcaToDmrgType, VaryingGeometryType> {

	typedef typename DcaToDmrgType::RealType RealType;
	typedef typename DcaToDmrgType::InputNgType InputNgType;
	typedef LanczosPlusPlus::ModelBase<RealType,VaryingGeometryType,DcaToDmrgType> ModelType;
	typedef LanczosPlusPlus::BasisBase<VaryingGeometryType> BasisType;
	typedef LanczosPlusPlus::DefaultSymmetry<VaryingGeometryType,BasisType> SpecialSymmetryType;
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

	DcaSolverLanczos(DcaToDmrgType& myInput,
	                 const VaryingGeometryType& geometry2,
	                 typename InputNgType::Readable& io)
	 : myInput_(myInput),
	   modelSelector_(myInput,geometry2),
	   model_(modelSelector_()),
	   engine_(model_,geometry2.numberOfSites(),myInput)
	{

		std::cout<<geometry2;
		RealType Eg = engine_.gsEnergy();
		std::cout.precision(8);
		std::cout<<"DcaSolverLanczos::Energy= "<<Eg<<"\n";
	}

	void solve(MatrixType& gf,const VectorSizeType& sites,const PlotParamsType& plotParams)
	{
		SizeType gfOp = LanczosPlusPlus::ProgramGlobals::operator2id("c");
		typename PsimagLite::Vector<PairType>::Type spins(1);
		spins[0] = PairType(0,0);

		if (sites.size() != 2) throw PsimagLite::RuntimeError("No sites in input file!\n");

		std::cout<<"#dca indexing: gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";
		VectorSizeType sitesLanczos(2);
		sitesLanczos[0] = myInput_.dcaIndexToDmrgIndex(sites[0]);
		sitesLanczos[1] = myInput_.dcaIndexToDmrgIndex(sites[1]);
		std::cout<<"#lanczos indexing (i="<<sitesLanczos[0]<<",j="<<sitesLanczos[1]<<")\n";

		ContinuedFractionCollectionType cfCollection;
		SizeType norbitals = maxOrbitals(model_);
		for (SizeType orb1 = 0;orb1 < norbitals; ++orb1) {
			for (SizeType orb2 = 0;orb2 < norbitals; ++orb2) {
				engine_.spectralFunction(cfCollection,
				                          gfOp,
				                          sitesLanczos[0],
				                          sitesLanczos[1],
				                          spins,
				                          std::pair<SizeType,SizeType>(orb1,orb2));
			}
		}

		plotAll(gf,sites,cfCollection,plotParams);
	}

private:

	void plotAll(MatrixType& gf,
	                  const VectorSizeType& sites,
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
		for (SizeType x = 0; x < gf.n_row(); ++x)
			gf(x,sites[0] + sites[1]*Nc) = v[x].second;
	}

	template<typename SomeLanczosModelType>
	SizeType maxOrbitals(const SomeLanczosModelType& model)
	{
		SizeType res=0;
		for (SizeType i=0;i<model.geometry().numberOfSites();i++) {
			if (res<model.orbitals(i)) res=model.orbitals(i);
		}

		return res;
	}

	DcaToDmrgType& myInput_;
	LanczosPlusPlus::ModelSelector<RealType,VaryingGeometryType,DcaToDmrgType> modelSelector_;
	const ModelType& model_;
	EngineType engine_;
};
}

#endif

