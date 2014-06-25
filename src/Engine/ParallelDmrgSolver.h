#ifndef PARALLEL_DCA_SOLVER_DMRG_H
#define PARALLEL_DCA_SOLVER_DMRG_H
#include "Mpi.h"
#include "Concurrency.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "IoSimple.h"
#include "ContinuedFractionCollection.h"
#include "../../dmrgpp/src/Engine/ParametersDmrgSolver.h"
#include "../../dmrgpp/src/Engine/VectorWithOffsets.h"
#include "../../dmrgpp/src/Engine/MatrixVectorOnTheFly.h"
#include "../../dmrgpp/src/Engine/LeftRightSuper.h"
#include "../../dmrgpp/src/Engine/TargetingCorrectionVector.h"
#include "../../dmrgpp/src/Engine/TargetingGroundState.h"
#include "../../dmrgpp/src/Engine/Operators.h"
#include "../../dmrgpp/src/Engine/BasisWithOperators.h"
#include "../../dmrgpp/src/Engine/ModelHelperLocal.h"
#include "../../dmrgpp/src/Engine/ModelBase.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"
#include "../../dmrgpp/src/Engine/ModelSelector.h"
#include "../../dmrgpp/src/Engine/WaveFunctionTransfFactory.h"
#include "../../dmrgpp/src/Engine/DmrgSolver.h"

namespace OpenDca {

struct Run {

	enum EnumType {TYPE_NORMAL, TYPE_DAGGER};

	Run(SizeType site1, SizeType index1, EnumType type1)
	: site(site1),omegaIndex(index1),dynamicDmrgType(type1)
	{}

	SizeType site;
	SizeType omegaIndex;
	EnumType dynamicDmrgType;
};

template<typename DcaSolverBaseType>
class ParallelDmrgSolver {

	typedef typename DcaSolverBaseType::PlotParamsType PlotParamsType;
	typedef typename DcaSolverBaseType::DcaToDmrgType DcaToDmrgType;
	typedef typename DcaSolverBaseType::VaryingGeometryType VaryingGeometryType;
	typedef typename DcaToDmrgType::InputNgType InputNgType;
	typedef typename DcaToDmrgType::RealType RealType;
	typedef Dmrg::ParametersDmrgSolver<RealType,
	                                   typename InputNgType::Readable> ParametersDmrgSolverType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef Dmrg::Basis<SparseMatrixType> BasisType;
	typedef Dmrg::Operators<BasisType> OperatorsType;
	typedef Dmrg::BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef Dmrg::LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef Dmrg::ModelHelperLocal<LeftRightSuperType> ModelHelperType;
	typedef Dmrg::ModelBase<ModelHelperType,
	                                  ParametersDmrgSolverType,
	                                  DcaToDmrgType,
	                                  VaryingGeometryType> ModelBaseType;
	typedef Dmrg::MatrixVectorOnTheFly<ModelBaseType> MatrixVectorType;
	typedef Dmrg::VectorWithOffsets<ComplexType> VectorWithOffsetType;
	typedef Dmrg::WaveFunctionTransfFactory<LeftRightSuperType,
	                                        VectorWithOffsetType> WaveFunctionTransfType;
	typedef Dmrg::TargetingCorrectionVector<PsimagLite::LanczosSolver,
	                                        MatrixVectorType,
	                                        WaveFunctionTransfType> TargetingType;
	typedef Dmrg::TargetingGroundState<PsimagLite::LanczosSolver,
	                                   MatrixVectorType,
	                                   WaveFunctionTransfType> TargetingGroundStateType;
	typedef typename TargetingType::MatrixVectorType::ModelType ModelType;
	typedef typename TargetingType::TargettingParamsType TargettingParamsType;
	typedef typename TargetingGroundStateType::TargettingParamsType GsParamsType;
	typedef Dmrg::DmrgSolver<TargetingType> SolverType;
	typedef Dmrg::DmrgSolver<TargetingGroundStateType> SolverGroundStateType;
	typedef typename DcaSolverBaseType::MatrixType MatrixType;
	typedef PsimagLite::TridiagonalMatrix<RealType> TridiagonalMatrixType;
	typedef PsimagLite::ContinuedFraction<TridiagonalMatrixType> ContinuedFractionType;
	typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
	ContinuedFractionCollectionType;

public:

	typedef Run RunType;
	typedef typename PsimagLite::Vector<RunType>::Type VectorRunType;

	static const SizeType freqDependent = 0;

	ParallelDmrgSolver(DcaToDmrgType& myInput,
	                              const VaryingGeometryType& geometry2,
	                              typename InputNgType::Readable& io,
	                              MatrixType& gf,
	                              const VectorRunType& runs,
	                              const PlotParamsType* plotParams)
	: myInput_(myInput),
	  geometry2_(geometry2),
	  paramsDmrg_(io),
	  gf_(gf),
	  runs_(runs),
	  plotParams_(plotParams),
	  modelSelector_(paramsDmrg_.model),
	  model_(modelSelector_(paramsDmrg_,myInput_,geometry2_))
	{
		paramsDmrg_.electronsUp = myInput_.electrons(DcaToDmrgType::SPIN_UP);
		paramsDmrg_.electronsDown = myInput_.electrons(DcaToDmrgType::SPIN_DOWN);
		GsParamsType tsp(myInput_,model_);
		SolverGroundStateType dmrgSolver(model_,tsp,myInput_);
		dmrgSolver.main(geometry2_);
		energy_ = dmrgSolver.energy();
	}

	void thread_function_(SizeType threadNum,
	                                    SizeType blockSize,
	                                    SizeType total,
	                                    pthread_mutex_t* myMutex)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = PsimagLite::Concurrency::npthreads;
		SizeType clusterSites = gf_.n_col();
		clusterSites = static_cast<SizeType>(sqrt(clusterSites));
		SizeType l = paramsDmrg_.filename.length();
		paramsDmrg_.filename = paramsDmrg_.filename.substr(0,l-4) +
		                         ttos(mpiRank) + ".txt";
		for (SizeType p=0;p<blockSize;p++) {
			SizeType px = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (px >= total || px >= runs_.size()) continue;
			const RunType& run = runs_[px];

			TargettingParamsType tsp(myInput_,model_);

			tsp.type( (run.dynamicDmrgType == RunType::TYPE_NORMAL) ? 0 : 1);
			if (run.dynamicDmrgType == RunType::TYPE_DAGGER)
				tsp.transposeConjugate(0);

			SizeType siteDmrg = myInput_.dcaIndexToDmrgIndex(run.site);
			tsp.setSite(0,siteDmrg);

			RealType omegaValue = (plotParams_) ? plotParams_->omega1 +
			                        plotParams_->deltaOmega*run.omegaIndex : 0.0;
			if (!freqDependent) omegaValue = 0;
			tsp.omega(omegaValue);

			SolverType dmrgSolver(model_,tsp,myInput_);
			dmrgSolver.main(geometry2_);

			energy_ = dmrgSolver.energy();

			accumulateGf(run,dmrgSolver);
		}

		PsimagLite::MPI::allReduce(gf_);
	}

	RealType energy() const
	{
		return energy_;
	}

private:

	void accumulateGf(const RunType& run,
	                  const SolverType& dmrgSolver)
	{
		if (gf_.n_row() == 0) return;

		SizeType clusterSites = gf_.n_col();

		if (freqDependent) {
			for (SizeType site2 = 0; site2 < clusterSites; ++site2) {
				SizeType site2Dmrg = myInput_.dcaIndexToDmrgIndex(site2);
				gf_(run.omegaIndex,
				    run.site + site2*clusterSites) += dmrgSolver.inSitu(site2Dmrg);
			}

			return;
		}

		PsimagLite::IoSimple::In io(myInput_.outputFile());
		ContinuedFractionCollectionType cf(io);
		typename ContinuedFractionCollectionType::PlotDataType v;
		cf.plot(v,*plotParams_);
		for (SizeType site2 = 0; site2 < clusterSites; ++site2) {
			for (SizeType x=0;x<v.size();x++)
				gf_(x,run.site + site2*clusterSites) += v[x].second;
		}
	}

	DcaToDmrgType& myInput_;
	const VaryingGeometryType& geometry2_;
	ParametersDmrgSolverType paramsDmrg_;
	MatrixType& gf_;
	const VectorRunType& runs_;
	const PlotParamsType* plotParams_;
	Dmrg::ModelSelector<ModelType> modelSelector_;
	const ModelType& model_;
	RealType energy_;
};

}

#endif

