#ifndef PARALLEL_DCA_SOLVER_DMRG_H
#define PARALLEL_DCA_SOLVER_DMRG_H
#include "Mpi.h"
#include "Concurrency.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "../../dmrgpp/src/Engine/ParametersDmrgSolver.h"
#include "../../dmrgpp/src/Engine/VectorWithOffsets.h"
#include "../../dmrgpp/src/Engine/MatrixVectorOnTheFly.h"
#include "../../dmrgpp/src/Engine/LeftRightSuper.h"
#include "../../dmrgpp/src/Engine/TargetingCorrectionVector.h"
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
	typedef Dmrg::ParametersDmrgSolver<RealType,typename InputNgType::Readable> ParametersDmrgSolverType;
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
	typedef Dmrg::WaveFunctionTransfFactory<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
	typedef Dmrg::TargetingCorrectionVector<PsimagLite::LanczosSolver,
	                                        MatrixVectorType,
	                                        WaveFunctionTransfType> TargetingType;
	typedef typename TargetingType::MatrixVectorType::ModelType ModelType;
	typedef typename TargetingType::TargettingParamsType TargettingParamsType;
	typedef Dmrg::DmrgSolver<TargetingType> SolverType;
	typedef typename DcaSolverBaseType::MatrixType MatrixType;
	
public:

	typedef Run RunType;
	typedef typename PsimagLite::Vector<RunType>::Type VectorRunType;

	ParallelDmrgSolver(DcaToDmrgType& myInput,
	                              const VaryingGeometryType& geometry2,
	                              typename InputNgType::Readable& io,
	                              MatrixType& gf,
	                              const VectorRunType& runs,
	                              const PlotParamsType& plotParams)
	: myInput_(myInput),
	  geometry2_(geometry2),
	  paramsDmrg_(io),
	  gf_(gf),
	  runs_(runs),
	  plotParams_(plotParams),
	  modelSelector_(paramsDmrg_.model),
	  model_(modelSelector_(paramsDmrg_,myInput_,geometry2_))
	{}

	void thread_function_(SizeType threadNum,
	                                    SizeType blockSize,
	                                    SizeType total,
	                                    pthread_mutex_t* myMutex)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = PsimagLite::Concurrency::npthreads;
		SizeType clusterSites = gf_.n_col();
		clusterSites = static_cast<SizeType>(sqrt(clusterSites));
		
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

			RealType omegaValue = plotParams_.omega1 + plotParams_.deltaOmega*run.omegaIndex;
			tsp.omega(omegaValue);

			SolverType dmrgSolver(model_,tsp,myInput_);
			dmrgSolver.main(geometry2_);

			for (SizeType site2 = 0; site2 < clusterSites; ++site2) {
				SizeType site2Dmrg = myInput_.dcaIndexToDmrgIndex(site2);
				gf_(run.omegaIndex,siteDmrg + site2Dmrg*clusterSites) += dmrgSolver.inSitu(site2Dmrg);
			}
		}

		PsimagLite::MPI::allReduce(gf_);
	}

private:

	DcaToDmrgType& myInput_;
	const VaryingGeometryType& geometry2_;
	ParametersDmrgSolverType paramsDmrg_;
	MatrixType& gf_;
	const VectorRunType& runs_;
	const PlotParamsType& plotParams_;
	Dmrg::ModelSelector<ModelType> modelSelector_;
	const ModelType& model_;
};

}

#endif

