#ifndef DCA_SOLVER_DMRG_H
#define DCA_SOLVER_DMRG_H
#include "DcaSolverBase.h"
#include "LanczosSolver.h"
#include "Parallelizer.h"
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
#include "ParallelDmrgSolver.h"

namespace OpenDca {

template<typename DcaToDmrgType,typename VaryingGeometryType>
class DcaSolverDmrg : public DcaSolverBase<DcaToDmrgType, VaryingGeometryType> {

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
	typedef Run RunType;
	typedef typename PsimagLite::Vector<RunType>::Type VectorRunType;

public:

	typedef DcaSolverBase<DcaToDmrgType, VaryingGeometryType> DcaSolverBaseType;
	typedef typename DcaSolverBaseType::VectorSizeType VectorSizeType;
	typedef typename DcaSolverBaseType::PlotParamsType PlotParamsType;
	typedef typename DcaSolverBaseType::MatrixType MatrixType;

	DcaSolverDmrg(DcaToDmrgType& myInput,
	                        const VaryingGeometryType& geometry2,
	                        typename InputNgType::Readable& io)
	 :  myInput_(myInput),
	  geometry2_(geometry2),
	  paramsDmrg_(io),
	  modelSelector_(paramsDmrg_.model),
	  model_(modelSelector_(paramsDmrg_,myInput,geometry2)),
	  tsp_(myInput,model_),
	  dmrgSolver_(model_,tsp_,myInput),
	  sitesDone_(0)
	{}

	void solve(MatrixType& gf,const VectorSizeType& sites,const PlotParamsType& plotParams)
	{
		if (sites.size() != 2) throw PsimagLite::RuntimeError("DcaSolverDmrg\n");
		
		std::cout<<"#dca indexing: gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";
		
		VectorSizeType sitesLanczos(2);
		sitesLanczos[0] = myInput_.dcaIndexToDmrgIndex(sites[0]);
		sitesLanczos[1] = myInput_.dcaIndexToDmrgIndex(sites[1]);
		std::cout<<"#lanczos indexing (i="<<sitesLanczos[0]<<",j="<<sitesLanczos[1]<<")\n";
		
		if (std::find(sitesDone_.begin(),sitesDone_.end(),sites[0]) != sitesDone_.end())
			return;

		VectorRunType runs;
		for (SizeType x = 0; x < gf.n_row(); ++x) {
			RunType run1(sites[0],x,RunType::TYPE_NORMAL);
			runs.push_back(run1);
			RunType run2(sites[0],x,RunType::TYPE_DAGGER);
			runs.push_back(run2);
		}
		
		solveForSite(gf,runs,plotParams);
		
		sitesDone_.push_back(sites[0]);
	}

private:

	void solveForSite(MatrixType& gf,const VectorRunType& runs,const PlotParamsType& plotParams)
	{
		typedef ParallelDmrgSolver<PlotParamsType,MatrixType> ParallelDmrgSolverType;
		typedef PsimagLite::Parallelizer<ParallelDmrgSolverType> ParallelizerType;
		ParallelizerType threadedSolver(PsimagLite::Concurrency::npthreads,PsimagLite::MPI::COMM_WORLD);

		ParallelDmrgSolverType helperSolver(gf,runs,plotParams);

		threadedSolver.loopCreate(runs.size(),helperSolver);
	}

private:
	
	DcaToDmrgType& myInput_;
	const VaryingGeometryType& geometry2_;
	ParametersDmrgSolverType paramsDmrg_;
	Dmrg::ModelSelector<ModelType> modelSelector_;
	const ModelType& model_;
	TargettingParamsType tsp_;
	SolverType dmrgSolver_;
	VectorSizeType sitesDone_;
};
}

#endif

