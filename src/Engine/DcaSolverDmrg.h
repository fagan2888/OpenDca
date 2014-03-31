#ifndef DCA_SOLVER_DMRG_H
#define DCA_SOLVER_DMRG_H
#include "DcaSolverBase.h"
#include "LanczosSolver.h"
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

public:

	typedef DcaSolverBase<DcaToDmrgType, VaryingGeometryType> DcaSolverBaseType;
	typedef typename DcaSolverBaseType::VectorSizeType VectorSizeType;
	typedef typename DcaSolverBaseType::PlotParamsType PlotParamsType;
	typedef typename DcaSolverBaseType::MatrixType MatrixType;

	DcaSolverDmrg(DcaToDmrgType& myInput,
	              const VaryingGeometryType& geometry2,
	              typename InputNgType::Readable& io)
	 : geometry2_(geometry2),
	  paramsDmrg_(io),
	  modelSelector_(paramsDmrg_.model),
	  model_(modelSelector_(paramsDmrg_,myInput,geometry2)),
	  tsp_(myInput,model_),
	  dmrgSolver_(model_,tsp_,myInput)
	{
	}

	void solve(MatrixType& gf,const VectorSizeType& sites,const PlotParamsType& plotParams)
	{
		dmrgSolver_.main(geometry2_);
	}

private:

	const VaryingGeometryType& geometry2_;
	ParametersDmrgSolverType paramsDmrg_;
	Dmrg::ModelSelector<ModelType> modelSelector_;
	const ModelType& model_;
	TargettingParamsType tsp_;
	SolverType dmrgSolver_;
};
}

#endif

