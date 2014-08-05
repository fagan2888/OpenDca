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

	Run(SizeType site1, SizeType orbital1, SizeType index1, EnumType type1)
	: site(site1),orbital(orbital1),omegaIndex(index1),dynamicDmrgType(type1)
	{}

	SizeType site;
	SizeType orbital;
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
	                                   DcaToDmrgType> ParametersDmrgSolverType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef Dmrg::Basis<SparseMatrixType> BasisType;
	typedef Dmrg::Operators<BasisType> OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
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

public:

	typedef Run RunType;
	typedef typename PsimagLite::Vector<RunType>::Type VectorRunType;

	static const SizeType freqDependent = 0;

	ParallelDmrgSolver(DcaToDmrgType& myInput,
	                   const VaryingGeometryType& geometry2,
	                   MatrixType& gf,
	                   const VectorRunType& runs,
	                   const PlotParamsType* plotParams)
	: myInput_(myInput),
	  geometry2_(geometry2),
	  paramsDmrg_(myInput),
	  gf_(gf),
	  runs_(runs),
	  plotParams_(plotParams),
	  modelSelector_(paramsDmrg_.model),
	  model_(modelSelector_(paramsDmrg_,myInput_,geometry2_)),
	  tsp_(myInput_,model_)
	{
		if (geometry2.label(0) != "star")
			throw PsimagLite::RuntimeError("ParallelDmrgSolver: only geometry star\n");

		paramsDmrg_.electronsUp = myInput_.electrons(SPIN_UP);
		paramsDmrg_.electronsDown = myInput_.electrons(SPIN_DOWN);
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


			MatrixType matC = model_.naturalOperator("c",0,run.orbital);
			SizeType n = matC.n_row();
			MatrixType matZero(n,n);
			if (run.dynamicDmrgType == RunType::TYPE_DAGGER)
				matC = transposeConjugate(matC);
			Su2RelatedType su2Related1;
			Su2RelatedType su2Related2;
			OperatorType opZero(matZero,
			                    -1,
			                    PairType(0,0),
			                    1,
			                    su2Related1);
			OperatorType opC(matC,
			                 -1,
			                 PairType(0,0),
			                 1,
			                 su2Related2);

			tsp_.setOperator(0,1,opZero);
			tsp_.setOperator(1,0,opC);

			tsp_.type( (run.dynamicDmrgType == RunType::TYPE_NORMAL) ? 0 : 1);

			assert(run.site == 0);
			//SizeType siteDmrg = myInput_.dcaIndexToDmrgIndex(run.site);
			//tsp_.setSite(1,siteDmrg);

			RealType omegaValue = (plotParams_) ? plotParams_->omega1 +
			                        plotParams_->deltaOmega*run.omegaIndex : 0.0;
			if (!freqDependent) omegaValue = 0;
			tsp_.omega(omegaValue);

			SolverType dmrgSolver(model_,tsp_,myInput_);
			dmrgSolver.main(geometry2_);

			energy_ = dmrgSolver.energy();

			accumulateGf(run,dmrgSolver,paramsDmrg_.filename);
		}

		PsimagLite::MPI::allReduce(gf_);
	}

	RealType energy() const
	{
		return energy_;
	}

private:

	void accumulateGf(const RunType& run,
	                  const SolverType& dmrgSolver,
	                  PsimagLite::String filename)
	{
		if (gf_.n_row() == 0) return;

		SizeType Nc = 1;

		if (freqDependent) {
			SizeType site2Dmrg = myInput_.dcaIndexToDmrgIndex(run.site);
			gf_(run.omegaIndex,run.site + run.site*Nc + run.orbital*Nc*Nc) +=
			        dmrgSolver.inSitu(site2Dmrg);

			return;
		}

		PsimagLite::String str = "#Avector";
		PsimagLite::IoSimple::In io(filename);
		io.advance(str,PsimagLite::IoSimple::In::LAST_INSTANCE);
		ContinuedFractionType cf(io);
		typename ContinuedFractionType::PlotDataType v;
		cf.plot(v,*plotParams_);
		std::cout<<"ParallelDmrgSolver run.orbital = "<<run.orbital;
		std::cout<<" run.site = "<<run.site<<"\n";
		for (SizeType x=0;x<v.size();x++) {
			ComplexType tmp = dmrgFilterY(v[x].second);
			SizeType x2 = dmrgFilterX(x,gf_.n_row());
			gf_(x2,run.site + run.site*Nc + run.orbital*Nc*Nc) += tmp;
		}
	}

	SizeType dmrgFilterX(int x, int total) const
	{
		assert(x < total);
		return total - x - 1;
	}

	ComplexType dmrgFilterY(ComplexType z) const
	{
		z.real() *= (-1);
		z.imag() *= (1);
		return z;
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
	TargettingParamsType tsp_;
};

}

#endif

