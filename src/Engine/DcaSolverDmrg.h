#ifndef DCA_SOLVER_DMRG_H
#define DCA_SOLVER_DMRG_H
#include "DcaSolverBase.h"
#include "Parallelizer.h"
#include "ParallelDmrgSolver.h"

namespace OpenDca {

template<typename DcaToDmrgType,typename VaryingGeometryType>
class DcaSolverDmrg : public DcaSolverBase<DcaToDmrgType, VaryingGeometryType> {

	typedef typename DcaToDmrgType::InputNgType InputNgType;
	typedef typename DcaToDmrgType::RealType RealType;
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
	 : myInput_(myInput),
	   geometry2_(geometry2),
	   io_(io),
	   sitesDone_(0)
	{}

	void solve(MatrixType& gf,
	           const VectorSizeType& sites,
	           const PlotParamsType& plotParams)
	{
		if (sites.size() != 2) throw PsimagLite::RuntimeError("DcaSolverDmrg\n");

		std::cout<<"#dca indexing: gf(i="<<sites[0]<<",j="<<sites[1]<<")\n";

		VectorSizeType sitesLanczos(2);
		sitesLanczos[0] = myInput_.dcaIndexToDmrgIndex(sites[0]);
		sitesLanczos[1] = myInput_.dcaIndexToDmrgIndex(sites[1]);
		std::cout<<"#lanczos indexing (i="<<sitesLanczos[0];
		std::cout<<",j="<<sitesLanczos[1]<<")\n";

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

	RealType findLowestEnergy()
	{
		throw PsimagLite::RuntimeError("findLowestEnergy\n");
	}

private:

	void solveForSite(MatrixType& gf,
	                  const VectorRunType& runs,
	                  const PlotParamsType& plotParams)
	{
		typedef ParallelDmrgSolver<DcaSolverBaseType> ParallelDmrgSolverType;
		typedef PsimagLite::Parallelizer<ParallelDmrgSolverType> ParallelizerType;
		ParallelizerType threadedSolver(PsimagLite::Concurrency::npthreads,
		                                PsimagLite::MPI::COMM_WORLD);

		ParallelDmrgSolverType helperSolver(myInput_,geometry2_,io_,gf,runs,plotParams);

		threadedSolver.loopCreate(runs.size(),helperSolver);
	}

private:

	DcaToDmrgType& myInput_;
	const VaryingGeometryType& geometry2_;
	typename InputNgType::Readable& io_;
	VectorSizeType sitesDone_;
};
}

#endif

