#ifndef PARALLEL_DCA_SOLVER_DMRG_H
#define PARALLEL_DCA_SOLVER_DMRG_H
#include "Mpi.h"
#include "Concurrency.h"

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

template<typename PlotParamsType,typename MatrixType>
class ParallelDmrgSolver {

public:

	typedef Run RunType;
	typedef typename PsimagLite::Vector<RunType>::Type VectorRunType;

	ParallelDmrgSolver(MatrixType& gf,const VectorRunType& runs,const PlotParamsType& plotParams)
	{}

	void thread_function_(SizeType threadNum,
	                                    SizeType blockSize,
	                                    SizeType total,
	                                    pthread_mutex_t* myMutex)
	{
		throw PsimagLite::RuntimeError("testing");
	}

private:
};

}

#endif

