/*
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
	typedef ParallelDmrgSolver<DcaSolverBaseType> ParallelDmrgSolverType;

	static const SizeType freqDependent =  ParallelDmrgSolverType::freqDependent;

	DcaSolverDmrg(DcaToDmrgType& myInput,
	              const VaryingGeometryType& geometry2,
	              typename InputNgType::Readable& io)
	 : myInput_(myInput),
	   geometry2_(geometry2),
	   io_(io),
	   sitesDone_(0),
	   energy_(0.0)
	{
		MatrixType gf(0,0);
		VectorRunType runs;
		DcaToDmrgType myInput2 = myInput_;
		ParallelDmrgSolverType helperSolver(myInput2,geometry2_,gf,runs,0);
		energy_ = helperSolver.energy();
	}

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

		SizeType orbitals = 1;
		io_.readline(orbitals,"Orbitals=");
		VectorRunType runs;
		SizeType total = (freqDependent) ? gf.n_row() : 1;
		for (SizeType orb = 0; orb < orbitals; ++orb) {
			for (SizeType x = 0; x < total; ++x) {
				RunType run1(sites[0],orb,x,RunType::TYPE_NORMAL);
				runs.push_back(run1);
				RunType run2(sites[0],orb,x,RunType::TYPE_DAGGER);
				runs.push_back(run2);
			}
		}

		solveForSite(gf,runs,&plotParams);

		sitesDone_.push_back(sites[0]);
	}

	RealType findLowestEnergy() const
	{
		return energy_;
	}

	RealType density(SizeType total1, SizeType total2) const
	{
		PsimagLite::String msg("DcaSolverDmrg::");
		msg += "density: unimplemented\n";
		throw PsimagLite::RuntimeError(msg);
	}

private:

	void solveForSite(MatrixType& gf,
	                  const VectorRunType& runs,
	                  const PlotParamsType* plotParamsPtr)
	{
		typedef PsimagLite::Parallelizer<ParallelDmrgSolverType> ParallelizerType;
		ParallelizerType threadedSolver(PsimagLite::Concurrency::npthreads,
		                                PsimagLite::MPI::COMM_WORLD);

		DcaToDmrgType myInput = myInput_;
		ParallelDmrgSolverType helperSolver(myInput,geometry2_,gf,runs,plotParamsPtr);

		threadedSolver.loopCreate(runs.size(),helperSolver);

		energy_ = helperSolver.energy();
	}

private:

	DcaToDmrgType& myInput_;
	const VaryingGeometryType& geometry2_;
	typename InputNgType::Readable& io_;
	VectorSizeType sitesDone_;
	RealType energy_;
};
}

#endif

