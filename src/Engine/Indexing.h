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
#ifndef OPENDCA_INDEXING_H
#define OPENDCA_INDEXING_H
#include "Vector.h"
#include "Matrix.h"
#include "DcaLoopGlobals.h"

namespace OpenDca {

template<typename ParametersType_,typename GeometryType>
class Indexing {

	typedef typename ParametersType_::RealType RealType_;

public:

	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef std::pair<int,int> PairIntIntType;
	typedef ParametersType_ ParametersType;

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		return 0;
	}

	Indexing(const ParametersType& params, const GeometryType& geometry)
	    : params_(params),geometry_(geometry)
	{}

	SizeType dcaToDmrg(SizeType i) const
	{
		if (params_.largeKs == 1) {
			SizeType orb = i % params_.orbitals;
			SizeType site = static_cast<SizeType>(i/params_.orbitals);
			return site + orb*geometry_.numberOfSites();
		} else {
			if (params_.orbitals > 1) {
				PsimagLite::String str("Nc>1 and orbitals>1 not supported\n");
				throw PsimagLite::RuntimeError(str);
			}
		}

		SizeType Nc = params_.largeKs;
		SizeType NcOver2 = static_cast<SizeType>(Nc/2);
		SizeType NcLy = 2;
		SizeType nBath = params_.nofPointsInBathPerClusterPoint;
		SizeType r = 0;

		if (i < Nc) {
			r = clusterIndexIndexing(i,NcLy);
			return r + NcOver2*nBath;
		}

		SizeType alpha = 0;
		getBathPoint(r,alpha,i,Nc,nBath);
		SizeType r2 = clusterIndexIndexing(r,NcLy);

		if (r2 < NcOver2)
			return r2 + alpha*nBath;

		return r2 -NcOver2 + alpha*nBath + Nc + NcOver2*nBath;
	}

	PairIntIntType bathCluster(SizeType alpha, SizeType, SizeType ind) const
	{
		SizeType bathOrbital = alpha % params_.orbitals;
		SizeType clusterOrbital = ind %  params_.orbitals;
		if (bathOrbital != clusterOrbital) return PairIntIntType(-1,-1);

		SizeType csno = static_cast<SizeType>(ind/params_.orbitals);
		SizeType bathSite = static_cast<SizeType>(alpha/params_.orbitals);
		SizeType largeKs2 = params_.largeKs * params_.largeKs;
		assert(csno < params_.largeKs);
		SizeType index = csno + csno*params_.largeKs + clusterOrbital*largeKs2;
		return PairIntIntType(bathSite,index);
	}

	SizeType lambda(SizeType ind,
	                SizeType jnd,
	                SizeType alpha,
	                SizeType lambdaSize) const
	{
		SizeType csno1 = static_cast<SizeType>(ind/params_.orbitals);
		SizeType csno2 = static_cast<SizeType>(jnd/params_.orbitals);
		SizeType largeKs2 = params_.largeKs * params_.largeKs;
		SizeType nBath = static_cast<SizeType>(lambdaSize/params_.orbitals);
		SizeType alphaCorrected = orbitalBathToBathOrbital(alpha, nBath);
		SizeType index = csno1 + csno2*params_.largeKs + alphaCorrected* largeKs2;
		assert(index < lambdaSize);
		return index;
	}

	SizeType clusterIndexIndexing(SizeType ind,SizeType ly) const
	{
		SizeType dim = geometry_.dimension();
		VectorRealType rvector(dim);
		geometry_.index2Rvector(ind,rvector);
		SizeType x = static_cast<SizeType>(rvector[0]);
		if (dim == 1)
			return x;

		assert(dim == 2);
		return static_cast<SizeType>(rvector[1]) + x*ly;

	}

	// Formula is i=clusterSize + Nb*r + alpha, given Nc, Nb and i determines alpha and r
	void getBathPoint(SizeType& r,
	                  SizeType& alpha,
	                  SizeType i,
	                  SizeType clusterSize,
	                  SizeType nBath) const
	{
		assert(i >= clusterSize);
		i-=clusterSize;
		r = static_cast<SizeType>(i/nBath);
		alpha = i-r*nBath;
	}

	bool isInBath(SizeType ind) const
	{
		bool b = (ind != 4 && ind != 5);
		b &= (ind != 6 && ind != 7);
		return b;
	}

	SizeType orbitalBathToBathOrbital(SizeType alpha, SizeType nBath) const
	{
		div_t q = div(alpha,params_.orbitals);
		return q.rem * nBath + q.quot;
	}

	const ParametersType& params_;
	const GeometryType& geometry_;
}; // class Indexing

}

#endif

