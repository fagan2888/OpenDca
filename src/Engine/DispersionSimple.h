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
#ifndef DCA_DISPERSION_SIMPLE_H
#define DCA_DISPERSION_SIMPLE_H
#include "Matrix.h"
#include "Vector.h"
#include "Concurrency.h"
#include "DcaParameters.h"

namespace OpenDca {

template<typename RealType,typename GeometryType_>
class DispersionSimple {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	typedef GeometryType_ GeometryType;
	typedef DcaParameters<RealType> ParametersType;

	DispersionSimple(const ParametersType& params, const GeometryType& geometry)
	: params_(params),geometry_(geometry)
	{}

	RealType operator()(SizeType coarseIndex,
	                    SizeType gamma1,
	                    SizeType gamma2,
	                    SizeType fineIndex) const
	{
		if (gamma1 != gamma2) return 0.0;

		VectorRealType kvector(geometry_.dimension());
		geometry_.index2Kvector(coarseIndex,kvector);

		VectorRealType kvector2(geometry_.dimension());
		geometry_.getMeshVector(kvector2,fineIndex);

		for (SizeType i=0;i<kvector.size();++i)
			kvector[i] += kvector2[i];

		return getDispersion(kvector)
		       + onSitePotential(coarseIndex,gamma1,gamma2);
	}

private:

	RealType getDispersion(const VectorRealType& kvector) const
	{
		SizeType ndim=geometry_.dimension();
		RealType kplus = 0;

		for (SizeType k=0;k<ndim;++k)
			kplus += cos(kvector[k]);

		return -2.0 * kplus;
	}

	RealType onSitePotential(SizeType K, SizeType gamma1, SizeType gamma2) const
	{
		if (gamma1 != gamma2) return 0;
		SizeType index = K + gamma1*geometry_.numberOfSites();
		assert(index< params_.potentialV.size());
		return params_.potentialV[index];
	}

	const ParametersType& params_;
	const GeometryType& geometry_;
};

} // namespace OpenDca

#endif // DCA_DISPERSION_SIMPLE_H

