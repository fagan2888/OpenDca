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

