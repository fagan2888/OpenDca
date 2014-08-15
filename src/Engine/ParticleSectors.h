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
#ifndef OPENDCA_PARTICLE_SECTORS_H
#define OPENDCA_PARTICLE_SECTORS_H
#include "Vector.h"
#include "String.h"
#include "Matrix.h"
#include "DcaLoopGlobals.h"

namespace OpenDca {

template<typename ParametersType_,typename InputNgType_>
class ParticleSectors {

	typedef typename ParametersType_::RealType RealType_;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	typedef ParametersType_ ParametersType;
	typedef InputNgType_ InputNgType;

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		return 0;
	}

	ParticleSectors(const ParametersType& params,
	                typename InputNgType::Readable& io)
	    : params_(params),
	      io_(io),
	      electronsUp_(0),
	      electronsDown_(0)
	{
		io_.readline(electronsUp_,"TargetElectronsUp=");
		io_.readline(electronsDown_,"TargetElectronsDown=");

		SizeType totalSites = params_.largeKs*(1+ params_.nofPointsInBathPerClusterPoint);
		SizeType onp1 = 2*params_.orbitals * totalSites + 1;
		muFeatureOffset_.resize(onp1);
		SizeType sum = 0;
		for (SizeType i = 0; i < onp1; ++i) {
			sum += (i+1);
			muFeatureOffset_[i] = sum;
		}
	}

	SizeType sectors() const
	{
		if (params_.dcaOptions.find("noparticlesectors") != PsimagLite::String::npos)
			return 0;

		SizeType onp1 = muFeatureOffset_.size();
		SizeType onp2 = onp1 + 1;
		return static_cast<SizeType>(onp1*onp2*0.5);
	}

	void set(SizeType ind)
	{
		SizeType onp1 = muFeatureOffset_.size();
		for (SizeType i = 0; i < onp1; ++i) {
			if (ind < muFeatureOffset_[i]) return set(ind,i);
		}
	}

	SizeType electrons(SpinEnum spin) const
	{
		return (spin == SPIN_UP) ? electronsUp_ : electronsDown_;
	}

private:

	void set(SizeType ind, SizeType o)
	{
		SizeType ind2 = (o == 0) ? 0 : ind - muFeatureOffset_[o-1];
		assert(o == 0 || ind >= muFeatureOffset_[o-1]);
		SizeType total = o;
		electronsUp_ = ind2;
		electronsDown_ = total - ind2;
		SizeType onp1 = muFeatureOffset_.size();
		assert(onp1 > 0);
		onp1--;
		assert(!(onp1 & 1));
		onp1 = static_cast<SizeType>(onp1*0.5);
		if (electronsUp_ > onp1 ||
		        electronsDown_ > onp1 ||
		        electronsDown_ > electronsUp_) {
			electronsUp_ = electronsDown_ = 0;
		}
	}

	const ParametersType& params_;
	typename InputNgType::Readable io_;
	SizeType electronsUp_;
	SizeType electronsDown_;
	VectorSizeType muFeatureOffset_;
}; // class ParticleSectors

}

#endif

