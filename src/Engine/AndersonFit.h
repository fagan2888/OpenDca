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
#ifndef ANDERSON_FIT_H
#define ANDERSON_FIT_H
#include "Vector.h"
#include "MersenneTwister.h"
#include "MultiMin.h"

namespace OpenDca {

// Anderson distance
template<class ComplexType,typename ParametersType>
class AndersonDistance
{
	typedef typename PsimagLite::Real<ComplexType>::Type RealType_;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;

public:

	typedef RealType_ RealType;
	typedef RealType FieldType;

	AndersonDistance(const ParametersType& params,const VectorType& gf)
	: params_(params),
	  gf_(gf),
	  matsubaras_(params_.numberOfMatsubaras)
	{
		SizeType nwn = static_cast<SizeType>(params_.numberOfMatsubaras/2);
		for (SizeType i=0;i<params_.numberOfMatsubaras;++i)
			matsubaras_[i]=(2.0*(i-nwn)+1.0)*M_PI/params_.beta;
	}

	RealType operator()(const RealType* data, SizeType n1) const
	{
		fillVector(data,n1);
		RealType tmp=0.0;
		SizeType totaln=matsubaras_.size();
		SizeType effectiveNs = 0;

		for (SizeType n=0;n<totaln;++n) {
			if (!omegaCondition(n)) continue;
			effectiveNs++;
			ComplexType ctmp = gf_[n];
			ComplexType ctmp2=andersonG0(p_,matsubaras_[n]);
			ComplexType z = ctmp - ctmp2;
			RealType tmp2 = std::real(z)*std::real(z) + std::imag(z)*std::imag(z);
			tmp += tmp2;
		}

		tmp /= static_cast<RealType>(effectiveNs+1.0);

		return tmp;
	}

	RealType function(const RealType* data, SizeType n1) const
	{
		return operator()(data,n1);
	}

	void gradient(RealType* result,
	              SizeType n0,
	              const RealType* data,
	              SizeType n1) const
	{
		fillVector(data,n1);
		assert(n0 == n1);
		SizeType totaln=matsubaras_.size();
		SizeType effectiveNs = 0;

		for (SizeType i = 0; i < n1; ++i) {
			RealType tmp=0.0;
			for (SizeType n=0;n<totaln;++n) {
				if (!omegaCondition(n)) continue;
				if (i == 0) effectiveNs++;
				ComplexType ctmp = gf_[n];
				ComplexType ctmp2=andersonG0(p_,matsubaras_[n]);
				ComplexType tmp2 = ctmp - ctmp2;
				ComplexType g0Gradient = andersonG0gradient(p_,
				                                            matsubaras_[n],
				                                            i);
				tmp += (std::real(tmp2) * std::real(g0Gradient) +
				        std::imag(tmp2) * std::imag(g0Gradient));
			}

			tmp *= -2.0;
			tmp /= static_cast<RealType>(effectiveNs+1.0);
			result[i] = tmp;
		}
	}

	SizeType size() const { return 2*params_.nofPointsInBathPerClusterPoint; }

private:

	void fillVector(const RealType* data, SizeType n) const
	{
		if (p_.size() != n) p_.resize(n);

		for (SizeType i = 0; i < n; ++i)
			p_[i] = data[i];
	}

	// This is the function that defines the fitting
	ComplexType andersonG0(const VectorRealType& p,RealType wn) const
	{
		SizeType n = static_cast<SizeType>(p.size()/2);
		ComplexType ctmp=0;

		for (SizeType i = 0;i < n; ++i)
			ctmp += p[i]*p[i]/ComplexType(-p[i+n],wn);

		return ctmp;
	}

	// This is the function that defines the fitting
	ComplexType andersonG0gradient(const VectorRealType& p,
	                               RealType wn,
	                               SizeType ind) const
	{
		SizeType n = static_cast<SizeType>(p.size()/2);

		if (ind < n) return 2.0*p[ind]/ComplexType(-p[ind+n],wn);

		SizeType ind2 = ind - n;
		assert(ind2 < p.size());
		ComplexType ctmp = ComplexType(-p[ind],wn);
		return p[ind2]*p[ind2]/(ctmp*ctmp);
	}

	bool omegaCondition(SizeType n) const
	{
		return (fabs(matsubaras_[n]) < params_.andersonFitCutoff);
	}

	const ParametersType& params_;
	const VectorType& gf_;
	VectorRealType matsubaras_;
	mutable VectorRealType p_;
};

template<typename VectorType,typename ParametersType>
class AndersonFit {

	typedef typename VectorType::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef PsimagLite::MersenneTwister RngType;
	typedef AndersonDistance<ComplexType,ParametersType> AndersonDistanceType;
	typedef PsimagLite::MultiMin<AndersonDistanceType> MinimizerType;

public:

	AndersonFit(const ParametersType& params)
	: params_(params),rng_(0),randomize_(false)
	{
		PsimagLite::String opt = params_.dcaOptions;
		randomize_ = isOption("FittingRandomize");
		if (!randomize_) {
			randomize_ = isOption("RandomizeFitting");
		}
	}

	// Project some Gf to get the Anderson PArameters, V_p and e_p
	// p[0] to p[n-1] contains V_p
	// p[n] to p[2n-1] contains e_p
	void fit(VectorRealType& p,const VectorType& Gf)
	{
		bool converged = false;

		 for (RealType tolerance=1e-12;tolerance<1;tolerance *=10) {
			std::fill(p.begin(), p.end(), 0.0);
			int err = aux(p,Gf,tolerance);
			if (err > 0) {
				converged = true;
				break;
			} else {
				std::cerr<<"AndersonFit: did not converge, trying with tolerance = ";
				std::cerr<<tolerance<<"\n";
			}
		}

		 if (converged) return;

		PsimagLite::String str("AndersonFit: Too many iterations\n");
		throw PsimagLite::RuntimeError(str);
	}

private:

	// Project some Gf to get the Anderson PArameters, V_p and e_p
	// p[0] to p[n-1] contains V_p
	// p[n] to p[2n-1] contains e_p
	int aux(VectorRealType& p,const VectorType& Gf,RealType tolerance)
	{
		SizeType maxIter = params_.andersonFitMaxIter;
		RealType delta = params_.andersonFitDelta;
		RealType maxGradient = params_.andersonFitMaxGradient;
		AndersonDistanceType andersonDistance(params_,Gf);
		MinimizerType minimizer(andersonDistance,maxIter,maxGradient);

		setInitialP(p);

		int err=minimizer.conjugateFr(p,delta,tolerance);

		if (err > 0) {
			std::cerr<<"getAndersonParameters: converged after ";
			std::cerr<<err<<" iterations and ";
			std::cerr<<" tolerance="<<tolerance<<std::endl;
			std::cerr<<"pInitial\n";
			std::cerr<<p;
		}

		return err;
	}

	void setInitialP(VectorRealType& p)
	{
		for (SizeType k=0;k<p.size()/2;++k)
			p[k]=1.0;

		for (SizeType k=p.size()/2;k<p.size();++k)
			p[k] = 0.0;

		if (!randomize_) return;

		for (SizeType k=0;k<p.size()/2;++k)
			p[k] = 2.0 * (rng_() - 0.5);

		for (SizeType k=p.size()/2;k<p.size();++k)
			p[k] = rng_() - 0.5;
	}

	bool isOption(PsimagLite::String what) const
	{
		return (params_.dcaOptions.find(what) != PsimagLite::String::npos);
	}

	const ParametersType& params_;
	RngType rng_;
	bool randomize_;
};

}

#endif

