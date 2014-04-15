#ifndef ANDERSON_FIT_H
#define ANDERSON_FIT_H
#include "Vector.h"
#include "String.h"
#include "RandomForTests.h"
#include "Minimizer.h"

namespace OpenDca {

// Anderson distance
template<class ComplexType,typename ParametersType>
class AndersonDistance
{
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;

public:

	typedef RealType FieldType;

	AndersonDistance(const ParametersType& params,const VectorType& gf)
	: params_(params),gf_(gf),matsubaras_(params_.numberOfMatsubaras)
	{
		SizeType nwn = static_cast<SizeType>(params_.numberOfMatsubaras/2);
		for (SizeType i=0;i<params_.numberOfMatsubaras;++i)
			matsubaras_[i]=(2.0*(i-nwn)+1.0)*M_PI/params_.beta;
	}

	RealType operator()(RealType* data, SizeType n1)
	{
		fillVector(data,n1);
		RealType tmp=0.0;
		SizeType totaln=matsubaras_.size();

		for (SizeType n=0;n<totaln;++n) {
			ComplexType ctmp = gf_[n];
			ComplexType ctmp2=andersonG0(p_,matsubaras_[n],params_.mu);
			RealType tmp2 = std::norm(ctmp-ctmp2);
			tmp += tmp2*tmp2;
		}

		tmp /= static_cast<RealType>(totaln+1.0);

		return tmp;
	}

	SizeType size() const { return 2*params_.nofPointsInBathPerClusterPoint; }

private:

	void fillVector(RealType* data, SizeType n)
	{
		if (p_.size() != n) p_.resize(n);

		for (SizeType i = 0; i < n; ++i)
			p_[i] = data[i];
	}

	// This is the function that defines the fitting
	ComplexType andersonG0(const VectorRealType& p,RealType wn,RealType mu)
	{
		SizeType n = static_cast<SizeType>(p.size()/2);
		ComplexType ctmp=0;

		for (SizeType i = 0;i < n; ++i)
			ctmp += p[i]*p[i]/ComplexType(-p[i+n],wn);

		return ctmp;
	}

	const ParametersType& params_;
	const VectorType& gf_;
	VectorRealType matsubaras_;
	VectorRealType p_;
};

template<typename VectorType,typename ParametersType>
class AndersonFit {

	typedef typename VectorType::value_type ComplexType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef PsimagLite::RandomForTests<RealType> RngType;
	typedef AndersonDistance<ComplexType,ParametersType> AndersonDistanceType;
	typedef PsimagLite::Minimizer<RealType,AndersonDistanceType> MinimizerType;

public:

	AndersonFit(const ParametersType& params)
	: params_(params),rng_(0)
	{}

	// Project some Gf to get the Anderson PArameters, V_p and e_p
	// p[0] to p[n-1] contains V_p
	// p[n] to p[2n-1] contains e_p
	void fit(VectorRealType& p,const VectorType& Gf,RealType integral)
	{
		bool converged = false;

		 for (RealType tolerance=1e-6;tolerance<1e-3;tolerance *=10) {
			std::fill(p.begin(), p.end(), 0.0);
			int err = aux(p,Gf,tolerance,integral);
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
	int aux(VectorRealType& p,const VectorType& Gf,RealType tolerance,RealType integral)
	{
		SizeType maxIter = 1000;
		RealType delta = 1e-3;
		AndersonDistanceType andersonDistance(params_,Gf);
		MinimizerType minimizer(andersonDistance,maxIter);
		int err = 1;

		for (SizeType i=0;i<20;++i) {
			for (SizeType k=0;k<p.size()/2;++k)
				p[k]=2.0*fabs(integral)/(M_PI*p.size());

			for (SizeType k=p.size()/2;k<p.size();++k)
				p[k] = rng_();

			err=minimizer.simplex(p,delta,tolerance);

			if (err > 0) break;
		}

		if (err > 0) {
			std::cerr<<"getAndersonParameters: converged after ";
			std::cerr<<err<<" iterations and ";
			std::cerr<<" tolerance="<<tolerance<<std::endl;
			std::cerr<<"pInitial\n";
			std::cerr<<p;
		}

		return err;
	}

	const ParametersType& params_;
	RngType rng_;
};

}

#endif

