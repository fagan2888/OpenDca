#ifndef MULTIMIN_H_
#define MULTIMIN_H_

namespace PsimagLite {

#ifdef USE_GSL

#include <gsl/gsl_multimin.h>
template<typename FunctionType>
class MultiMin {

	typedef typename FunctionType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	MultiMin(const FunctionType& f, SizeType maxIter)
	: f_(f),maxIter_(maxIter),s_(0)
	{
		t_ = gsl_multimin_fdfminimizer_conjugate_fr;
		s_ = gsl_multimin_fdfminimizer_alloc (t_, f_.size());
	}

	~MultiMin()
	{
		if (s_)
			gsl_multimin_fdfminimizer_free (s_);
	}

	RealType conjugateFr(VectorRealType& p,RealType delta,RealType tolerance) const
	{
		gsl_multimin_function_fdf myFunc;
		//~ fdf.f = ;
		//~ fdf.df=;
		//~ fdf.fdf=;
		//~ fdf.n=;
		//~ fdf.params=; 
		gsl_vector* x = gsl_vector_alloc(myFunc.n);
		for (SizeType i = 0; i < p.size(); ++i)   gsl_vector_set (x, i, p[i]);

		gsl_multimin_fdfminimizer_set (s_, &myFunc, x, delta, tolerance);
		SizeType iter = 0;
		int status = 0;
		do {
			iter++;
			status = gsl_multimin_fdfminimizer_iterate (s_);

			if (status) break;

			status = gsl_multimin_test_gradient (s_->gradient, 1e-3);
			if (status == GSL_SUCCESS) std::cout<<"Minimum found\n";

		} while (status == GSL_CONTINUE && iter < maxIter_);

		for (SizeType i = 0; i < p.size(); ++i)  p[i] = gsl_vector_get (s_->x, i);
		gsl_vector_free(x);

		return (status) ? -1 : iter;
	}

private:

	//~ double f (const gsl_vector * x, void * params)
	//~ {
	//~ }

	//~ void  df(const gsl_vector * x, void * params, gsl_vector * g)
	//~ {
	//~ }

	//~ void fdf (const gsl_vector * x, void * params, double * f, gsl_vector * g)
	//~ {
	//~ }

	const FunctionType& f_;
	SizeType maxIter_;
	const gsl_multimin_fdfminimizer_type *t_;
	gsl_multimin_fdfminimizer *s_;

};
#else
template<typename FunctionType>
class MultiMin {

	typedef typename FunctionType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	MultiMin(const FunctionType& f, SizeType maxIter)
	: f_(f),maxIter_(maxIter)
	{
		throw PsimagLite::RuntimeError("MultiMin needs USE_GSL\n");
	}

	RealType conjugateFr(VectorRealType& p,RealType delta,RealType tolerance) const
	{
		throw PsimagLite::RuntimeError("MultiMin needs USE_GSL\n");
	}

private:

	const FunctionType& f_;
	SizeType maxIter_;

};
#endif

}
#endif

