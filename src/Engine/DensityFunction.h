#ifndef OPENDCA_DENSITY_FUNC_H
#define OPENDCA_DENSITY_FUNC_H

template<typename MatrixType, typename ParamsType>
class DensityFunction {

	typedef typename MatrixType::value_type ComplexType;

public:

	typedef typename PsimagLite::Real<ComplexType>::Type RealType;

	DensityFunction(const MatrixType& interacting,
	                const ParamsType& params,
	                RealType targetDensity)
	    : interacting_(interacting), params_(params), targetDensity_(targetDensity)
	{}

	// G(\tau=0) = T \sum_n G_interacting (iwn)
	RealType operator()(RealType mu) const
	{
		ComplexType sum = 0;
		for (SizeType i = 0; i < interacting_.n_row(); ++i) {
			for (SizeType j = 0; j < interacting_.n_col(); ++j) {
				sum += interacting_(i,j);
			}
		}

		RealType density = -std::imag(sum)/params_.beta;
		return density - targetDensity_;
	}

private:

	const MatrixType& interacting_;
	const ParamsType& params_;
	RealType targetDensity_;
}; // DensityFunction

#endif // OPENDCA_DENSITY_FUNC_H

