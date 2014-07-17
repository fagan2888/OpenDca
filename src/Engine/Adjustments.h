#ifndef OPENDCA_ADJUSTMENTS_H
#define OPENDCA_ADJUSTMENTS_H

namespace OpenDca {

template<typename FunctionsType>
class Adjustments {

	typedef typename FunctionsType::RealType RealType;
	typedef typename FunctionsType::ParametersType ParametersType;
	typedef std::pair<RealType,RealType> PairRealRealType;

public:

	Adjustments(FunctionsType& functions, const ParametersType& params)
	    : functions_(functions), params_(params)
	{}

	void adjChemPot() const
	{
		RealType mu = adjChemPot_();
		std::cout<<"Old mu= "<<params_.mu<<" ";
		params_.mu = mu;
		functions_.updatePotentialV();
		std::cout<<"New mu= "<<params_.mu<<"\n";
	}

	RealType adjChemPot_() const
	{
		PairRealRealType aAndB = findAandB();
		for (RealType tolerance = 1e-3; tolerance < 1; tolerance *= 2) {
			try {
				RealType mu = adjChemPot_(aAndB, tolerance);
				return mu;
			} catch (std::exception& e) {}
		}

		throw PsimagLite::RuntimeError("adjChemPot failed\n");
	}

	RealType adjChemPot_(const PairRealRealType& aAndB, RealType tol) const
	{
		typedef PsimagLite::RootFindingBisection<FunctionsType> RootFindingType;
		RootFindingType  rootFinding(functions_,aAndB.first, aAndB.second,1000,tol);
		RealType mu = params_.mu;
		rootFinding(mu);
		return mu;
	}

	PairRealRealType findAandB() const
	{
		for (RealType value = 1.0; value < 100.0; value++) {
			RealType value2 = functions_(value) * functions_(-value);
			if (value2 < 0) return PairRealRealType(-value,value);
		}

		throw PsimagLite::RuntimeError("RootFinding init failed\n");
	}

	bool isOption(PsimagLite::String what) const
	{
		return (params_.dcaOptions.find(what) != PsimagLite::String::npos);
	}

private:

	FunctionsType& functions_;
	const ParametersType& params_;
}; // class Adjustments

} // namespace OpenDca

#endif // OPENDCA_ADJUSTMENTS_H

