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
		static SizeType counter = 0;
		bool printNvsMu = true;
		bool printThings = (counter > 0) & printNvsMu;

		if (printThings) {
			for (RealType mu = -4.0; mu <= 0; mu += 0.1) {
				params_.mu = mu;
				functions_.updatePotentialV();
				RealType fOfMu = functions_(mu);
				std::cout<<mu<<" "<<fOfMu<<" NVSMU\n";
			}

			throw PsimagLite::RuntimeError("stop for testing\n");
		}

		counter++;

		std::cout<<"Old mu= "<<params_.mu<<" ";
		RealType mu = adjChemPot_();
		params_.mu = mu;
		functions_.updatePotentialV();

		std::cout<<"New mu= "<<params_.mu<<" ";
		if (!printThings) {
			std::cout<<"\n";
			return;
		}

		std::cout<<" n(mu)-n(target) = "<<functions_(mu)<<" ";
		std::cout<<"n(target) = "<<params_.targetDensity<<" ";
		std::cout<<" lanczos electrons Up = "<<functions_.electrons(SPIN_UP)<<" ";
		std::cout<<" lanczos electrons Down = "<<functions_.electrons(SPIN_DOWN)<<"\n";
	}

	RealType adjChemPot_() const
	{
		PairRealRealType aAndB = findAandB();
		typedef PsimagLite::RootFindingBisection<FunctionsType> RootFindingType;
		RealType tol = params_.muTolerance;
		SizeType iter = params_.muMaxIter;
		RootFindingType  rootFinding(functions_,aAndB.first, aAndB.second,iter,tol);
		RealType mu = rootFinding();
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

