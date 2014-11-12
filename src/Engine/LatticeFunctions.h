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
#ifndef OPENDCA_LATTICE_FUNC_H
#define OPENDCA_LATTICE_FUNC_H
#include "../../../PsimagLite/src/FreqEnum.h"
#include "DcaLoopGlobals.h"

namespace OpenDca {

template<typename MatrixType,
         typename VectorType,
         typename DispersionType>
class LatticeFunctions {

	typedef typename MatrixType::value_type ComplexType;
	typedef typename DispersionType::GeometryType GeometryType;

public:

	enum PrintEnum {PRINT_YES, PRINT_NO};

	typedef typename DispersionType::ParametersType ParametersType;
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;

	LatticeFunctions(const ParametersType& params,
	                const GeometryType& geometry,
	                const DispersionType& dispersion)
	    : freqEnum_(PsimagLite::FREQ_MATSUBARA),
	      params_(params),
	      geometry_(geometry),
	      gf_(params_.numberOfMatsubaras,
	          params_.largeKs*params_.orbitals*params_.orbitals),
	      sigma_(gf_.n_row(),gf_.n_col()),
	      dispersion_(dispersion)
	{
		makeGf(freqEnum_,PRINT_YES);
	}

	// G(\tau=0) = T \sum_n G_interacting (iwn,orbital1, orbital2, site)
	// DIAGONAL IN SITE AND DIAGONAL IN ORBITAL
	RealType operator()(RealType mu) const
	{
		params_.mu = mu;
		makeGf(freqEnum_, PRINT_NO);

		ComplexType sum = 0;
		for (SizeType i = 0; i < gf_.n_row(); ++i) {
			for (SizeType j = 0; j < gf_.n_col(); ++j) { // FIXME
				sum += gf_(i,j) + correction(i,j);
			}
		}

		sum += integratedCorrection();
		RealType density = 0.5*params_.orbitals + std::real(sum)/params_.beta;
		return density - params_.targetDensity;
	}

	RealType makeSigma(const MatrixType& interacting,
	                   const MatrixType& nonInteracting,
	                   PsimagLite::FreqEnum)
	{
		for (SizeType i = 0;i < sigma_.n_row(); ++i) {
			for (SizeType j=0;j < sigma_.n_col(); ++j) {
				// j = clusterK + orb1*largeKs + orb2*largeKs*orbitals
				SizeType clusterK = j % params_.largeKs;
				SizeType tmp = static_cast<SizeType>(j/params_.largeKs);
				SizeType orb1 = static_cast<SizeType>(tmp / params_.orbitals);
				SizeType orb2 = tmp % params_.orbitals;
				if (orb1 != orb2) continue;
				SizeType jj = clusterK + orb1 * params_.largeKs;
				sigma_(i,j) = 1.0/nonInteracting(i,jj) - 1.0/interacting(i,jj);
			}
		}

		std::cout<<"#SIGMA\n";
		std::cout<<sigma_;

		return calcRealSigma(sigma_);
	}

	void setFreqType(PsimagLite::FreqEnum freqEnum)
	{
		freqEnum_ = freqEnum;
	}

	void updatePotentialV() {}

	SizeType omegaSize(PsimagLite::FreqEnum freqEnum) const
	{
		if (freqEnum == PsimagLite::FREQ_REAL)
			return params_.omegas;
		else
			return params_.numberOfMatsubaras;
	}

	ComplexType omegaValue(SizeType omegaIndex,PsimagLite::FreqEnum freqEnum) const
	{
		if (freqEnum == PsimagLite::FREQ_REAL)
			return ComplexType(params_.omegaBegin + params_.omegaStep*omegaIndex,
			                   params_.delta);
		else
			return ComplexType(0, matsubara(omegaIndex));
	}

	RealType matsubara(int ind) const
	{
		int halfNs = static_cast<int>(params_.numberOfMatsubaras*0.5);
		RealType factor = 2.0*M_PI/params_.beta;
		int ind2 = ind - halfNs;
		if (ind2 >= 0) return factor*(ind2 + 1);
		return factor*ind2;
	}

	const MatrixType& gf() const { return gf_; }

	const MatrixType& sigma() const { return sigma_; }

	void makeGf(PsimagLite::FreqEnum freqEnum, PrintEnum printEnum) const
	{
		VectorType gckf(omegaSize(freqEnum));
		SizeType Nc = params_.largeKs;
		SizeType norb = params_.orbitals;

		for (SizeType k = 0; k < Nc; ++k) {
			for (SizeType gamma1 = 0; gamma1 < norb; ++gamma1) {
				for (SizeType gamma2 = 0; gamma2 < norb; ++gamma2) {
					SizeType index = k+gamma1*Nc+gamma2*Nc*norb;
					makeGf(gckf,k,gamma1,gamma2,freqEnum);
					for (SizeType omegaIndex=0; omegaIndex<gckf.size();++omegaIndex) {
						gf_(omegaIndex,index) = gckf[omegaIndex]/
						              (1.0+sigma_(omegaIndex,index)*gckf[omegaIndex]);
					}
				}
			}
		}

		if (printEnum == PRINT_YES) {
			std::cout<<"#GCKFSC\n";
			std::cout<<gf_;
		}
	}

	SizeType electrons(SpinEnum) const
	{
		throw PsimagLite::RuntimeError("LatticeFunctions has no integer electrons\n");
	}

private:

	void makeGf(VectorType& gckf,
	            SizeType K,
	            SizeType gamma1,
	            SizeType gamma2,
	            PsimagLite::FreqEnum freqEnum) const
	{
		SizeType meshPoints = geometry_.sizeOfMesh();
		SizeType norb = params_.orbitals;
		SizeType Nc = params_.largeKs;
		SizeType index = K+gamma1*Nc+gamma2*Nc*norb;

		for (SizeType omegaIndex = 0; omegaIndex < gckf.size(); ++omegaIndex) {
			gckf[omegaIndex] = 0.0;
			if (gamma1 != gamma2) continue;
			for (SizeType ktilde = 0; ktilde < meshPoints; ++ktilde) {
				ComplexType tmp = -dispersion_(K,gamma1,gamma2,ktilde)
				                  -sigma_(omegaIndex,index);
				tmp += (params_.mu + omegaValue(omegaIndex,freqEnum));
				gckf[omegaIndex] += 1.0/tmp;
			}

			gckf[omegaIndex] /= meshPoints;
		}
	}

	RealType calcRealSigma(const MatrixType& sigma) const
	{
		RealType tmp =0.0;
		int indexOfZero=static_cast<int>((sigma.n_row()-1)*0.5);
		for (SizeType i = 0; i < params_.largeKs; ++i)
			tmp += fabs(std::real(sigma(indexOfZero,i)));

		return tmp/params_.largeKs;
	}


	// FIXME: Correct for orbitals
	RealType correction(SizeType i, SizeType j) const
	{
		RealType num = params_.potentialV[0] - params_.mu + sigmaHartree(j);
		RealType wn = matsubara(i);
		return num/(wn*wn);
	}

	// No freq. dependency here
	RealType sigmaHartree(SizeType) const
	{
		return 0.0; // FIXME
	}

	RealType integratedCorrection() const
	{
		return 0.0; // FIXME
	}

	PsimagLite::FreqEnum freqEnum_;
	const ParametersType& params_;
	const GeometryType& geometry_;
	mutable MatrixType gf_;
	MatrixType sigma_;
	const DispersionType& dispersion_;
}; // LatticeFunctions

} // namespace OpenDca

#endif // OPENDCA_LATTICE_FUNC_H

