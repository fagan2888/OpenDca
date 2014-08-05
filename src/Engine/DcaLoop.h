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
#ifndef DCA_LOOP_H
#define DCA_LOOP_H
#include "Matrix.h"
#include "Vector.h"
#include "Concurrency.h"
#include "EffectiveHamiltonian.h"
#include "../../../PsimagLite/src/FreqEnum.h"
#include "LatticeFunctions.h"
#include "RootFindingBisection.h"
#include "Adjustments.h"

namespace OpenDca {

template<typename DispersionType,typename InputNgType>
class DcaLoop {

	typedef typename DispersionType::ParametersType ParametersType;
	typedef typename ParametersType::RealType RealType;
	typedef typename DispersionType::GeometryType GeometryType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef EffectiveHamiltonian<ParametersType,
	                             GeometryType,
	                             InputNgType> EffectiveHamiltonianType;
	typedef typename EffectiveHamiltonianType::DcaToDmrgType DcaToDmrgType;
	typedef LatticeFunctions<MatrixType,
	                        VectorType,
	                        DispersionType> LatticeFunctionsType;
	typedef Adjustments<LatticeFunctionsType> AdjustmentsType;

public:

	DcaLoop(const ParametersType& params,
	        typename InputNgType::Readable& io)
	: params_(params),
	  io_(io),
	  geometry_(io,false,params_.smallKs),
	  dispersion_(params,geometry_),
	  fTCoefsR2K_(params.largeKs,params.largeKs),
	  latticeFunctions_(params_,geometry_,dispersion_),
	  adjustments_(latticeFunctions_,params_)
	{
		SizeType Nc = params_.largeKs;
		SizeType dim = geometry_.dimension();
		VectorRealType tmpVec(dim);
		VectorRealType tmpVec2(dim);

		for (SizeType ic=0;ic<Nc;ic++) {
			geometry_.index2Rvector(ic,tmpVec);
			for (SizeType ik=0;ik<Nc;ik++) {
				geometry_.index2Kvector(ik,tmpVec2);
				RealType tmp = scalarProduct(tmpVec2,tmpVec);
				fTCoefsR2K_(ik,ic)=ComplexType(cos(tmp),sin(tmp));
			}
		}

		disableMpiInDmrgPlusPlus();
	}

	void main(PsimagLite::FreqEnum freqEnum,SizeType iterations)
	{
		SizeType largeKs = params_.largeKs;
		SizeType norb = params_.orbitals;
		VectorRealType barEpsilon(largeKs*norb*norb);
		coarseDispersion(barEpsilon);
		std::cout<<"barEpsilon\n";
		std::cout<<barEpsilon;
		bool lanczosReal = adjustments_.isOption("lanczosreal");
		bool adjustMuLattice = (!adjustments_.isOption("adjustmucluster") &
		                        !adjustments_.isOption("noadjustmu"));
		SizeType omegasRealOrImag = (lanczosReal) ? params_.omegas :
		                                            params_.numberOfMatsubaras;
		MatrixType gfcluster(omegasRealOrImag,largeKs*largeKs*params_.orbitals);
		MatrixType gammaOmegaRealOrImag(latticeFunctions_.omegaSize(freqEnum),largeKs*norb);
		MatrixType G0inverse(omegasRealOrImag,largeKs*norb);
		MatrixType gfclusterK(omegasRealOrImag,params_.largeKs*params_.orbitals);
		RealType sigmaNorm = 0;
		latticeFunctions_.setFreqType(freqEnum);

		for (SizeType i = 0; i < iterations; ++i) {
			latticeFunctions_.makeGf(freqEnum, LatticeFunctionsType::PRINT_YES);
			if (adjustMuLattice) adjustments_.adjChemPot();

			diagUpdate(gfcluster,gammaOmegaRealOrImag,barEpsilon,freqEnum);
			std::cout<<"#gfcluster\n";
			std::cout<<gfcluster;

			// transform from real to k space
			ft4(gfclusterK,gfcluster,fTCoefsR2K_,params_.largeKs);

			makeG0(G0inverse,gammaOmegaRealOrImag,barEpsilon,freqEnum);

			RealType dcaError = sigmaNorm;
			sigmaNorm = latticeFunctions_.makeSigma(gfclusterK,G0inverse,freqEnum);
			dcaError -= sigmaNorm;

			std::cout<<"sigma\n";
			std::cout<<latticeFunctions_.sigma();
			std::cout<<"Dca iteration= "<<i<<" error in sigma= "<<fabs(dcaError)<<"\n";
		}
	}

private:

	void diagUpdate(MatrixType& gfCluster,
	                MatrixType& gammaOmegaRealOrImag,
	                const VectorRealType& ekbar,
	                PsimagLite::FreqEnum freqEnum)
	{
		const MatrixType& gckfsc = latticeFunctions_.gf();
		bool lanczosReal = adjustments_.isOption("lanczosreal");
		VectorRealType integral(ekbar.size());
		MatrixType deltaOmega(gckfsc.n_row(),gckfsc.n_col());
		EffectiveHamiltonianType effectiveHamiltonian(params_,geometry_,io_);

		std::cerr<<"gamma omega...\n";
		getDeltaOmega(deltaOmega,integral,ekbar,freqEnum);

		MatrixType gammaOmega(params_.numberOfMatsubaras,gckfsc.n_col());
		getGammaOmega(gammaOmega,deltaOmega,freqEnum);
		std::cerr<<"make h params...\n";
		effectiveHamiltonian.build(gammaOmega,ekbar,integral,freqEnum);
		std::cout<<effectiveHamiltonian;

		//throw PsimagLite::RuntimeError("testing\n");
		std::cerr<<"lanczos started...\n";
		effectiveHamiltonian.solve(gfCluster);
		std::cerr<<"lanczos done\n";

		MatrixType* gfClusterMatsubara = 0;

		if (lanczosReal) {
			gfClusterMatsubara = new MatrixType(gckfsc.n_row(),gckfsc.n_col());
			hilbertTransfFromReal(*gfClusterMatsubara,gfCluster);
		} else {
			gfClusterMatsubara = &gfCluster;
		}

		std::cout<<"#gfMatsubara\n";
		for (SizeType i=0;i<gfClusterMatsubara->n_row();++i) {
			RealType wn = latticeFunctions_.matsubara(i);
			std::cout<<wn<<" ";
			for (SizeType j=0;j<gfClusterMatsubara->n_col();++j)
				std::cout<<gfClusterMatsubara->operator()(i,j)<<" ";
			std::cout<<"\n";
		}

		if (lanczosReal) delete gfClusterMatsubara;

		const MatrixType& p = effectiveHamiltonian.andersonParameters();

		getGammaKOmega(gammaOmegaRealOrImag,p,freqEnum);
	}

	void coarseDispersion(VectorRealType& barEpsilon)
	{
		SizeType meshPoints = geometry_.sizeOfMesh();
		SizeType largeKs = params_.largeKs;
		SizeType norb = params_.orbitals;
		assert(barEpsilon.size() == largeKs*norb*norb);

		for (SizeType bigK = 0; bigK < largeKs; ++bigK) {
			for (SizeType gamma1 = 0; gamma1 < norb; ++gamma1) {
				for (SizeType gamma2 = 0; gamma2 < norb; ++gamma2) {
					SizeType index = bigK + gamma1*largeKs + gamma2*largeKs*norb;
					barEpsilon[index] = 0.0;
					for (SizeType ktilde = 0; ktilde < meshPoints; ++ktilde)
						barEpsilon[index] += dispersion_(bigK,gamma1,gamma2,ktilde);
					barEpsilon[index] /= meshPoints;
				}
			}
		}
	}

	//! Delta_k(omega) = omega + ekbar(k) -sigma(omega,k) - 1.0/gckfsc_(omega,k)
	void getDeltaOmega(MatrixType& deltaOmega,
	                   VectorRealType& integral,
	                   const VectorRealType& ekbar,
	                   PsimagLite::FreqEnum freqEnum)
	{
		SizeType norb = params_.orbitals;
		SizeType largeKs = params_.largeKs;
		const MatrixType& gckfsc = latticeFunctions_.gf();

		for (SizeType j=0;j<gckfsc.n_col();j++) integral[j]=0.0;

		for (SizeType i=0;i<gckfsc.n_row();++i) { // loop over freq.
			ComplexType omega = latticeFunctions_.omegaValue(i,freqEnum);
			for (SizeType bigK = 0; bigK < largeKs; ++bigK) {
				for (SizeType gamma1 = 0; gamma1 < norb; ++gamma1) {
					for (SizeType gamma2 = 0; gamma2 < norb; ++gamma2) {
						SizeType index = bigK + gamma1*largeKs + gamma2*largeKs*norb;
						deltaOmega(i,index) = 0.0;
						if (gamma1 != gamma2) continue;
						assert(std::norm(gckfsc(i,index)) > 1e-20);
						deltaOmega(i,index)= omega + params_.mu - ekbar[index]
						      -latticeFunctions_.sigma()(i,index) - 1.0/gckfsc(i,index);
						integral[index] += std::imag(deltaOmega(i,index));
					}
				}
			}
		}

		std::cout<<"#DELTAOMEGA\n";
		for (SizeType i=0;i<gckfsc.n_row();++i) { // loop over freq.
			ComplexType omega = latticeFunctions_.omegaValue(i,freqEnum);
			std::cout<<std::real(omega)<<" ";
			for (SizeType j=0;j<gckfsc.n_col();++j)
				std::cout<<deltaOmega(i,j)<<" ";
			std::cout<<"\n";
		}

		for (SizeType j=0;j<gckfsc.n_col();++j) {
			integral[j]=integral[j]*(params_.omegaStep);
			std::cerr<<"integral["<<j<<"]="<<integral[j]<<"\n";
		}
	}

	void getGammaOmega(MatrixType& gammaOmega,
	                   const MatrixType& deltaOmega,
	                   PsimagLite::FreqEnum freqEnum)
	{
		if (freqEnum == PsimagLite::FREQ_MATSUBARA) {
			gammaOmega = deltaOmega;
			return;
		}

		hilbertTransfFromReal(gammaOmega,deltaOmega);

		std::cout<<"#GAMMA\n";
		for (SizeType i=0;i<gammaOmega.n_row();++i) {
			RealType wn = latticeFunctions_.matsubara(i);
			std::cout<<wn<<" ";
			for (SizeType j=0;j<gammaOmega.n_col();++j)
				std::cout<<gammaOmega(i,j)<<" ";
			std::cout<<"\n";
		}
	}

	void hilbertTransfFromReal(MatrixType& imagOmega,
	                           const MatrixType& realOmega)
	{
		RealType factor = -params_.omegaStep/M_PI;

		for (SizeType i=0;i<params_.numberOfMatsubaras;++i) {
			RealType wn = latticeFunctions_.matsubara(i);

			for (SizeType k=0;k<realOmega.n_col();++k) {
				ComplexType sum = 0.0;
				for (SizeType j=0;j<realOmega.n_row();++j) {
					RealType rOmega = params_.omegaBegin + params_.omegaStep * j;
					sum += factor * std::imag(realOmega(j,k))/ComplexType(-rOmega,wn);
				}

				imagOmega(i,k) = sum;
			}
		}
	}

	void makeG0(MatrixType& G0inverse,
	            const MatrixType& gammakomega,
	            const VectorRealType& epsbar,
	            PsimagLite::FreqEnum freqEnum)
	{
		assert(freqEnum == PsimagLite::FREQ_MATSUBARA);
		bool lanczosReal = adjustments_.isOption("lanczosreal");

		std::cout<<"#G0\n";
		for (SizeType i = 0;i < gammakomega.n_row(); ++i) {
			ComplexType omega = latticeFunctions_.omegaValue(i,(lanczosReal) ?
			                    PsimagLite::FREQ_REAL : PsimagLite::FREQ_MATSUBARA);
			for (SizeType j = 0;j < epsbar.size(); ++j) {
				// j = clusterK + orb1*largeKs + orb2*largeKs*orbitals
				SizeType clusterK = j % params_.largeKs;
				SizeType tmp = static_cast<SizeType>(j/params_.largeKs);
				SizeType orb1 = static_cast<SizeType>(tmp / params_.orbitals);
				SizeType orb2 = tmp % params_.orbitals;
				if (orb1 != orb2) continue;
				SizeType jj = clusterK + orb1 * params_.largeKs;
				G0inverse(i,jj) = omega + params_.mu -epsbar[j] - gammakomega(i,jj);
				ComplexType tmp2 = 1.0/G0inverse(i,jj);
				std::cout<<tmp2<<" ";
			}

			std::cout<<"\n";
		}
	}

	void ft4(MatrixType& gfdest,
	         const MatrixType& gfsource,
	         const MatrixType& ftCoeffs,
	         SizeType Nc) const
	{
		for (SizeType orb = 0; orb < params_.orbitals; ++orb)
			ft4(gfdest,gfsource,ftCoeffs,Nc,orb);
	}

	void ft4(MatrixType& gfdest,
	         const MatrixType& gfsource,
	         const MatrixType& ftCoeffs,
	         SizeType Nc,
	         SizeType orb) const
	{
		ComplexType csum = 0.0;

		for (SizeType ik=0;ik<Nc;ik++) {
			for (SizeType n=0;n<gfsource.n_row();n++) {
				csum=0.0;
				for (SizeType ic=0;ic<Nc;ic++)
					for (SizeType jc=0;jc<Nc;jc++)
						csum+= gfsource(n,ic+jc*Nc+orb*Nc*Nc)*
						      conj(ftCoeffs(ik,ic))*ftCoeffs(ik,jc);
				gfdest(n,ik+orb*Nc)=csum; //! no factor of 1/Nc
			}
		}
	}

	void disableMpiInDmrgPlusPlus()
	{
		typename PsimagLite::Vector<PsimagLite::String>::Type v;

		v.push_back("BlockMatrix");
		v.push_back("HamiltonianConnection");
		v.push_back("KronConnections");
		v.push_back("Operators");
		v.push_back("ParallelDensityMatrix");
		v.push_back("ParallelTriDiag");
		v.push_back("ParallelWft");

		for (SizeType i = 0; i < v.size(); ++i)
			PsimagLite::Concurrency::mpiDisable(v[i]);
	}

	void getGammaKOmega(MatrixType& gammaFreq,
	                    const MatrixType& p,
	                    PsimagLite::FreqEnum freqEnum)
	{
		for (SizeType k=0;k<gammaFreq.n_col();++k) {
			for (SizeType j=0;j<latticeFunctions_.omegaSize(freqEnum);++j) {
				ComplexType z = latticeFunctions_.omegaValue(j,freqEnum);
				gammaFreq(j,k)=andersonG0(p,k,z);
			}
		}

		std::cout<<"#gammaFreq\n";
		std::cout<<gammaFreq;
	}

	ComplexType andersonG0(const MatrixType& p, SizeType k,ComplexType z) const
	{
		SizeType n=static_cast<SizeType>(p.n_row()/2);
		ComplexType ctmp=0.0;

		for (SizeType i=0;i<n;++i)
			ctmp += std::conj(p(i,k))*p(i,k)/(-p(i+n,k)+z);

		return ctmp;
	}

	RealType onSitePotential(SizeType K, SizeType gamma1, SizeType gamma2) const
	{
		if (gamma1 != gamma2) return 0;
		SizeType index = K + gamma1*geometry_.numberOfSites();
		assert(index< params_.potentialV.size());
		return params_.potentialV[index];
	}

	const ParametersType& params_;
	typename InputNgType::Readable& io_;
	GeometryType geometry_;
	DispersionType dispersion_;
	MatrixType fTCoefsR2K_;
	LatticeFunctionsType latticeFunctions_;
	AdjustmentsType adjustments_;
}; // class DcaLoop

} // namespace OpenDca

#endif // DCA_LOOP_H

