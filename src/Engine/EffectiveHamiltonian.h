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
#ifndef EFFECTIVE_HAMILTONIAN_H
#define EFFECTIVE_HAMILTONIAN_H

#include "Vector.h"
#include "NoPthreads.h"
#include "Parallelizer.h"
#include "DcaToDiag.h"
#include "InputNg.h"
#include "AndersonFit.h"
#include "FreqEnum.h"
#include "RootFindingBisection.h"
#include "ClusterFunctions.h"
#include "Adjustments.h"

namespace OpenDca {

template<typename ParametersType, typename GeometryType, typename InputNgType>
class EffectiveHamiltonian {

	typedef typename ParametersType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef Dmrg::Basis<SparseMatrixType> BasisType;
	typedef Dmrg::Operators<BasisType> OperatorsType;
	typedef Dmrg::BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef Dmrg::LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef Dmrg::ModelHelperLocal<LeftRightSuperType> ModelHelperType;
	typedef DcaToDiag<ParametersType,GeometryType,InputNgType> DcaToDiagType_;
	typedef ClusterFunctions<DcaToDiagType_> ClusterFunctionsType;
	typedef typename ClusterFunctionsType::MatrixType MatrixType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef AndersonFit<VectorType,ParametersType> AndersonFitType;
	typedef Adjustments<ClusterFunctionsType> AdjustmentsType;

public:

	typedef DcaToDiagType_ DcaToDiagType;

	EffectiveHamiltonian(const ParametersType& params,
	                                 const GeometryType& geometry,
	                                 typename InputNgType::Readable& io)
	: params_(params),
	  geometry_(geometry),
	  io_(io),
	  p_(2*params_.nofPointsInBathPerClusterPoint,params_.largeKs*params_.orbitals),
	  gammaRealFreq_(params_.omegas,params_.largeKs*params_.orbitals),
	  garbage_(0),
	  clusterFunctions_(params_,io),
	  adjustments_(clusterFunctions_,params_)
	{}

	~EffectiveHamiltonian()
	{
		delete garbage_;
	}

	// Get LanczosParams from Gf
	void build(MatrixType& gf,
	           const VectorRealType& ekbar,
	           PsimagLite::FreqEnum)
	{
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		SizeType largeKs = params_.largeKs;

		VectorRealType p(2*nBath);
		VectorType gfTmp(gf.n_row());
		AndersonFitType andersonFit(params_);

		for (SizeType site = 0; site < largeKs; ++site) {
			for (SizeType orb = 0; orb < params_.orbitals; ++orb) {
				SizeType jj = site + orb * largeKs + orb * largeKs * params_.orbitals;
				std::cout<<"#GFTMP"<<site<<" "<<orb<<"\n";
				for (SizeType i=0;i<gf.n_row();++i) {
					gfTmp[i]=gf(i,jj);
					std::cout<<gfTmp[i]<<"\n";
				}

				andersonFit.fit(p,gfTmp);
				std::cout<<"#PANDERSON"<<site<<" "<<orb<<"\n";
				for (SizeType i=0;i<p.size();++i) std::cout<<i<<" "<<p[i]<<"\n";

				saveAndersonParameters(p,site + orb*largeKs);
			}
		}

		inversionSymmetry(p_,0,nBath);
		inversionSymmetry(p_,nBath,p_.n_row());

		bool allbathSitesEqual = adjustments_.isOption("AllBathSitesEqual");
		if (allbathSitesEqual) bathSitesSymmetry(p_);

		bool allOrbitalsEqual = adjustments_.isOption("AllOrbitalsEqual");
		if (allOrbitalsEqual) orbitalsSymmetry(p_);

		std::cout<<"#REALGAMMA\n";
		std::cout<<gammaRealFreq_.n_row()<<" "<<gammaRealFreq_.n_col()<<"\n";
		for (SizeType i=0;i<gammaRealFreq_.n_row();++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j=0;j<gammaRealFreq_.n_col();++j)
				std::cout<<gammaRealFreq_(i,j)<<" ";
			std::cout<<"\n";
		}

		std::cout<<"#\n";

		std::cout<<"andersonVp and andersonEp\n";
		std::cout<<p_;

		// Calculate "hoppings"
		garbage_ = makeHubbardParams(ekbar,io_);

		clusterFunctions_.build(garbage_);
	}

	void solve(MatrixType& gfCluster)
	{
		if (!garbage_) {
			PsimagLite::String str("solve must call build first\n");
			throw PsimagLite::RuntimeError("EffectiveHamiltonian::" + str);
		}

		bool adjustMuCluster = adjustments_.isOption("adjustmucluster");

		if (adjustMuCluster) adjustments_.adjChemPot();

		clusterFunctions_.sweepParticleSectors(true);

		clusterFunctions_.findGf(gfCluster);
	}

	const MatrixType& andersonParameters() const
	{
		return p_;
	}

private:

	void saveAndersonParameters(const VectorRealType&src,SizeType k)
	{
		SizeType n=src.size();
		assert(p_.n_row() == n);
		assert(k < p_.n_col());
		for (SizeType i=0;i<n;++i) {
			p_(i,k) = src[i];
		}
	}

	DcaToDiagType* makeHubbardParams(const VectorRealType& ekbar,
	                                 typename InputNgType::Readable& io)
	{
		SizeType Nc=params_.largeKs;
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		MatrixType tCluster(Nc*params_.orbitals,Nc*params_.orbitals);
		MatrixType tBathCluster(nBath,Nc*Nc*params_.orbitals);
		VectorType vbar(Nc*Nc*nBath*params_.orbitals);

		ft(tCluster,ekbar);
		ft5(tBathCluster,p_);
		calcBathVbar(vbar,p_);

		std::cout<<"tCluster\n";
		std::cout<<tCluster;

		std::cout<<"tBathCluster\n";
		std::cout<<tBathCluster;

		return new DcaToDiagType(params_,tCluster,tBathCluster,vbar,geometry_,io);
	}

	void ft(MatrixType& dest,const VectorRealType& src) const
	{
		for (SizeType gamma = 0; gamma < params_.orbitals; ++gamma)
			ft(dest,src,gamma);
	}

	void ft(MatrixType& dest,const VectorRealType& src,SizeType gamma) const
	{
		SizeType largeKs = params_.largeKs;

		for (SizeType ii=0;ii<largeKs;++ii) {
			SizeType i = ii + gamma*largeKs;
			for (SizeType jj=0;jj<largeKs;++jj) {
				SizeType j = jj + gamma*largeKs;
				dest(i,j)=0.0;
				for (SizeType k=0;k<largeKs;++k) {
					RealType tmp = krProduct(k,i) - krProduct(k,j);
					dest(i,j) += ComplexType(cos(tmp),sin(tmp))*src[k];
				}

				dest(i,j) = std::real(dest(i,j))/largeKs;
				//if (std::norm(dest(i,j))<1e-8) dest(i,j)=0.0;
			}
		}
	}

	void ft5(MatrixType& dest,const MatrixType& src) const
	{
		for (SizeType gamma = 0; gamma < params_.orbitals; ++gamma)
			ft5(dest,src,gamma);
	}

	void ft5(MatrixType& dest,const MatrixType& src,SizeType gamma) const
	{
		SizeType Nc = params_.largeKs;
		SizeType nBath = dest.n_row();

		for (SizeType b=0;b<nBath;++b) { // bath
			for (SizeType i=0;i<Nc;++i) { //cluster real space
				for (SizeType j=0;j<Nc;++j) { //cluster real space
					dest(b,i+j*Nc+gamma*Nc*Nc) = 0.0;
					for (SizeType k=0;k<Nc;++k) { //cluster k-space
						RealType tmp = krProduct(k,i,j);
						dest(b,i+j*Nc+gamma*Nc*Nc)+=
						ComplexType(cos(tmp),sin(tmp))*src(b,k+gamma*Nc);
					}

					dest(b,i+j*Nc+gamma*Nc*Nc)=std::real(dest(b,i+j*Nc+gamma*Nc*Nc))/Nc;
				}
			}
		}
	}

	void calcBathVbar(VectorType& vbar,const MatrixType& src) const
	{
		for (SizeType gamma = 0; gamma < params_.orbitals; ++gamma)
			calcBathVbar(vbar,src,gamma);
	}

	// src is (2 * nbath, Nc * orbitals)
	void calcBathVbar(VectorType& vbar,
	                  const MatrixType& src,
	                  SizeType gamma) const
	{
		SizeType Nc=params_.largeKs;
		SizeType nBath = static_cast<SizeType>(0.5*src.n_row());
		assert(vbar.size() == Nc*Nc*nBath*params_.orbitals);

		for (SizeType i=0;i<nBath;i++) { // i is in the bath
			SizeType ii = (i + nBath*gamma)*Nc*Nc;
			for (SizeType j=0;j<Nc;j++) { // j is in the cluster (real space)
				for (SizeType s=0;s<Nc;s++) { // s is in the cluster (real space)
					vbar[j+s*Nc+ii]=0.0;
					for (SizeType k=0;k<Nc;k++) { // k is in the cluster (k-space)
						RealType tmp = krProduct(k,j,s);
						vbar[j+s*Nc+ii]+=
						ComplexType(cos(tmp),sin(tmp))*src(i + nBath,k+gamma*Nc);
					}

					vbar[j+s*Nc+ii] /= Nc;
				}
			}
		}
	}

	RealType krProduct(SizeType k, SizeType r) const
	{
		SizeType ndim=geometry_.dimension();
		VectorRealType rvector(ndim),kvector(ndim);

		geometry_.index2Rvector(r,rvector);
		geometry_.index2Kvector(k,kvector);
		return scalarProduct(kvector,rvector);
	}

	RealType krProduct(SizeType k, SizeType r, SizeType r2) const
	{
		SizeType ndim=geometry_.dimension();
		VectorRealType rvector(ndim),kvector(ndim),rvector2(ndim);

		geometry_.index2Rvector(r,rvector);
		geometry_.index2Rvector(r2,rvector2);
		geometry_.index2Kvector(k,kvector);

		for (SizeType i=0;i<rvector.size();++i)
			rvector[i]-=rvector2[i];

		return scalarProduct(kvector,rvector);
	}

	void inversionSymmetry(MatrixType& p,int start, int end)
	{
		if (params_.largeKs == 1) return;

		int nBath = end - start;
		MatrixType tempMatrix(nBath,p.n_col());

		for (int i=start;i<end;++i) {
			for (SizeType j=0;j<p.n_col();++j) {
				tempMatrix(i-start,j)=0.0;
				for (SizeType op=0;op<geometry_.nGroupK();++op) {
					SizeType k = geometry_.ickequ(j,op);
					tempMatrix(i-start,j) += p(i,k);
				}
				tempMatrix(i-start,j) /= static_cast<RealType>(geometry_.nGroupK());
			}
		}

		for (int i=start;i<end;++i)
			for (SizeType j=0;j<p.n_col();++j)
				p(i,j)=tempMatrix(i-start,j); // andersonp = tempMatrix
	}

	void bathSitesSymmetry(MatrixType& p)
	{
		SizeType nBath = static_cast<SizeType>(0.5*p.n_row());
		for (SizeType bathSite = 1; bathSite < nBath; ++bathSite) {
			for (SizeType j = 0; j < p_.n_col(); ++j) {
				p(bathSite, j) = p(0, j);
				p(bathSite + nBath, j) = p(nBath, j);
			}
		}
	}

	void orbitalsSymmetry(MatrixType& p)
	{
		SizeType norbitals = params_.orbitals;

		for (SizeType orbital = 1; orbital < norbitals; ++orbital)
			for (SizeType j = 0; j < p_.n_row(); ++j)
				for (SizeType k = 0; k < params_.largeKs; ++k)
					p(j, k + orbital*params_.largeKs) = p(j, k);
	}

	const ParametersType& params_;
	const GeometryType& geometry_;
	typename InputNgType::Readable& io_;
	MatrixType p_;
	MatrixType gammaRealFreq_;
	DcaToDiagType* garbage_;
	ClusterFunctionsType clusterFunctions_;
	AdjustmentsType adjustments_;
};

template<typename RealType, typename GeometryType,typename InputNgType>
std::ostream& operator<<(std::ostream& os,
                         EffectiveHamiltonian<RealType,GeometryType,InputNgType>&)
{
	os<<"EffectiveHamiltonian\n";
	return os;
}

} // namespace DynamicClusterApprox

#endif // EFFECTIVE_HAMILTONIAN_H

