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
#ifndef DCA_TO_DIAG_H
#define DCA_TO_DIAG_H
#include "Vector.h"
#include "String.h"
#include "Matrix.h"
#include "DcaLoopGlobals.h"

namespace OpenDca {

template<typename RealType>
struct HubbardParams {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;

	MatrixRealType hoppings;
	VectorRealType potentialV;
};

template<typename RealType>
std::ostream& operator<<(std::ostream& os,const HubbardParams<RealType>& hp)
{
	std::cout<<"hoppings\n";
	std::cout<<hp.hoppings;
	std::cout<<"potentialV\n";
	std::cout<<hp.potentialV;
	return os;
}

template<typename ParametersType_,typename GeometryType,typename InputNgType_>
class DcaToDiag {

	typedef typename ParametersType_::RealType RealType_;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::complex<RealType_> ComplexType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef PsimagLite::Matrix<RealType_> MatrixRealType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef HubbardParams<RealType_> HubbardParamsType;

public:

	typedef ParametersType_ ParametersType;
	typedef InputNgType_ InputNgType;
	typedef RealType_ RealType;

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		return 0;
	}

	DcaToDiag(const ParametersType& params,
	          const MatrixType& tCluster,
	          const MatrixType& tBathCluster,
	          const VectorType& lambda,
	          const GeometryType& geometry,
	          typename InputNgType::Readable& io)
	 : params_(params),
	   tCluster_(tCluster),
	   tBathCluster_(tBathCluster),
	   geometry_(geometry),
	   io_(io),
	   lastTermSeen_(0),
	   connectorsCounter_(0),
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

		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		SizeType total=params_.largeKs*(1+nBath)*params_.orbitals;

		// FIXME: CHECK THAT tcluster and tBath have no diagonals

		// Assign "hoppings" and "hubbards"
		//hubbardParams_.exchangeJ.resize(total,total);
		hubbardParams_.potentialV.resize(2*total,0.0);
		hubbardParams_.hoppings.resize(total,total);
		SizeType r = 0;
		SizeType r2 = 0;
		SizeType alpha = 0;
		SizeType alpha2 = 0;

		originalV_.resize(2*total);
		io_.read(originalV_,"potentialV");

		SizeType clusterSize = params_.largeKs * params_.orbitals;
		SizeType nBathOrbitals = nBath * params_.orbitals;
		RealType mu = params_.mu;

		for (SizeType i=0;i<total;++i) {
			SizeType ind = dcaIndexToDmrgIndex(i);
			for (SizeType j=0;j<total;++j) {
				SizeType jnd = dcaIndexToDmrgIndex(j);
				//hubbardParams_.exchangeJ(i,j)=0.0;
				hubbardParams_.hoppings(ind,jnd)=0.0;
				if (i<clusterSize && j<clusterSize) { //we're in the cluster
					hubbardParams_.hoppings(ind,jnd)=std::real(tCluster(i,j));
					if (ind == jnd)
						hubbardParams_.potentialV[ind] =
						  hubbardParams_.potentialV[ind + total] = originalV_[ind] - mu;
				} else if (i<clusterSize && j>=clusterSize) { // we're inter cluster
					getBathPoint(r,alpha,j,clusterSize,nBathOrbitals);
					hubbardParams_.hoppings(ind,jnd)=tBathClusterCorrected(alpha,r,i);
				} else if (i>=clusterSize && j<clusterSize) { // we're inter cluster
					getBathPoint(r,alpha,i,clusterSize,nBathOrbitals);
					hubbardParams_.hoppings(ind,jnd)=tBathClusterCorrected(alpha,r,j);
				} else if (i>=clusterSize && j>=clusterSize) { //we're in the bath
					getBathPoint(r,alpha,i,clusterSize,nBathOrbitals);
					getBathPoint(r2,alpha2,j,clusterSize,nBathOrbitals);

					if (alpha==alpha2 && r!=r2) {
						hubbardParams_.hoppings(ind,jnd) = 0.0;
						  //lambdaCorrected(lambda,r,r2,alpha);
					}

					if (alpha==alpha2 && r==r2) {
						hubbardParams_.potentialV[ind] =
						  hubbardParams_.potentialV[ind+total]
						  = lambdaCorrected(lambda,r,r2,alpha);
					}
				}
			}
		}

		for (SizeType k = 0; k < hubbardParams_.hoppings.n_row(); ++k)
			hubbardParams_.hoppings(k,k) = 0.0;

		for (SizeType k = 0; k < params_.largeKs; ++k) {
			VectorRealType kvector(geometry_.dimension());
			geometry_.index2Rvector(k,kvector);
			std::cout<<"index2Kvector for k = "<<k<<"\n";
			for (SizeType i = 0; i < kvector.size(); ++i)
				std::cout<<kvector[i]<<" ";
			std::cout<<"\n";
		}

		std::cout<<hubbardParams_;
	}

	void updatePotentialV()
	{
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		SizeType total=params_.largeKs*(1+nBath)*params_.orbitals;
		SizeType clusterSize = params_.largeKs * params_.orbitals;

		for (SizeType i=0;i<total;++i) {
			SizeType ind = dcaIndexToDmrgIndex(i);
			for (SizeType j=0;j<total;++j) {
				SizeType jnd = dcaIndexToDmrgIndex(j);
				if (i>=clusterSize || j>=clusterSize)
					continue; //we're NOT in the cluster
				if (ind != jnd) continue; // we're NOT on-site
				hubbardParams_.potentialV[ind] =
				  hubbardParams_.potentialV[ind + total] = originalV_[ind] - params_.mu;
			}
		}
	}

	void read(VectorSizeType& x,PsimagLite::String label)
	{
		if (label == "TSPSites" || label == "TSPLoops" || label == "COOKED_EXTRA") {
			io_.read(x,label);
		} else {
			unimplemented("read",label);
		}
	}

	template<typename T>
	void readline(T& x,PsimagLite::String label)
	{
		if (label == "TotalNumberOfSites=") {
			x = static_cast<T>(hubbardParams_.hoppings.n_row()/params_.orbitals);
		} else if (label == "NumberOfTerms=") {
			x = geometry_.terms();
		} else if (label == "DegreesOfFreedom=" ||
		           label == "LadderLeg=" ||
		           label == "FERMIONSIGN=" ||
		           label == "DynamicDmrgType=" ||
		           label == "Orbitals=" ||
		           label == "FeAsMode=" ||
		           label == "LanczosStepsForEnergyConvergence=" ||
		           label == "LanczosSaveLanczosVectors=" ||
		           label == "LanczosSteps=" ||
		           label == "DynamicDmrgStepsForEnergyConvergence=" ||
		           label == "DynamicDmrgSaveLanczosVectors=" ||
		           label == "DynamicDmrgSteps=") {
			io_.readline(x,label);
		} else if (label == "TargetElectronsUp=") {
			x = electronsUp_;
		} else if (label == "TargetElectronsDown=") {
			x = electronsDown_;
		} else if (label == "BathSitesPerSite=") {
			x = params_.nofPointsInBathPerClusterPoint;
		} else {
			std::cerr<<"WARNING: readline unrecognized "<<label<<"\n";
			io_.readline(x,label);
		}
	}

	void readline(RealType& x,PsimagLite::String label)
	{
		io_.readline(x,label);
	}

	void read(VectorRealType& v, PsimagLite::String label)
	{
		if (label == "hubbardU" || label == "FiniteLoops") {
			io_.read(v,label);
			return;
		}

		if (label == "potentialV") {
			v = hubbardParams_.potentialV;
			return;
		}

		if (label != "Connectors")
			unimplemented("read",label);

		setHoppings(v, hubbardParams_.hoppings);

		connectorsCounter_++;
		if (connectorsCounter_ >= geometry_.directions(0))
			connectorsCounter_ = 0;
	}

	void readline(PsimagLite::String& x,
	              PsimagLite::String label,
	              bool clean = true)
	{
		if (label == "GeometryKind=") {
			x = geometry_.label(lastTermSeen_);
		} else {
			io_.readline(x,label);
		}
	}

	template<typename T>
	void readMatrix(T& t,PsimagLite::String label)
	{
		SizeType nBath = params_.nofPointsInBathPerClusterPoint;
		if (label == "RAW_MATRIX") {
			io_.readMatrix(t,label);
		} else if (label == "Connectors") {
			setHoppingsDmft(t,hubbardParams_.hoppings);
			connectorsCounter_++;
			if (connectorsCounter_ >= nBath)
				connectorsCounter_ = 0;
		} else {
			unimplemented("readMatrix",label);
		}
	}

	void readKnownSize(VectorSizeType& x,PsimagLite::String label)
	{
		if (label == "JMVALUES") {
			io_.readKnownSize(x,label);
		} else {
			unimplemented("readKnownSize",label);
		}
	}

	SizeType dcaIndexToDmrgIndex(SizeType i) const
	{
		if (params_.largeKs == 1) {
			SizeType orb = i % params_.orbitals;
			SizeType site = static_cast<SizeType>(i/params_.orbitals);
			return site + orb*geometry_.numberOfSites();
		} else {
			if (params_.orbitals > 1) {
				PsimagLite::String str("Nc>1 and orbitals>1 not supported\n");
				throw PsimagLite::RuntimeError(str);
			}
		}

		SizeType Nc = params_.largeKs;
		SizeType NcOver2 = static_cast<SizeType>(Nc/2);
		SizeType NcLy = 2;
		SizeType nBath = params_.nofPointsInBathPerClusterPoint;
		SizeType r = 0;

		if (i < Nc) {
			r = clusterIndexDcaToDiag(i,NcLy);
			return r + NcOver2*nBath;
		}

		SizeType alpha = 0;
		getBathPoint(r,alpha,i,Nc,nBath);
		SizeType r2 = clusterIndexDcaToDiag(r,NcLy);

		if (r2 < NcOver2)
			return r2 + alpha*nBath;

		return r2 -NcOver2 + alpha*nBath + Nc + NcOver2*nBath;
	}

	SizeType muFeatureSize() const
	{
		if (params_.dcaOptions.find("nomufeature") != PsimagLite::String::npos)
			return 0;

		SizeType onp1 = muFeatureOffset_.size();
		SizeType onp2 = onp1 + 1;
		return static_cast<SizeType>(onp1*onp2*0.5);
	}

	void muFeatureSet(SizeType ind)
	{
		SizeType onp1 = muFeatureOffset_.size();
		for (SizeType i = 0; i < onp1; ++i) {
			if (ind < muFeatureOffset_[i]) return muFeatureSet(ind,i);
		}
	}

	SizeType electrons(SpinEnum spin) const
	{
		return (spin == SPIN_UP) ? electronsUp_ : electronsDown_;
	}

private:

	void muFeatureSet(SizeType ind, SizeType o)
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

	RealType tBathClusterCorrected(SizeType alpha, SizeType r, SizeType ind) const
	{
		SizeType bathOrbital = alpha % params_.orbitals;
		SizeType clusterOrbital = ind %  params_.orbitals;
		if (bathOrbital != clusterOrbital) return 0.0;

		SizeType csno = static_cast<SizeType>(ind/params_.orbitals);
		SizeType bathSite = static_cast<SizeType>(alpha/params_.orbitals);
		SizeType largeKs2 = params_.largeKs * params_.largeKs;
		assert(csno < params_.largeKs);
		SizeType index = csno + csno*params_.largeKs + clusterOrbital*largeKs2;
		return std::real(tBathCluster_(bathSite,index));
	}

	RealType lambdaCorrected(const VectorType& lambda,
	                         SizeType ind,
	                         SizeType jnd,
	                         SizeType alpha) const
	{
		SizeType csno1 = static_cast<SizeType>(ind/params_.orbitals);
		SizeType csno2 = static_cast<SizeType>(jnd/params_.orbitals);
		SizeType largeKs2 = params_.largeKs * params_.largeKs;
		SizeType nBath = static_cast<SizeType>(lambda.size()/params_.orbitals);
		SizeType alphaCorrected = orbitalBathToBathOrbital(alpha, nBath);
		SizeType index = csno1 + csno2*params_.largeKs + alphaCorrected* largeKs2;
		assert(index < lambda.size());
		return std::real(lambda[index]);
	}

	SizeType clusterIndexDcaToDiag(SizeType ind,SizeType ly) const
	{
		SizeType dim = geometry_.dimension();
		VectorRealType rvector(dim);
		geometry_.index2Rvector(ind,rvector);
		SizeType x = static_cast<SizeType>(rvector[0]);
		if (dim == 1)
			return x;

		assert(dim == 2);
		return static_cast<SizeType>(rvector[1]) + x*ly;

	}

	// Formula is i=clusterSize + Nb*r + alpha, given Nc, Nb and i determines alpha and r
	void getBathPoint(SizeType& r,
	                  SizeType& alpha,
	                  SizeType i,
	                  SizeType clusterSize,
	                  SizeType nBath) const
	{
		assert(i >= clusterSize);
		i-=clusterSize;
		r = static_cast<SizeType>(i/nBath);
		alpha = i-r*nBath;
	}

	bool isInBath(SizeType ind) const
	{
		bool b = (ind != 4 && ind != 5);
		b &= (ind != 6 && ind != 7);
		return b;
	}

	SizeType orbitalBathToBathOrbital(SizeType alpha, SizeType nBath) const
	{
		div_t q = div(alpha,params_.orbitals);
		return q.rem * nBath + q.quot;
	}

	void setHoppings(VectorRealType& v, const MatrixRealType& hoppings) const
	{
		if (params_.largeKs == 1) {
			PsimagLite::String str("setHoppings\n");
			throw PsimagLite::RuntimeError(str);
		} else {
			if (params_.orbitals > 1) {
				PsimagLite::String str("Nc>1 and orbitals>1 not supported\n");
				throw PsimagLite::RuntimeError(str);
			}

			if (geometry_.label(0) != "ladderbath") {
				PsimagLite::String str("Only ladderbath supported for Nc>1\n");
				throw PsimagLite::RuntimeError(str);
			}
		}

		std::cerr<<"setHoppings: ladderLeg HARDWIRED to 2\n";
		SizeType ladderLeg = 2;
		SizeType lx = static_cast<SizeType>(params_.largeKs/ladderLeg);
		SizeType nBath = params_.nofPointsInBathPerClusterPoint;
		SizeType start = nBath * static_cast<SizeType>(params_.largeKs/2);;

		if (connectorsCounter_ == 0) {
			SizeType cx = 2*(lx-1);
			v.resize(cx);
			for (SizeType i = 0; i < cx; i+=2) {
				v[i] = hoppings(start,start+2);
				v[i+1] = hoppings(start+1,start+3);
				start+=2;
			}
		} else if (connectorsCounter_ == 1) {
			SizeType cx = 2*(lx-1);
			v.resize(cx);
			for (SizeType i = 0; i < cx; ++i) {
				v[i] = hoppings(start,start+1);
				start+=2;
			}
		} else if (connectorsCounter_ == 2) {
			SizeType total = hoppings.n_row();
			v.resize(params_.largeKs*nBath);
			for (SizeType i = 0; i < total; ++i) {
				if (!isInBath(i)) continue;
				for (SizeType j = 0; j < total; ++j) {
					if (!isInBath(j)) continue;
					SizeType handle = geometry_.handle(lastTermSeen_,i,j);
					assert(handle < v.size());
					v[handle] = hoppings(i,j);
				}
			}
		}
	}

	template<typename T>
	void setHoppingsDmft(T& m, const MatrixRealType& hoppings) const
	{
		throw PsimagLite::RuntimeError("setHoppingsDmft error\n");
	}

	void setHoppingsDmft(MatrixRealType& m, const MatrixRealType& hoppings) const
	{
		assert(geometry_.label(0) == "star");
		assert(params_.largeKs == 1);

		SizeType i = connectorsCounter_ + 1;
		m.resize(params_.orbitals,params_.orbitals);
		m.setTo(0.0);
		for (SizeType orb = 0; orb < params_.orbitals; ++orb)
			m(orb,orb) = hubbardParams_.hoppings(0+orb*geometry_.numberOfSites(),
			                                     i+orb*geometry_.numberOfSites());
	}

	void unimplemented(PsimagLite::String functionName,PsimagLite::String label)
	{
		PsimagLite::String str = "unimplemented " + functionName + "(...)";
		str += " for label " + label + "\n";
		std::cerr<<"ERROR: "<<str;
		throw PsimagLite::RuntimeError(str);
	}

	const ParametersType& params_;
	const MatrixType& tCluster_;
	const MatrixType& tBathCluster_;
	HubbardParamsType hubbardParams_;
	const GeometryType& geometry_;
	typename InputNgType::Readable io_;
	SizeType lastTermSeen_;
	SizeType connectorsCounter_;
	SizeType electronsUp_;
	SizeType electronsDown_;
	VectorRealType originalV_;
	VectorSizeType muFeatureOffset_;
}; // class DcaToDiag

}

#endif

