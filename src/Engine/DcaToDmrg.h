#ifndef DCA_TO_DMRG_H
#define DCA_TO_DMRG_H
#include "Vector.h"
#include "String.h"
#include "Matrix.h"

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

template<typename ParametersType,typename GeometryType,typename InputNgType_>
class DcaToDmrg {

	typedef typename ParametersType::RealType RealType_;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::complex<RealType_> ComplexType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef PsimagLite::Matrix<RealType_> MatrixRealType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef HubbardParams<RealType_> HubbardParamsType;

public:

	typedef InputNgType_ InputNgType;
	typedef RealType_ RealType;

	DcaToDmrg(const ParametersType& params,
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
	   connectorsCounter_(0)
	{
		SizeType Nc=params_.largeKs;
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		SizeType total=Nc*(1+nBath);

		// FIXME: CHECK THAT tcluster and tBath have no diagonals

		// Assign "hoppings" and "hubbards"
		//hubbardParams_.exchangeJ.resize(total,total);
		hubbardParams_.potentialV.resize(2*total,0.0);
		hubbardParams_.hoppings.resize(total,total);
		SizeType r = 0;
		SizeType r2 = 0;
		SizeType alpha = 0;
		SizeType alpha2 = 0;

		for (SizeType i=0;i<total;++i) {
			SizeType ind = dcaIndexToDmrgIndex(i);
			for (SizeType j=0;j<total;++j) {
				SizeType jnd = dcaIndexToDmrgIndex(j);
				//hubbardParams_.exchangeJ(i,j)=0.0;
				hubbardParams_.hoppings(ind,jnd)=0.0;
				if (i<params_.largeKs && j<params_.largeKs) { //we're in the cluster
					hubbardParams_.hoppings(ind,jnd)=std::real(tCluster(i,j));
				} else if (i<params_.largeKs && j>=params_.largeKs) { // we're inter cluster
					getBathPoint(r,alpha,j,Nc,nBath);
					hubbardParams_.hoppings(ind,jnd)=std::real(tBathCluster(alpha,r+i*Nc));
				} else if (i>=params_.largeKs && j<params_.largeKs) { // we're inter cluster
					getBathPoint(r,alpha,i,Nc,nBath);
					hubbardParams_.hoppings(ind,jnd)=std::real(tBathCluster(alpha,j+r*Nc));
				} else if (i>=params_.largeKs && j>=params_.largeKs) { //we're in the bath
					getBathPoint(r,alpha,i,Nc,nBath);
					getBathPoint(r2,alpha2,j,Nc,nBath);
					if (alpha==alpha2 && r!=r2)
						hubbardParams_.hoppings(ind,jnd)=std::real(lambda[r+r2*Nc+alpha*Nc*Nc]);
					if (alpha==alpha2 && r==r2)
						hubbardParams_.potentialV[ind]=hubbardParams_.potentialV[ind+total]
						                              =std::real(lambda[r+r*Nc+alpha*Nc*Nc]);
				}
			}
		}

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

	void read(VectorSizeType& x,PsimagLite::String label)
	{
		if (label == "TSPSites" || label == "TSPLoops") {
			io_.read(x,label);
		} else {
			unimplemented("read",label);
		}
	}

	template<typename T>
	void readline(T& x,PsimagLite::String label)
	{
		if (label == "TotalNumberOfSites=") {
			x = hubbardParams_.hoppings.n_row();
		} else if (label == "NumberOfTerms=") {
			x = geometry_.terms();
		} else if (label == "DegreesOfFreedom=" ||
		           label == "LadderLeg=" ||
		           label == "FERMIONSIGN=" ||
		           label == "DynamicDmrgType=") {
			io_.readline(x,label);
		} else if (label == "BathSitesPerSite=") {
			x = params_.nofPointsInBathPerClusterPoint;
		} else {
			unimplemented("readline",label);
		}
	}

	void readline(RealType& x,PsimagLite::String label)
	{
		io_.readline(x,label);
	}

	void read(VectorRealType& v, PsimagLite::String label)
	{
		if (label == "hubbardU") {
			io_.read(v,label);
			return;
		}

		if (label == "potentialV") {
			v = hubbardParams_.potentialV;
			return;
		}

		if (label != "Connectors")
			unimplemented("read",label);

		std::cerr<<"DcaToDmrg::read hardcoded for 2x2 with 2 bath sites\n";
		const MatrixRealType& hoppings = hubbardParams_.hoppings;

		if (connectorsCounter_ == 0) {
			v.resize(2);
			v[0] = hoppings(4,6);
			v[1] = hoppings(5,7);
		} else if (connectorsCounter_ == 1) {
			v.resize(2);
			v[0] = hoppings(4,5);
			v[1] = hoppings(6,7);
		} else if (connectorsCounter_ == 2) {
			SizeType total = hoppings.n_row();
			v.resize(8);
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

		connectorsCounter_++;
		if (connectorsCounter_ > 2) connectorsCounter_ = 0;
	}

	void readline(PsimagLite::String& x,PsimagLite::String label)
	{
		if (label == "GeometryKind=") {
			x = geometry_.label(lastTermSeen_);
		} else if (label == "GeometryOptions=" ||
		           label == "TSPProductOrSum=" ||
		           label == "TSPOperator=" ||
		           label == "CorrectionVectorAlgorithm=") {
			io_.readline(x,label);
		} else {
			unimplemented("readline",label);
		}
	}

	template<typename T>
	void readMatrix(T& t,PsimagLite::String label)
	{
		if (label == "RAW_MATRIX") {
			io_.readMatrix(t,label);
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
		SizeType Nc = params_.largeKs;
		SizeType NcOver2 = static_cast<SizeType>(Nc/2);
		SizeType NcLy = 2;
		SizeType nBath = params_.nofPointsInBathPerClusterPoint;
		SizeType r = 0;

		if (i < Nc) {
			r = clusterIndexDcaToDmrg(i,NcLy);
			return r + NcOver2*nBath;
		}

		SizeType alpha = 0;
		getBathPoint(r,alpha,i,Nc,nBath);
		SizeType r2 = clusterIndexDcaToDmrg(r,NcLy);

		if (r2 < NcOver2)
			return r2 + alpha*nBath;

		return r2 -NcOver2 + alpha*nBath + Nc + NcOver2*nBath;
	}

private:

	SizeType clusterIndexDcaToDmrg(SizeType ind,SizeType ly) const
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

	void unimplemented(PsimagLite::String functionName,PsimagLite::String label)
	{
		PsimagLite::String str = "unimplemented " + functionName + "(...)";
		str += " for label " + label + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	// Formula is i=Nc + Nb*r + alpha, given Nc, Nb and i determines alpha and r
	void getBathPoint(SizeType& r,SizeType& alpha,SizeType i,SizeType Nc,SizeType nBath) const
	{
		assert(i >= Nc);
		i-=Nc;
		r = static_cast<SizeType>(i/nBath);
		alpha = i-r*nBath;
	}

	bool isInBath(SizeType ind) const
	{
		bool b = (ind != 4 && ind != 5);
		b &= (ind != 6 && ind != 7);
		return b;
	}

	const ParametersType& params_;
	const MatrixType& tCluster_;
	const MatrixType& tBathCluster_;
	HubbardParamsType hubbardParams_;
	const GeometryType& geometry_;
	typename InputNgType::Readable& io_;
	SizeType lastTermSeen_;
	SizeType connectorsCounter_;
};

}

#endif

