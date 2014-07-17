#ifndef EFFECTIVE_HAMILTONIAN_H
#define EFFECTIVE_HAMILTONIAN_H

#include "Vector.h"
#include "NoPthreads.h"
#include "Parallelizer.h"
#include "DcaToDmrg.h"
#include "InputNg.h"
#include "AndersonFit.h"
#include "FreqEnum.h"
#include "RootFindingBisection.h"
#include "ClusterFunctions.h"

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
	typedef DcaToDmrg<ParametersType,GeometryType,InputNgType> DcaToDmrgType_;
	typedef ClusterFunctions<DcaToDmrgType_> ClusterFunctionsType;
	typedef typename ClusterFunctionsType::MatrixType MatrixType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef AndersonFit<VectorType,ParametersType> AndersonFitType;

	typedef std::pair<RealType,RealType> PairRealRealType;

public:

	typedef DcaToDmrgType_ DcaToDmrgType;

	EffectiveHamiltonian(const ParametersType& params,
	                                 const GeometryType& geometry,
	                                 typename InputNgType::Readable& io)
	: params_(params),
	  geometry_(geometry),
	  io_(io),
	  p_(2*params_.nofPointsInBathPerClusterPoint,params_.largeKs*params_.orbitals),
	  gammaRealFreq_(params_.omegas,params_.largeKs*params_.orbitals),
	  garbage_(0),
	  clusterFunctions_(params_,io)
	{}

	~EffectiveHamiltonian()
	{
		delete garbage_;
	}

	// Get LanczosParams from Gf
	void build(MatrixType& gf,
	           const VectorRealType& ekbar,
	           const VectorRealType& integral,
	           PsimagLite::FreqEnum freqEnum)
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

				andersonFit.fit(p,gfTmp,integral[jj]);
				std::cout<<"#PANDERSON"<<site<<" "<<orb<<"\n";
				for (SizeType i=0;i<p.size();++i) std::cout<<i<<" "<<p[i]<<"\n";

				saveAndersonParameters(p,site + orb*largeKs);
			}
		}

		inversionSymmetry(p_,0,nBath);
		inversionSymmetry(p_,nBath,p_.n_row());

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
	}

	void solve(MatrixType& gfCluster)
	{
		if (!garbage_) {
			PsimagLite::String str("solve must call build first\n");
			throw PsimagLite::RuntimeError("EffectiveHamiltonian::" + str);
		}

		//bool adjustMuCluster = (!isOption("noadjustmu") & !isOption("adjustmulattice"));

		//if (adjustMuCluster) adjChemPot();

		clusterFunctions_.sweepParticleSectors(*garbage_);

		clusterFunctions_.findGf(gfCluster,*garbage_);
	}

	const MatrixType& andersonParameters() const
	{
		return p_;
	}

private:

	void adjChemPot() const
	{
		RealType mu = adjChemPot_();
		std::cout<<"Old mu= "<<params_.mu<<" ";
		params_.mu = mu;
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
		typedef PsimagLite::RootFindingBisection<ClusterFunctionsType> RootFindingType;
		RootFindingType  rootFinding(clusterFunctions_,aAndB.first, aAndB.second,1000,tol);
		RealType mu = params_.mu;
		rootFinding(mu);
		return mu;
	}

	PairRealRealType findAandB() const
	{
		for (RealType value = 1.0; value < 100.0; value++) {
			RealType value2 = clusterFunctions_(value) * clusterFunctions_(-value);
			if (value2 < 0) return PairRealRealType(-value,value);
		}

		throw PsimagLite::RuntimeError("RootFinding init failed\n");
	}

	void saveAndersonParameters(const VectorRealType&src,SizeType k)
	{
		SizeType n=src.size();
		assert(p_.n_row() == n);
		assert(k < p_.n_col());
		for (SizeType i=0;i<n;++i) {
			p_(i,k) = src[i];
		}
	}

	DcaToDmrgType* makeHubbardParams(const VectorRealType& ekbar,
	                                 typename InputNgType::Readable& io)
	{
		SizeType Nc=params_.largeKs;
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		MatrixType tCluster(Nc*params_.orbitals,Nc*params_.orbitals);
		MatrixType tBathCluster(nBath,Nc*Nc*params_.orbitals);
		VectorType lambda(Nc*Nc*nBath*params_.orbitals);

		ft(tCluster,ekbar);
		ft5(tBathCluster,p_,0);
		calcBathLambda(lambda,p_,nBath);

		std::cout<<"tCluster\n";
		std::cout<<tCluster;

		std::cout<<"tBathCluster\n";
		std::cout<<tBathCluster;

		return new DcaToDmrgType(params_,tCluster,tBathCluster,lambda,geometry_,io);
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

	void ft5(MatrixType& dest,const MatrixType& src,SizeType offset) const
	{
		for (SizeType gamma = 0; gamma < params_.orbitals; ++gamma)
			ft5(dest,src,offset,gamma);
	}

	void ft5(MatrixType& dest,const MatrixType& src,SizeType offset,SizeType gamma) const
	{
		SizeType Nc = params_.largeKs;

		for (SizeType b=0;b<dest.n_row();++b) { // bath
			for (SizeType i=0;i<Nc;++i) { //cluster real space
				for (SizeType j=0;j<Nc;++j) { //cluster real space
					dest(b,i+j*Nc+gamma*Nc*Nc) = 0.0;
					for (SizeType k=0;k<Nc;++k) { //cluster k-space
						RealType tmp = krProduct(k,i,j);
						dest(b,i+j*Nc+gamma*Nc*Nc)+=
						ComplexType(cos(tmp),sin(tmp))*src(b+offset,k+gamma*Nc);
					}

					dest(b,i+j*Nc+gamma*Nc*Nc)=std::real(dest(b,i+j*Nc+gamma*Nc*Nc))/Nc;
				}
			}
		}
	}

	void calcBathLambda(VectorType& lambda,const MatrixType& src,SizeType offset) const
	{
		for (SizeType gamma = 0; gamma < params_.orbitals; ++gamma)
			calcBathLambda(lambda,src,offset,gamma);
	}

	void calcBathLambda(VectorType& lambda,
	                    const MatrixType& src,
	                    SizeType offset,
	                    SizeType gamma) const
	{
		SizeType Nc=params_.largeKs;
		SizeType nBath = static_cast<SizeType>(0.5*src.n_row());
		assert(lambda.size() == Nc*Nc*nBath*params_.orbitals);

		for (SizeType i=0;i<nBath;i++) { // i is in the bath
			SizeType ii = (i + nBath*gamma)*Nc*Nc;
			for (SizeType j=0;j<Nc;j++) { // j is in the cluster (real space)
				for (SizeType s=0;s<Nc;s++) { // j is in the cluster (real space)
					lambda[j+s*Nc+ii]=0.0;
					for (SizeType k=0;k<Nc;k++) { // k is in the cluster (k-space)
						RealType tmp = krProduct(k,j,s);
						lambda[j+s*Nc+ii]+=
						ComplexType(cos(tmp),sin(tmp))*src(i+offset,k+gamma*Nc);
					}

					lambda[j+s*Nc+ii] /= Nc;
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

	bool isOption(PsimagLite::String what) const
	{
		return (params_.dcaOptions.find(what) != PsimagLite::String::npos);
	}

	const ParametersType& params_;
	const GeometryType& geometry_;
	typename InputNgType::Readable& io_;
	MatrixType p_;
	MatrixType gammaRealFreq_;
	DcaToDmrgType* garbage_;
	ClusterFunctionsType clusterFunctions_;
};

template<typename RealType, typename GeometryType,typename InputNgType>
std::ostream& operator<<(std::ostream& os,
                         EffectiveHamiltonian<RealType,GeometryType,InputNgType>& params)
{
	os<<"EffectiveHamiltonian\n";
	return os;
}

} // namespace DynamicClusterApprox

#endif // EFFECTIVE_HAMILTONIAN_H

