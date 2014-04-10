#ifndef EFFECTIVE_HAMILTONIAN_H
#define EFFECTIVE_HAMILTONIAN_H

#include "Vector.h"
#include "NoPthreads.h"
#include "Parallelizer.h"
#include "DcaParameters.h"
#include "DcaToDmrg.h"
#include "Geometry/Geometry.h"
#include "InputNg.h"
#include "AndersonFit.h"
#include "DcaSolverLanczos.h"
#include "DcaSolverDmrg.h"

namespace OpenDca {

template<typename RealType, typename InputNgType>
class EffectiveHamiltonian {

	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::CrsMatrix<ComplexType> SparseMatrixType;
	typedef Dmrg::Basis<SparseMatrixType> BasisType;
	typedef Dmrg::Operators<BasisType> OperatorsType;
	typedef Dmrg::BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef Dmrg::LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef Dmrg::ModelHelperLocal<LeftRightSuperType> ModelHelperType;
	typedef PsimagLite::Geometry<ComplexType,
	                             typename InputNgType::Readable,
	                             Dmrg::ProgramGlobals> GeometryType_;
	typedef DcaParameters<RealType> ParametersType_;
	typedef DcaToDmrg<ParametersType_,GeometryType_,InputNgType> DcaToDmrgType_;
	typedef PsimagLite::Geometry<RealType,DcaToDmrgType_,Dmrg::ProgramGlobals> VaryingGeometryType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef AndersonFit<VectorType,ParametersType_> AndersonFitType;
	typedef DcaSolverLanczos<DcaToDmrgType_, VaryingGeometryType> DcaSolverLanczosType;
	typedef DcaSolverDmrg<DcaToDmrgType_, VaryingGeometryType>  DcaSolverDmrgType;
	typedef typename DcaSolverDmrgType::DcaSolverBaseType DcaSolverBaseType;
	typedef typename DcaSolverBaseType::VectorSizeType VectorSizeType;
	typedef typename DcaSolverBaseType::PlotParamsType PlotParamsType;

public:

	typedef GeometryType_ GeometryType;
	typedef ParametersType_ ParametersType;
	typedef DcaToDmrgType_ DcaToDmrgType;

	EffectiveHamiltonian(const ParametersType& params,
	                                 const GeometryType& geometry,
	                                 typename InputNgType::Readable& io)
	: params_(params),
	  geometry_(geometry),
	  io_(io),
	  andersonVp_(params_.nofPointsInBathPerClusterPoint,params_.largeKs),
	  andersonEp_(params_.nofPointsInBathPerClusterPoint,params_.largeKs),
	  gammaRealFreq_(params_.omegas,params_.largeKs),
	  garbage_(0)
	{}

	~EffectiveHamiltonian()
	{
		delete garbage_;
	}

	// Get LanczosParams from Gf
	void build(MatrixType& gf,
	           const VectorRealType& ekbar,
	           const VectorRealType& integral)
	{
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;

		VectorRealType p(2*nBath);
		VectorType gfTmp(gf.n_row());
		AndersonFitType andersonFit(params_);

		for (SizeType j=0;j<params_.largeKs;++j) {
			// first find Anderson parameters by fitting procedure (min. of distance)
			std::cout<<"#GFTMP"<<j<<"\n";
			for (SizeType i=0;i<gf.n_row();++i) {
				gfTmp[i]=gf(i,j);
				std::cout<<i<<" "<<std::real(gfTmp[i])<<" "<<std::imag(gfTmp[i])<<"\n";
			}

			andersonFit.fit(p,gfTmp,integral[j]);
			std::cout<<"#PANDERSON"<<j<<"\n";
			for (SizeType i=0;i<p.size();++i) std::cout<<i<<" "<<p[i]<<"\n";

			// save Anderson parameters for this j of the cluster
			saveAndersonParameters(p,j);
		}

		inversionSymmetry(andersonVp_);
		inversionSymmetry(andersonEp_);

		getGammaKOmegaReal();

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

		std::cout<<"andersonVp\n";
		std::cout<<andersonVp_;

		std::cout<<"andersonEp\n";
		std::cout<<andersonEp_;

		// Calculate "hoppings"	
		garbage_ = makeHubbardParams(ekbar,io_);
	}

	void solve(MatrixType& gfCluster)
	{
		if (!garbage_)
			throw PsimagLite::RuntimeError("EffectiveHamiltonian::solve must call build first\n");

		DcaToDmrgType& myInput = *garbage_;
		VaryingGeometryType geometry2(myInput,false,params_.smallKs);

		PsimagLite::String dcaSolver = params_.dcaSolver;
		DcaSolverBaseType* solver = 0;
		if (dcaSolver == "Dmrg") {
			solver = new DcaSolverDmrgType(myInput,geometry2,io_);
		} else if (dcaSolver == "Lanczos") {
			solver = new DcaSolverLanczosType(myInput,geometry2,io_);
		} else {
			PsimagLite::String str("makeHubbardParams(): Unknown solver ");
			str += dcaSolver;
			throw PsimagLite::RuntimeError(str);
		}

		RealType omegaEnd = params_.omegas * params_.omegaStep + params_.omegaBegin;

		PlotParamsType plotParams(params_.omegaBegin,
		                                            omegaEnd,
		                                            params_.omegaStep,
		                                            params_.delta);

		SizeType Nc = params_.largeKs;
		VectorSizeType sites(2);
		for (SizeType i = 0; i < Nc; ++i) {
			for (SizeType j = i; j < Nc; ++j) {
				sites[0] = i;
				sites[1] = j; // dca indexing
				solver->solve(gfCluster,sites,plotParams);
			}
		}

		for (SizeType i = 0; i < Nc; ++i)
			for (SizeType j = 0; j < i; ++j)
				for (SizeType x = 0; x < gfCluster.n_row(); ++x)
					gfCluster(x,i+j*Nc) = gfCluster(x,j+i*Nc);

		delete solver;
	}

	const ComplexType& andersonG0(SizeType omegaIndex, SizeType K) const
	{
		return gammaRealFreq_(omegaIndex,K);
	}

private:

	void saveAndersonParameters(const VectorRealType&src,SizeType k)
	{
		SizeType n=static_cast<SizeType>(src.size()/2);

		for (SizeType i=0;i<n;++i) {
			andersonVp_(i,k) = src[i];
			andersonEp_(i,k) = src[i+n];
		}
	}

	void getGammaKOmegaReal()
	{
		SizeType nBath=andersonVp_.n_row();
		VectorRealType p(2*nBath);

		for (SizeType k=0;k<params_.largeKs;++k) {
			for (SizeType i=0;i<nBath;++i) { //i over bath
				p[i]=std::real(andersonVp_(i,k));
				p[i+nBath]=std::real(andersonEp_(i,k));
			}

			getGammaKOmegaReal(k,p);
		}
	}

	void getGammaKOmegaReal(SizeType k,const VectorRealType& p)
	{
		for (SizeType i=0;i<params_.omegas;++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			ComplexType z(realOmega,params_.delta);
			gammaRealFreq_(i,k)=andersonG0Real(p,z,params_.mu);
		}
	}

	ComplexType andersonG0Real(const VectorRealType&p, ComplexType z,RealType mu) const
	{
		SizeType n=static_cast<SizeType>(p.size()/2);
		ComplexType ctmp=0.0;

		for (SizeType i=0;i<n;++i) ctmp += p[i]*p[i]/(-p[i+n]+z);

		return ctmp;
	}

	DcaToDmrgType* makeHubbardParams(const VectorRealType& ekbar,typename InputNgType::Readable& io)
	{
		SizeType Nc=params_.largeKs;
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;
		MatrixType tCluster(Nc,Nc);
		MatrixType tBathCluster(nBath,Nc*Nc);
		VectorType lambda(Nc*Nc*andersonEp_.n_row());

		ft(tCluster,ekbar);
		ft5(tBathCluster,andersonVp_);
		calcBathLambda(lambda,andersonEp_);

		std::cout<<"tCluster\n";
		std::cout<<tCluster;

		std::cout<<"tBathCluster\n";
		std::cout<<tBathCluster;

		return new DcaToDmrgType(params_,tCluster,tBathCluster,lambda,geometry_,io);
	}

	void ft(MatrixType& dest,const VectorRealType& src) const
	{
		for (SizeType i=0;i<dest.n_row();++i) {
			for (SizeType j=0;j<dest.n_col();++j) {
				dest(i,j)=0.0;
				for (SizeType k=0;k<src.size();++k) {
					RealType tmp = krProduct(k,i) - krProduct(k,j);
					dest(i,j) += ComplexType(cos(tmp),sin(tmp))*src[k];
				}

				dest(i,j) = std::real(dest(i,j))/RealType(geometry_.numberOfSites());
				if (std::norm(dest(i,j))<1e-8) dest(i,j)=0.0;
			}
		}
	}

	void ft5(MatrixType& dest,const MatrixType& src) const
	{
		SizeType Nc = src.n_col();

		for (SizeType b=0;b<src.n_row();++b) { // bath
			for (SizeType i=0;i<Nc;++i) { //cluster real space
				for (SizeType j=0;j<Nc;++j) { //cluster real space
					dest(b,i+j*Nc) = 0.0;
					for (SizeType k=0;k<Nc;++k) { //cluster k-space
						RealType tmp = krProduct(k,i,j);
						dest(b,i+j*Nc)+=ComplexType(cos(tmp),sin(tmp))*src(b,k);
					}

					dest(b,i+j*Nc)=std::real(dest(b,i+j*Nc))/static_cast<RealType>(geometry_.numberOfSites());
				}
			}
		}
	}

	void calcBathLambda(VectorType& lambda,const MatrixType& src) const
	{
		SizeType Nc=params_.largeKs;
		assert(lambda.size() == Nc*Nc*src.n_row());

		for (SizeType i=0;i<src.n_row();i++) { // i is in the bath
			for (SizeType j=0;j<Nc;j++) { // j is in the cluster (real space)
				for (SizeType s=0;s<Nc;s++) { // j is in the cluster (real space)
					lambda[j+s*Nc+i*Nc*Nc]=0.0;
					for (SizeType k=0;k<Nc;k++) { // k is in the cluster (k-space)
						RealType tmp = krProduct(k,j,s);
						lambda[j+s*Nc+i*Nc*Nc]+= ComplexType(cos(tmp),sin(tmp))*src(i,k);
					}

					lambda[j+s*Nc+i*Nc*Nc] /= static_cast<RealType>(geometry_.numberOfSites());
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

	void inversionSymmetry(MatrixType& andersonp)
	{
		MatrixType tempMatrix(andersonp.n_row(),andersonp.n_col());

		for (SizeType i=0;i<andersonp.n_row();++i) {
			for (SizeType j=0;j<andersonp.n_col();++j) {
				tempMatrix(i,j)=0.0;
				for (SizeType op=0;op<geometry_.nGroupK();++op) {
					SizeType k = geometry_.ickequ(j,op);
					tempMatrix(i,j) += andersonp(i,k);
				}
				tempMatrix(i,j) /= static_cast<RealType>(geometry_.nGroupK());
			}
		}

		for (SizeType i=0;i<andersonp.n_row();++i)
			for (SizeType j=0;j<andersonp.n_col();++j)
				andersonp(i,j)=tempMatrix(i,j); // andersonp = tempMatrix		
	}

	const ParametersType& params_;
	const GeometryType& geometry_;
	typename InputNgType::Readable& io_;
	MatrixType andersonVp_;
	MatrixType andersonEp_;
	MatrixType gammaRealFreq_;
	DcaToDmrgType* garbage_;
};

template<typename RealType, typename InputNgType>
std::ostream& operator<<(std::ostream& os,EffectiveHamiltonian<RealType,InputNgType>& params)
{
	os<<"EffectiveHamiltonian\n";
	return os;
}

} // namespace DynamicClusterApprox

#endif // EFFECTIVE_HAMILTONIAN_H

