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
#include "FreqEnum.h"

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
	  p_(2*params_.nofPointsInBathPerClusterPoint,params_.largeKs*params_.orbitals),
	  gammaRealFreq_(params_.omegas,params_.largeKs*params_.orbitals),
	  garbage_(0)
	{}

	~EffectiveHamiltonian()
	{
		delete garbage_;
	}

	// Get LanczosParams from Gf
	void build(MatrixType& gf,
	           const VectorRealType& ekbar,
	           const VectorRealType& integral,
	           FreqEnum freqEnum)
	{
		SizeType nBath=params_.nofPointsInBathPerClusterPoint;

		VectorRealType p(2*nBath);
		VectorType gfTmp(gf.n_row());
		AndersonFitType andersonFit(params_);

		assert(p_.n_col()*p_.n_col() == gf.n_col());
		for (SizeType j=0;j<p_.n_col();++j) {
			SizeType jj = j+j*params_.orbitals;
			// first find Anderson parameters by fitting procedure (min. of distance)
			std::cout<<"#GFTMP"<<j<<"\n";
			for (SizeType i=0;i<gf.n_row();++i) {
				gfTmp[i]=gf(i,jj);
				std::cout<<i<<" "<<std::real(gfTmp[i])<<" "<<std::imag(gfTmp[i])<<"\n";
			}

			andersonFit.fit(p,gfTmp,integral[jj]);
			std::cout<<"#PANDERSON"<<j<<"\n";
			for (SizeType i=0;i<p.size();++i) std::cout<<i<<" "<<p[i]<<"\n";

			// save Anderson parameters for this j of the cluster
			saveAndersonParameters(p,j);
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

	const MatrixType& andersonParameters() const
	{
		return p_;
	}

private:

	void saveAndersonParameters(const VectorRealType&src,SizeType k)
	{
		SizeType n=src.size();
		assert(p_.n_row() == n);
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

	const ParametersType& params_;
	const GeometryType& geometry_;
	typename InputNgType::Readable& io_;
	MatrixType p_;
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

