#ifndef DCA_LOOP_H
#define DCA_LOOP_H
#include "Matrix.h"
#include "Vector.h"
#include "EffectiveHamiltonian.h"

namespace OpenDca {

template<typename ParamsType, typename GeometryType>
class Dispersion {

	typedef typename ParamsType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	Dispersion(const ParamsType& params, const GeometryType& geometry)
	: params_(params),geometry_(geometry)
	{}

	RealType operator()(SizeType coarseIndex,SizeType fineIndex) const
	{
		VectorRealType kvector(geometry_.dimension());
		geometry_.index2Kvector(coarseIndex,kvector);

		VectorRealType kvector2(geometry_.dimension());
		geometry_.getMeshVector(kvector2,fineIndex);

		for (SizeType i=0;i<kvector.size();++i)
			kvector[i] += kvector2[i];

		return getDispersion(kvector);
	}

private:

	RealType getDispersion(const VectorRealType& kvector) const
	{
		SizeType ndim=geometry_.dimension();
		RealType kplus = 0;

		for (SizeType k=0;k<ndim;++k)
			kplus += cos(kvector[k]);

		return -0.5*kplus;
	}

	const ParamsType& params_;
	const GeometryType& geometry_;
};

template<typename RealType,typename InputNgType>
class DcaLoop {

	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorType;
	typedef EffectiveHamiltonian<RealType,InputNgType> EffectiveHamiltonianType;
	typedef typename EffectiveHamiltonianType::GeometryType GeometryType;
	typedef typename EffectiveHamiltonianType::ParametersType ParametersType_;
	typedef Dispersion<ParametersType_,GeometryType> DispersionType;
	typedef typename EffectiveHamiltonianType::DcaToDmrgType DcaToDmrgType;

public:

	typedef ParametersType_ ParametersType;

	DcaLoop(const ParametersType& params,
	               typename InputNgType::Readable& io)
	: params_(params),
	  io_(io),
	  geometry_(io,false,params_.smallKs),
	  dispersion_(params,geometry_),
	  sigma_(params.omegas,params.largeKs),
	  gckfsc_(params.omegas,params.largeKs),
	  fTCoefsR2K_(params.largeKs,params.largeKs)
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
	}

	void main()
	{
		SizeType largeKs = params_.largeKs;
		SizeType totalOmegas = params_.omegas;
		VectorRealType barEpsilon(largeKs);
		coarseDispersion(barEpsilon);
		std::cout<<"barEpsilon\n";
		std::cout<<barEpsilon;
		MatrixType gfcluster(totalOmegas,largeKs*largeKs);
		MatrixType gammaOmega(params_.numberOfMatsubaras,gckfsc_.n_col());
		RealType sigmaNorm = 0;

		for (SizeType i = 0; i < params_.iterations; ++i) {
			makeGf();
			diagUpdate(gfcluster,gammaOmega,barEpsilon);
			std::cout<<"gfcluster\n";
			std::cout<<gfcluster;
			RealType dcaError = sigmaNorm;
			sigmaNorm = makeSigma(gfcluster,gammaOmega,barEpsilon);
			dcaError -= sigmaNorm;
			std::cout<<"sigma\n";
			std::cout<<sigma_;
			std::cout<<"Dca iteration= "<<i<<" error in sigma= "<<fabs(dcaError)<<"\n";
		}
	}

	void makeGf()
	{
		VectorType gckf(params_.omegas);
		SizeType largeKs = params_.largeKs;
		SizeType totalOmegas = params_.omegas;

		for (SizeType k = 0; k < largeKs; ++k) {
			makeGf(gckf,k);
			for (SizeType omegaIndex = 0; omegaIndex < totalOmegas; ++omegaIndex)
				gckfsc_(omegaIndex,k) = gckf[omegaIndex]/(1.0+sigma_(omegaIndex,k)*gckf[omegaIndex]);
		}
	}

	void diagUpdate(MatrixType& gfCluster,
	                           MatrixType& gammaOmega,
	                          const VectorRealType& ekbar)
	{
		VectorRealType integral(ekbar.size());
		MatrixType deltaOmega(gckfsc_.n_row(),gckfsc_.n_col());
		MatrixType gfLesser;
		EffectiveHamiltonianType effectiveHamiltonian(params_,geometry_,io_);

		std::cerr<<"gamma omega...\n";
		getDeltaOmega(deltaOmega,integral,ekbar);
		getGammaOmega(gammaOmega,deltaOmega);
		std::cerr<<"make h params...\n";
		effectiveHamiltonian.build(gammaOmega,ekbar,integral);
		std::cout<<effectiveHamiltonian;

		std::cerr<<"lanczos started...\n";
		effectiveHamiltonian.solve(gfCluster);
		std::cerr<<"lanczos done\n";
	}

private:

	void makeGf(VectorType& gckf,SizeType K)
	{
		SizeType totalOmegas = params_.omegas;
		SizeType meshPoints = geometry_.sizeOfMesh();
		for (SizeType omegaIndex = 0; omegaIndex < totalOmegas; ++omegaIndex) {
			gckf[omegaIndex] = 0.0;
			for (SizeType ktilde = 0; ktilde < meshPoints; ++ktilde) {
				RealType omega = params_.omegaBegin + params_.omegaStep * omegaIndex;
				ComplexType tmp = -dispersion_(K,ktilde) - params_.mu
				                                         + ComplexType(omega,params_.delta)
				                                         - sigma_(omegaIndex,K);
				gckf[omegaIndex] += 1.0/tmp;
			}
		}

		std::cout<<"gckf for k = "<<K<<"\n";
		std::cout<<gckf;
	}

	void coarseDispersion(VectorRealType& barEpsilon)
	{
		SizeType meshPoints = geometry_.sizeOfMesh();
		SizeType largeKs = params_.largeKs;
		assert(barEpsilon.size() == largeKs);

		for (SizeType bigK = 0; bigK < largeKs; ++bigK) {
			barEpsilon[bigK] = 0.0;
			for (SizeType ktilde = 0; ktilde < meshPoints; ++ktilde) {
				barEpsilon[bigK] += dispersion_(bigK,ktilde);
			}
		}
	}

	//! Delta_k(omega) = omega + ekbar(k) -sigma(omega,k) - 1.0/gckfsc_(omega,k)
	void getDeltaOmega(MatrixType& deltaOmega,
	                                  VectorRealType& integral,
	                                  const VectorRealType& ekbar)
	{
		for (SizeType j=0;j<gckfsc_.n_col();j++) integral[j]=0.0;

		for (SizeType i=0;i<gckfsc_.n_row();++i) { // loop over freq.
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			for (SizeType j=0;j<gckfsc_.n_col();++j) { // loop over k
				deltaOmega(i,j)= realOmega - ekbar[j] -sigma_(i,j) - 1.0/gckfsc_(i,j);
				if (std::imag(deltaOmega(i,j))>0)
					deltaOmega(i,j)=std::real(deltaOmega(i,j));
				integral[j] += std::imag(deltaOmega(i,j));
			}
		}

		std::cout<<"#DELTAOMEGA\n";
		for (SizeType i=0;i<gckfsc_.n_row();++i) { // loop over freq.
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j=0;j<gckfsc_.n_col();++j)
				std::cout<<std::imag(deltaOmega(i,j))<<" ";
			std::cout<<"\n";
		}

		for (SizeType j=0;j<gckfsc_.n_col();++j) {
			integral[j]=integral[j]*(params_.omegaStep);
			std::cerr<<"integral["<<j<<"]="<<integral[j]<<"\n";
		}
	}

	void getGammaOmega(MatrixType& gammaOmega,
	                   const MatrixType& deltaOmega)
	{
		RealType factor = -params_.omegaStep/M_PI;

		for (SizeType i=0;i<params_.numberOfMatsubaras;++i) {
			RealType wn = matSubara(i);

			for (SizeType k=0;k<deltaOmega.n_col();++k) {
				ComplexType sum = 0.0;
				for (SizeType j=0;j<deltaOmega.n_row();++j) {
					RealType realOmega = params_.omegaBegin + params_.omegaStep * j;
					sum += factor * std::imag(deltaOmega(j,k))/ComplexType(-realOmega,wn);
				}

				gammaOmega(i,k) = sum/M_PI;
			}
		}

		std::cout<<"#GAMMA\n";
		for (SizeType i=0;i<gammaOmega.n_row();++i) {
			RealType wn = matSubara(i);
			std::cout<<wn<<" ";
			for (SizeType j=0;j<params_.largeKs;++j)
				std::cout<<gammaOmega(i,j)<<" ";
			std::cout<<"\n";
		}
	}

	RealType matSubara(int ind) const
	{
		int halfNs = static_cast<int>(params_.numberOfMatsubaras*0.5);
		return (2.0*(ind-halfNs)+1.0)*M_PI/params_.beta;
	}

	RealType makeSigma(const MatrixType& gfCluster,
	                   const MatrixType& gammakomega,
	                   const VectorRealType& epsbar)
	{
		MatrixType& sigma = sigma_;
		MatrixType data(gfCluster.n_row(),params_.largeKs);

		// transform from real to k space
		ft4(data,gfCluster,fTCoefsR2K_,params_.largeKs);

		std::cout<<"#DATAIMAG\n";
		for (SizeType i = 0;i < data.n_row(); ++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j = 0;j < data.n_col(); ++j) {
				if (std::imag(data(i,j))>0) data(i,j)=std::real(data(i,j));
				std::cout<<imag(data(i,j))<<" ";
			}

			std::cout<<"\n";
		}

		std::cout<<"#DATAREAL\n";
		for (SizeType i = 0;i < data.n_row(); ++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j = 0;j < data.n_col(); ++j)
				std::cout<<real(data(i,j))<<" ";
			std::cout<<"\n";
		}

		std::cout<<"#ONEOVERDATAREAL\n";
		for (SizeType i = 0;i < data.n_row(); ++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j = 0;j < data.n_col(); ++j)
				std::cout<<std::real(1.0/data(i,j))<<" ";
			std::cout<<"\n";
		}

		std::cout<<"#ONEOVERDATAIMAG\n";
		for (SizeType i = 0; i < data.n_row(); ++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j = 0;j < data.n_col(); ++j)
				std::cout<<std::imag(1.0/data(i,j))<<" ";
			std::cout<<"\n";
		}

		for (SizeType i = 0;i < sigma.n_row(); ++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			for (SizeType j=0;j < sigma.n_col(); ++j) {
				sigma(i,j) = realOmega - epsbar[j] - gammakomega(i,j)
				                   - static_cast<RealType>(params_.largeKs)/data(i,j);
				if (std::imag(sigma(i,j))>0) sigma(i,j) = std::real(sigma(i,j));
			}
		}

		std::cout<<"#SIGMAIMAG\n";
		for (SizeType i=0;i<sigma.n_row();++i) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j=0;j<sigma.n_col();j++) std::cout<<imag(sigma(i,j))<<" ";
			std::cout<<"\n";
		}
		std::cout<<"#SIGMAREAL\n";
		for (SizeType i=0;i<sigma.n_row();i++) {
			RealType realOmega = params_.omegaBegin + params_.omegaStep * i;
			std::cout<<realOmega<<" ";
			for (SizeType j=0;j<sigma.n_col();j++)
				std::cout<<real(sigma(i,j))<<" ";
			std::cout<<"\n";
		}

		std::cout<<"#\n";

		// test for convergence by inspecting Im[sigma(0,ic) ]           
		return calcRealSigma(sigma);
	}

	RealType calcRealSigma(const MatrixType& sigma) const
	{
		RealType tmp =0.0;
		int indexOfZero=static_cast<int>((sigma.n_row()-1)*0.5);
		for (SizeType i = 0; i < params_.largeKs; ++i)
			tmp += fabs(std::real(sigma(indexOfZero,i)));

		return tmp/params_.largeKs;
	}

	void ft4(MatrixType& gfdest,
	         const MatrixType& gfsource,
	         const MatrixType& ftCoeffs,
	         SizeType Nc) const
	{
		ComplexType csum = 0.0;

		for (SizeType ik=0;ik<Nc;ik++) {
			for (SizeType n=0;n<gfsource.n_row();n++) {
				csum=0.0;
				for (SizeType ic=0;ic<Nc;ic++)
					for (SizeType jc=0;jc<Nc;jc++)
						csum+= gfsource(n,ic+jc*Nc)*conj(ftCoeffs(ik,ic))*ftCoeffs(ik,jc);
				gfdest(n,ik)=csum; //! no factor of 1/Nc
			}
		}
	}

	const ParametersType& params_;
	typename InputNgType::Readable& io_;
	GeometryType geometry_;
	DispersionType dispersion_;
	MatrixType sigma_;
	MatrixType gckfsc_;
	MatrixType fTCoefsR2K_;
}; // class DcaLoop

} // namespace OpenDca

#endif // DCA_LOOP_H

