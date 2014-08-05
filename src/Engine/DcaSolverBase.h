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
#ifndef DCA_SOLVER_BASE_H
#define DCA_SOLVER_BASE_H
#include "Vector.h"
#include "TridiagonalMatrix.h"
#include "ContinuedFractionCollection.h"
#include "ContinuedFraction.h"

namespace OpenDca {

template<typename DcaToDmrgType_,typename VaryingGeometryType_>
class DcaSolverBase {

public:

	typedef DcaToDmrgType_ DcaToDmrgType;
	typedef VaryingGeometryType_ VaryingGeometryType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename DcaToDmrgType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::Matrix<ComplexType> MatrixType;
	typedef PsimagLite::TridiagonalMatrix<RealType> TridiagonalMatrixType;
	typedef PsimagLite::ContinuedFraction<TridiagonalMatrixType> ContinuedFractionType;
	typedef PsimagLite::ContinuedFractionCollection<ContinuedFractionType>
	                    ContinuedFractionCollectionType;
	typedef typename ContinuedFractionType::PlotParamsType PlotParamsType;

	virtual ~DcaSolverBase() {}

	virtual void solve(MatrixType& gf,
	                   const VectorSizeType& sites,
	                   const PlotParamsType& plotParams) = 0;

	virtual RealType findLowestEnergy() const = 0;

	virtual RealType density(SizeType total1, SizeType total2) const = 0;

};
}

#endif

