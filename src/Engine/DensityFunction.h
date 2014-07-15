#ifndef OPENDCA_DENSITY_FUNC_H
#define OPENDCA_DENSITY_FUNC_H

class DensityFunction {

public:

	typedef double RealType;

	RealType operator()(RealType mu) const
	{
		return 0.0;
	}

}; // DensityFunction

#endif // OPENDCA_DENSITY_FUNC_H

