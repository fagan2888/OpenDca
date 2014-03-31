#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "Version.h" // do not commit this file, created dynamically

class Provenance {

public:

}; // Provenance

std::ostream& operator<<(std::ostream& os,const Provenance &prov)
{
	os<<"OpenDca: revision: "<<dynamicClusterApproxRevision<<"\n";
	os<<"OpenDca: diff: "<<dynamicClusterApproxDiff<<"\n";
	os<<"PsimagLite: revision: "<<psimagLiteRevision<<"\n";
	os<<"PsimagLite: diff: "<<psimagLiteDiff<<"\n";
	return os;
}

#endif // PROVENANCE_H

