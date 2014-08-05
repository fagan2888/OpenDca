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
#include "GitRevision.h"

int main(int argc,char *argv[])
{
	PsimagLite::GitRevision gitrev("./","openDca");
	std::cout<<gitrev;
	PsimagLite::GitRevision gitrev2("../../PsimagLite/","psimagLite");
	std::cout<<gitrev2;

}

