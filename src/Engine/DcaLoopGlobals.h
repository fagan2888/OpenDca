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

/** \ingroup DMRG */
/*@{*/

/*! \file DcaLoopGlobals.h
 *
 *
 *
 */
#ifndef DCA_LOOP_GLOBALS_H
#define DCA_LOOP_GLOBALS_H

#include "Vector.h"

namespace OpenDca {

enum SpinEnum {SPIN_UP, SPIN_DOWN};

struct DcaLoopGlobals {

	static const PsimagLite::String license;

}; // DcaLoopGlobals

const PsimagLite::String DcaLoopGlobals::license=
"Copyright (c) 2014, UT-Battelle, LLC\n"
"\n"
"OpenDca, Version 1.0\n"
"\n"
"OpenDca is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
"GNU General Public License for more details.\n"
"\n";

} // namespace OpenDca
/*@}*/
#endif

