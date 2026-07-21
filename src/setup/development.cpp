// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "development.h"

namespace development_ns {

// subroutine dummy
//
// $Id$
//
// Purpose: this is a reserved location for any development code.
// Experimental code may be placed after this routine. For example,
// if a new version of a subroutine is to be tested, it may be placed
// below. When compiled, the version below will be used and the
// original version ignored (at least this is true for most compilers).
// The exception to this is the main routine - any changes to main
// must be made in main.F itself.
//
// Note: this module (development.F) should always contain at least
// the following lines (without the comment characters):
//
//     subroutine dummy
//#  ifdef CMPLX
//     implicit complex(a-h,o-z)
//#  endif
//     return
//     end

void dummy()
{
    return;
}

} // namespace development_ns
