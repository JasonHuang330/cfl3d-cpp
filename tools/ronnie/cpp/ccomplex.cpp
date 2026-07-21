// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ccomplex.cpp — translation of ccomplex.F (CFL3D ronnie).
//
// ccomplex.F provides complex-arithmetic helper FUNCTIONs (ccabs, ccsign,
// ccmax, ccmin, cctanh, ccatan2, ...) that only exist in the CMPLX build.
// In the SEQUENTIAL REAL build every one reduces to a trivial real intrinsic
// under its #else branch, e.g.:
//     function ccabs(a)
//     ccabs = abs(a)
//
// A grep of the entire ronnie source tree (all .F except ccomplex.F) shows
// that the ONLY ccomplex symbol referenced anywhere in the real-build call
// graph is `ccabs` (30 call sites).  `ccabs` is already provided as an inline
// in ron_common.h:
//     inline double ccabs(double x){ return std::fabs(x); }
// which is exactly the real-build definition (abs(a)).
//
// None of the other functions (ccsign, ccsignrc, ccmax, ccmax4, ccmaxcr,
// ccmaxrc, ccmin, ccmin4, ccmin8, ccmincr, ccminrc, ccvmgt, cctanh, ccatan2,
// ccatan, ccacos, ccasin, cctan, cclog10, ccdim, ccerf) are called outside
// ccomplex.F itself, so none are emitted.
//
// This translation unit is therefore intentionally empty (no symbols kept):
// ccabs stays an inline in ron_common.h; everything else is CMPLX-only dead
// code and is omitted.
#include "ron_common.h"
