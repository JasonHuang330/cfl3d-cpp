// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "init_trim.h"
#include "ccomplex.h"
#include "runtime/common_blocks.h"
#include "runtime/fortran_array.h"
#include <complex>
#include <cstdint>

namespace init_trim_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q,
          FortranArray4DRef<double> qj0,
          FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0,
          FortranArray4DRef<double> tj0,
          FortranArray4DRef<double> tk0,
          FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol,
          FortranArray3DRef<double> volj0,
          FortranArray3DRef<double> volk0,
          FortranArray3DRef<double> voli0,
          int& nummem,
          FortranArray3DRef<double> x,
          FortranArray3DRef<double> y,
          FortranArray3DRef<double> z,
          FortranArray1DRef<int> nou,
          FortranArray2DRef<char[120]> bou,
          int& nbuf,
          int& ibufdim,
          int& iflagprnt)
{
    (void)nbl; (void)jdim; (void)kdim; (void)idim;
    (void)q; (void)qj0; (void)qk0; (void)qi0;
    (void)tj0; (void)tk0; (void)ti0;
    (void)vol; (void)volj0; (void)volk0; (void)voli0;
    (void)nummem;
    (void)x; (void)y; (void)z;
    (void)nou; (void)bou;
    (void)nbuf; (void)ibufdim;
    (void)iflagprnt;
}

void init_trim()
{
    // Aliases for /trim/ COMMON block fields
    float& dmtrmn   = cmn_trim.dmtrmn;
    float& dmtrmnm  = cmn_trim.dmtrmnm;
    float& dlcln    = cmn_trim.dlcln;
    float& dlclnm   = cmn_trim.dlclnm;
    float& trtol    = cmn_trim.trtol;
    float& dtrmsmx  = cmn_trim.dtrmsmx;
    float& dtrmsmn  = cmn_trim.dtrmsmn;
    float& dalfmx   = cmn_trim.dalfmx;
    float& ddtmx    = cmn_trim.ddtmx;
    int32_t& itrmt  = cmn_trim.itrmt;
    float& epstr    = cmn_trim.epstr;
    // fp and tp as flat float[16] column-major arrays
    float* fp_flat  = cmn_trim.fp;
    float* tp_flat  = cmn_trim.tp;

    // Aliases for /trim1/ COMMON block fields
    float& a11      = cmn_trim1.a11;
    float& a12      = cmn_trim1.a12;
    float& a22      = cmn_trim1.a22;
    float& r11      = cmn_trim1.r11;
    float& r22      = cmn_trim1.r22;
    float& r33      = cmn_trim1.r33;
    float& r44      = cmn_trim1.r44;

    //
    //     The following relaxation coefficients and limits are hardwired:
    //
    dtrmsmx =  0.50f;
    dtrmsmn = -0.50f;
    dmtrmn  =  0.f;
    dmtrmnm =  0.f;
    dlcln   =  0.f;
    dlclnm  =  0.f;
    ddtmx   =  0.0005f;
    dalfmx  =  0.0005f;
    trtol   =  0.0000002f;
    itrmt   =  0;
    epstr   =  1.e-7f;
    a11     =  0.f;
    a12     =  0.f;
    a22     =  0.f;
    r11     =  0.f;
    r22     =  0.f;
    r33     =  0.f;
    r44     =  0.f;

    // Fortran: dtr = tp(1,1)*tp(2,2)-tp(1,2)*tp(2,1)
    // tp(i,j) column-major: tp[(j-1)*4 + (i-1)]
    // tp(1,1) -> tp_flat[0], tp(2,2) -> tp_flat[5]
    // tp(1,2) -> tp_flat[4], tp(2,1) -> tp_flat[1]
    double dtr = (double)tp_flat[0] * (double)tp_flat[5]
               - (double)tp_flat[4] * (double)tp_flat[1];

    // Fortran: dtr = dtr + ccsignrc(1e-30,dtr)
    // ccsignrc takes double& args
    double ccsignrc_a = 1.e-30;
    double ccsignrc_b = dtr;
    dtr = dtr + ccomplex_ns::ccsignrc(ccsignrc_a, ccsignrc_b);

    // Fortran: fp(1,1) = tp(2,2)/dtr
    // fp(i,j) column-major: fp[(j-1)*4 + (i-1)]
    // fp(1,1) -> fp_flat[0]
    fp_flat[0] = (float)((double)tp_flat[5] / dtr);

    // Fortran: fp(1,2) = -tp(1,2)/dtr
    // fp(1,2) -> fp_flat[4]
    fp_flat[4] = (float)(-(double)tp_flat[4] / dtr);

    // Fortran: fp(2,1) = -tp(2,1)/dtr
    // fp(2,1) -> fp_flat[1]
    fp_flat[1] = (float)(-(double)tp_flat[1] / dtr);

    // Fortran: fp(2,2) = tp(1,1)/dtr
    // fp(2,2) -> fp_flat[5]
    fp_flat[5] = (float)((double)tp_flat[0] / dtr);

    return;
}

} // namespace init_trim_ns
