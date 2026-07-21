// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pre_cblki.h"
#include "termn8.h"
#include <algorithm>
#include <cstdio>

namespace pre_cblki_ns {

void pre_cblki(int& idimr, int& jdimr, int& kdimr, FortranArray2DRef<int> limblk, FortranArray2DRef<int> isva, int& it, int& ir, FortranArray3DRef<double> xr, FortranArray3DRef<double> yr, FortranArray3DRef<double> zr, FortranArray2DRef<double> xyzn, int& maxdims, int& ij, int& myid, int& ibufdim, int& nbuf, FortranArray2DRef<char[120]> bou, FortranArray1DRef<int> nou)
{
    int jst, jet, jinct;
    int kst, ket, kinct;
    double eps;
    int isr, ier, jsr, jer, ksr, ker;
    int jincr, iincr, kincr;
    int kloc1r, kloc2r, kr;
    int jloc1r, jloc2r, jr;
    int iloc1r, iloc2r, irn;
    int icount, jcount, kcount;
    int ilocr, jlocr, klocr;
    int termn8_err;

    jst = limblk(it,2);
    jet = limblk(it,5);
    if (jst == jet) {
        jinct = 1;
    } else {
        jinct = (jet-jst)/std::abs(jet-jst);
    }

    kst = limblk(it,3);
    ket = limblk(it,6);
    if (kst == ket) {
        kinct = 1;
    } else {
        kinct = (ket-kst)/std::abs(ket-kst);
    }

    eps = 0.;

    isr = limblk(ir,1);
    ier = limblk(ir,4);
    jsr = limblk(ir,2);
    jer = limblk(ir,5);
    ksr = limblk(ir,3);
    ker = limblk(ir,6);

    //
    // determine the side of the q array to transfer from
    //
    // k = constant side
    //
    if (isva(ir,1)+isva(ir,2) == 3) {
        if (ksr == 1) {
            kloc1r = 1;
            kloc2r = 2;
        } else {
            kloc1r = kdimr-1;
            kloc2r = kdimr-2;
        }

        if (kdimr == 2) {
            kloc1r = 1;
            kloc2r = 1;
        }

        if (jer == jsr) {
            jincr = 1;
        } else {
            jincr = (jer-jsr)/std::abs(jer-jsr);
        }

        if (ier == isr) {
            iincr = 1;
        } else {
            iincr = (ier-isr)/std::abs(ier-isr);
        }

        kr = kloc1r;
        if (kr > 1) kr = kdimr;
        if ((isva(ir,1) == isva(it,1)) ||
            (isva(ir,2) == isva(it,2))) {
            //
            // i varies with k     and     j varies with j
            //
            ij = 0;
            icount = -1;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                icount = icount + 1;
                jcount = -1;
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    jcount = jcount + 1;
                    ilocr  = isr + iincr*icount;
                    jlocr  = jsr + jincr*jcount;
                    ij = ij + 1;
                    if (ij > maxdims) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " stopping in pre_cblki ij.gt.maxdims");
                        termn8_err = -1;
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    xyzn(ij,1) = .25*( xr(jlocr  ,kr,ilocr) +
                                       xr(jlocr  ,kr,ilocr+1) +
                                       xr(jlocr+1,kr,ilocr) +
                                       xr(jlocr+1,kr,ilocr+1) );
                    xyzn(ij,2) = .25*( yr(jlocr  ,kr,ilocr) +
                                       yr(jlocr  ,kr,ilocr+1) +
                                       yr(jlocr+1,kr,ilocr) +
                                       yr(jlocr+1,kr,ilocr+1) );
                    xyzn(ij,3) = .25*( zr(jlocr  ,kr,ilocr) +
                                       zr(jlocr  ,kr,ilocr+1) +
                                       zr(jlocr+1,kr,ilocr) +
                                       zr(jlocr+1,kr,ilocr+1) );
                }
            }
        } else {
            //
            // j varies with k     and     i varies with j
            //
            ij = 0;
            jcount = -1;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                jcount = jcount + 1;
                icount = -1;
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    icount = icount + 1;
                    ilocr  = isr + iincr*icount;
                    jlocr  = jsr + jincr*jcount;
                    ij = ij + 1;
                    if (ij > maxdims) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " stopping in pre_cblki ij.gt.maxdims");
                        termn8_err = -1;
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    xyzn(ij,1) = .25*( xr(jlocr  ,kr,ilocr) +
                                       xr(jlocr  ,kr,ilocr+1) +
                                       xr(jlocr+1,kr,ilocr) +
                                       xr(jlocr+1,kr,ilocr+1) );
                    xyzn(ij,2) = .25*( yr(jlocr  ,kr,ilocr) +
                                       yr(jlocr  ,kr,ilocr+1) +
                                       yr(jlocr+1,kr,ilocr) +
                                       yr(jlocr+1,kr,ilocr+1) );
                    xyzn(ij,3) = .25*( zr(jlocr  ,kr,ilocr) +
                                       zr(jlocr  ,kr,ilocr+1) +
                                       zr(jlocr+1,kr,ilocr) +
                                       zr(jlocr+1,kr,ilocr+1) );
                }
            }
        }
    //
    // j = constant side
    //
    } else if (isva(ir,1)+isva(ir,2) == 4) {
        if (jsr == 1) {
            jloc1r = 1;
            jloc2r = 2;
        } else {
            jloc1r = jdimr-1;
            jloc2r = jdimr-2;
        }

        if (jdimr == 2) {
            jloc1r = 1;
            jloc2r = 1;
        }

        if (ier == isr) {
            iincr = 1;
        } else {
            iincr = (ier-isr)/std::abs(ier-isr);
        }

        if (ker == ksr) {
            kincr = 1;
        } else {
            kincr = (ker-ksr)/std::abs(ker-ksr);
        }

        jr = jloc1r;
        if (jr > 1) jr = jdimr;

        if ((isva(ir,1) == isva(it,1)) ||
            (isva(ir,2) == isva(it,2))) {
            //
            // i varies with j    and    k varies with k
            //
            ij = 0;
            kcount = -1;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                kcount = kcount + 1;
                icount = -1;
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    icount = icount + 1;
                    ilocr  = isr + iincr*icount;
                    klocr  = ksr + kincr*kcount;
                    ij = ij + 1;
                    if (ij > maxdims) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " stopping in pre_cblki ij.gt.maxdims");
                        termn8_err = -1;
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    xyzn(ij,1) = .25*( xr(jr,klocr  ,ilocr) +
                                       xr(jr,klocr  ,ilocr+1) +
                                       xr(jr,klocr+1,ilocr) +
                                       xr(jr,klocr+1,ilocr+1) );
                    xyzn(ij,2) = .25*( yr(jr,klocr  ,ilocr) +
                                       yr(jr,klocr  ,ilocr+1) +
                                       yr(jr,klocr+1,ilocr) +
                                       yr(jr,klocr+1,ilocr+1) );
                    xyzn(ij,3) = .25*( zr(jr,klocr  ,ilocr) +
                                       zr(jr,klocr  ,ilocr+1) +
                                       zr(jr,klocr+1,ilocr) +
                                       zr(jr,klocr+1,ilocr+1) );
                }
            }
        } else {
            //
            // i varies with k    and    k varies with j
            //
            ij = 0;
            icount = -1;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                icount = icount + 1;
                kcount = -1;
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    kcount = kcount + 1;
                    ilocr  = isr + iincr*icount;
                    klocr  = ksr + kincr*kcount;
                    ij = ij + 1;
                    if (ij > maxdims) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " stopping in pre_cblki ij.gt.maxdims");
                        termn8_err = -1;
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    xyzn(ij,1) = .25*( xr(jr,klocr  ,ilocr) +
                                       xr(jr,klocr  ,ilocr+1) +
                                       xr(jr,klocr+1,ilocr) +
                                       xr(jr,klocr+1,ilocr+1) );
                    xyzn(ij,2) = .25*( yr(jr,klocr  ,ilocr) +
                                       yr(jr,klocr  ,ilocr+1) +
                                       yr(jr,klocr+1,ilocr) +
                                       yr(jr,klocr+1,ilocr+1) );
                    xyzn(ij,3) = .25*( zr(jr,klocr  ,ilocr) +
                                       zr(jr,klocr  ,ilocr+1) +
                                       zr(jr,klocr+1,ilocr) +
                                       zr(jr,klocr+1,ilocr+1) );
                }
            }
        }
    //
    // i = constant side
    //
    } else if (isva(ir,1)+isva(ir,2) == 5) {
        if (isr == 1) {
            iloc1r = 1;
            iloc2r = 2;
        } else {
            iloc1r = idimr-1;
            iloc2r = idimr-2;
        }

        if (idimr == 2) {
            iloc1r = 1;
            iloc2r = 1;
        }

        if (jer == jsr) {
            jincr = 1;
        } else {
            jincr = (jer-jsr)/std::abs(jer-jsr);
        }

        if (ker == ksr) {
            kincr = 1;
        } else {
            kincr = (ker-ksr)/std::abs(ker-ksr);
        }

        irn = iloc1r;
        if (irn > 1) irn = idimr;
        if ((isva(ir,1) == isva(it,1)) ||
            (isva(ir,2) == isva(it,2))) {
            //
            // k varies with k    and    j varies with j
            //
            ij = 0;
            kcount = -1;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                kcount = kcount + 1;
                jcount = -1;
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    jcount = jcount + 1;
                    jlocr  = jsr + jincr*jcount;
                    klocr  = ksr + kincr*kcount;
                    ij = ij + 1;
                    if (ij > maxdims) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " stopping in pre_cblki ij.gt.maxdims");
                        termn8_err = -1;
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    xyzn(ij,1) = .25*( xr(jlocr  ,klocr,irn) +
                                       xr(jlocr  ,klocr+1,irn) +
                                       xr(jlocr+1,klocr,irn) +
                                       xr(jlocr+1,klocr+1,irn) );
                    xyzn(ij,2) = .25*( yr(jlocr  ,klocr,irn) +
                                       yr(jlocr  ,klocr+1,irn) +
                                       yr(jlocr+1,klocr,irn) +
                                       yr(jlocr+1,klocr+1,irn) );
                    xyzn(ij,3) = .25*( zr(jlocr  ,klocr,irn) +
                                       zr(jlocr  ,klocr+1,irn) +
                                       zr(jlocr+1,klocr,irn) +
                                       zr(jlocr+1,klocr+1,irn) );
                }
            }
        } else {
            //
            // j varies with k    and    k varies with j
            //
            ij = 0;
            jcount = -1;
            for (int k = kst; (kinct > 0) ? (k <= ket) : (k >= ket); k += kinct) {
                jcount = jcount + 1;
                kcount = -1;
                for (int j = jst; (jinct > 0) ? (j <= jet) : (j >= jet); j += jinct) {
                    kcount = kcount + 1;
                    jlocr  = jsr + jincr*jcount;
                    klocr  = ksr + kincr*kcount;
                    ij = ij + 1;
                    if (ij > maxdims) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " stopping in pre_cblki ij.gt.maxdims");
                        termn8_err = -1;
                        termn8_ns::termn8(myid, termn8_err, ibufdim, nbuf, bou, nou);
                    }
                    xyzn(ij,1) = .25*( xr(jlocr  ,klocr,irn) +
                                       xr(jlocr  ,klocr+1,irn) +
                                       xr(jlocr+1,klocr,irn) +
                                       xr(jlocr+1,klocr+1,irn) );
                    xyzn(ij,2) = .25*( yr(jlocr  ,klocr,irn) +
                                       yr(jlocr  ,klocr+1,irn) +
                                       yr(jlocr+1,klocr,irn) +
                                       yr(jlocr+1,klocr+1,irn) );
                    xyzn(ij,3) = .25*( zr(jlocr  ,klocr,irn) +
                                       zr(jlocr  ,klocr+1,irn) +
                                       zr(jlocr+1,klocr,irn) +
                                       zr(jlocr+1,klocr+1,irn) );
                }
            }
        }
    }

    return;
}

} // namespace pre_cblki_ns
