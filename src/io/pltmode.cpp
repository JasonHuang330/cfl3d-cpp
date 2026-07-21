// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "pltmode.h"
#include "lead.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

namespace pltmode_ns {

void pltmode(int& nm, int& iaes, FortranArray1DRef<int> nblg, int& ngrid, int& maxgr, int& maxbl,
             FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> iaesurf,
             FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
             FortranArray4DRef<int> ibcinfo, FortranArray1DRef<int> nbcj0,
             FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbck0,
             FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> nbci0,
             FortranArray1DRef<int> nbcidim, int& maxseg, int& maxsegdg,
             FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
             FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk,
             FortranArray5DRef<double> xmdi, int& jdim1, int& kdim1, int& idim1,
             int& maxaes, int& nmds, int& nbl1, int& iunitw, int& iflag)
{
    // COMMON block references
    int& jdim   = cmn_ginfo.jdim;
    int& kdim   = cmn_ginfo.kdim;
    int& idim   = cmn_ginfo.idim;
    int& ialph  = cmn_igrdtyp.ialph;

    // Local variables
    char titlemode[22];
    int jdm[501], kdm[501], idm[501];  // 1-based access via [1..500]

    int len1;
    int i, j, k, m, ns, igrid, nbl, iaesrf, ngd, l, nseg, lbc, iflg;
    int ist, ifn, kst, kfn, jst, jfn;

    // Build filename titlemode based on iaes and nm
    if (iaes > 9) {
        if (nm > 99) {
            len1 = 20;
            std::snprintf(titlemode, sizeof(titlemode), "aesurf%d_mode%d.plt", iaes, nm);
        } else if (nm > 9) {
            len1 = 19;
            std::snprintf(titlemode, sizeof(titlemode), "aesurf%d_mode%d.plt", iaes, nm);
        } else {
            len1 = 18;
            std::snprintf(titlemode, sizeof(titlemode), "aesurf%d_mode%d.plt", iaes, nm);
        }
    } else {
        if (nm > 99) {
            len1 = 19;
            std::snprintf(titlemode, sizeof(titlemode), "aesurf%d_mode%d.plt", iaes, nm);
        } else if (nm > 9) {
            len1 = 18;
            std::snprintf(titlemode, sizeof(titlemode), "aesurf%d_mode%d.plt", iaes, nm);
        } else {
            len1 = 17;
            std::snprintf(titlemode, sizeof(titlemode), "aesurf%d_mode%d.plt", iaes, nm);
        }
    }
    // Pad remaining characters with spaces (titlemode is 21 chars in Fortran)
    // The snprintf already null-terminates; we just need the first len1 chars
    // Build a null-terminated filename string of exactly len1 chars
    char filename[22];
    std::memcpy(filename, titlemode, len1);
    filename[len1] = '\0';

    // open(iunitw,file=titlemode(1:len1),form='formatted',status='unknown')
    fortran_open_unit(iunitw, filename, "w");

    // write(iunitw,23110)
    // format(' title="tecplot modeshape file"',/,' variables="x","y","z","phix","phiy","phiz"')
    FILE* funit = fortran_get_unit(iunitw);
    fprintf(funit, " title=\"tecplot modeshape file\"\n variables=\"x\",\"y\",\"z\",\"phix\",\"phiy\",\"phiz\"\n");

    // format 23111: ' ZONE I = ',i5,' J = ',i5
    // (used later in the iflag>0 branch)

    if (iflag == 0) {
        // get number of zones and zone dimensions for plot3d header
        // rewind(iunitw)
        rewind(funit);

        ngd = 0;
        for (igrid = 1; igrid <= ngrid; igrid++) {
            nbl = nblg(igrid);
            iaesrf = 0;
            for (int is2 = 1; is2 <= nsegdfrm(nbl); is2++) {
                iaesrf = iaesrf + iaesurf(nbl, is2);
            }
            if (iaesrf != 0) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                // j-direction faces (m=1: j=1, m=2: j=jdim)
                for (m = 1; m <= 2; m++) {
                    if (m == 1) {
                        j    = 1;
                        l    = 0;
                        nseg = nbcj0(nbl);
                    } else {
                        j    = jdim;
                        l    = 3;
                        nseg = nbcjdim(nbl);
                    }
                    for (ns = 1; ns <= nseg; ns++) {
                        lbc  = std::abs(jbcinfo(nbl, ns, 1, m));
                        iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
                                         *(lbc-2024)*(lbc-2034)*(lbc-2016);
                        if (iflg == 0) {
                            ngd++;
                            ist = jbcinfo(nbl, ns, 2, m);
                            ifn = jbcinfo(nbl, ns, 3, m);
                            kst = jbcinfo(nbl, ns, 4, m);
                            kfn = jbcinfo(nbl, ns, 5, m);
                            idm[ngd] = ifn - ist + 1;
                            kdm[ngd] = kfn - kst + 1;
                            jdm[ngd] = 1;
                        }
                    }
                }
                // k-direction faces (m=1: k=1, m=2: k=kdim)
                for (m = 1; m <= 2; m++) {
                    if (m == 1) {
                        k    = 1;
                        l    = 0;
                        nseg = nbck0(nbl);
                    } else {
                        k    = kdim;
                        l    = 3;
                        nseg = nbckdim(nbl);
                    }
                    for (ns = 1; ns <= nseg; ns++) {
                        lbc  = std::abs(kbcinfo(nbl, ns, 1, m));
                        iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
                                         *(lbc-2024)*(lbc-2034)*(lbc-2016);
                        if (iflg == 0) {
                            ngd++;
                            ist = kbcinfo(nbl, ns, 2, m);
                            ifn = kbcinfo(nbl, ns, 3, m);
                            jst = kbcinfo(nbl, ns, 4, m);
                            jfn = kbcinfo(nbl, ns, 5, m);
                            idm[ngd] = ifn - ist + 1;
                            kdm[ngd] = 1;
                            jdm[ngd] = jfn - jst + 1;
                        }
                    }
                }
                // i-direction faces (m=1: i=1, m=2: i=idim)
                for (m = 1; m <= 2; m++) {
                    if (m == 1) {
                        i    = 1;
                        l    = 0;
                        nseg = nbci0(nbl);
                    } else {
                        i    = idim;
                        l    = 3;
                        nseg = nbcidim(nbl);
                    }
                    for (ns = 1; ns <= nseg; ns++) {
                        lbc  = std::abs(ibcinfo(nbl, ns, 1, m));
                        iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
                                         *(lbc-2024)*(lbc-2034)*(lbc-2016);
                        if (iflg == 0) {
                            ngd++;
                            jst = ibcinfo(nbl, ns, 2, m);
                            jfn = ibcinfo(nbl, ns, 3, m);
                            kst = ibcinfo(nbl, ns, 4, m);
                            kfn = ibcinfo(nbl, ns, 5, m);
                            idm[ngd] = 1;
                            kdm[ngd] = kfn - kst + 1;
                            jdm[ngd] = jfn - jst + 1;
                        }
                    }
                }
            }
        }

    } else if (iflag > 0) {
        // write out x,y,z of modal surface
        ngd = 0;

        // j-direction faces
        for (m = 1; m <= 2; m++) {
            if (m == 1) {
                j    = 1;
                l    = 0;
                nseg = nbcj0(nbl1);
            } else {
                j    = jdim1;
                l    = 3;
                nseg = nbcjdim(nbl1);
            }
            for (ns = 1; ns <= nseg; ns++) {
                lbc  = std::abs(jbcinfo(nbl1, ns, 1, m));
                iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
                                  *(lbc-2024)*(lbc-2034)*(lbc-2016);
                if (iflg == 0) {
                    ngd++;
                    ist = jbcinfo(nbl1, ns, 2, m);
                    ifn = jbcinfo(nbl1, ns, 3, m);
                    kst = jbcinfo(nbl1, ns, 4, m);
                    kfn = jbcinfo(nbl1, ns, 5, m);
                    // write(iunitw,23111) ifn-ist+1,kfn-kst+1
                    fprintf(funit, " ZONE I = %5d J = %5d\n", ifn-ist+1, kfn-kst+1);
                    if (ialph == 0) {
                        for (k = kst; k <= kfn; k++) {
                            for (i = ist; i <= ifn; i++) {
                                fprintf(funit, " %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                                    x(j,k,i), y(j,k,i), z(j,k,i),
                                    (float)xmdj(k,i,l+1,nm,iaes),
                                    (float)xmdj(k,i,l+2,nm,iaes),
                                    (float)xmdj(k,i,l+3,nm,iaes));
                            }
                        }
                    } else {
                        for (k = kst; k <= kfn; k++) {
                            for (i = ist; i <= ifn; i++) {
                                fprintf(funit, " %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                                    x(j,k,i), z(j,k,i), -y(j,k,i),
                                    (float)xmdj(k,i,l+1,nm,iaes),
                                    (float)xmdj(k,i,l+3,nm,iaes),
                                    -(float)xmdj(k,i,l+2,nm,iaes));
                            }
                        }
                    }
                }
            }
        }

        // k-direction faces
        for (m = 1; m <= 2; m++) {
            if (m == 1) {
                k    = 1;
                l    = 0;
                nseg = nbck0(nbl1);
            } else {
                k    = kdim1;
                l    = 3;
                nseg = nbckdim(nbl1);
            }
            for (ns = 1; ns <= nseg; ns++) {
                lbc  = std::abs(kbcinfo(nbl1, ns, 1, m));
                iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
                                  *(lbc-2024)*(lbc-2034)*(lbc-2016);
                if (iflg == 0) {
                    ngd++;
                    ist = kbcinfo(nbl1, ns, 2, m);
                    ifn = kbcinfo(nbl1, ns, 3, m);
                    jst = kbcinfo(nbl1, ns, 4, m);
                    jfn = kbcinfo(nbl1, ns, 5, m);
                    // write(iunitw,23111) ifn-ist+1,jfn-jst+1
                    fprintf(funit, " ZONE I = %5d J = %5d\n", ifn-ist+1, jfn-jst+1);
                    if (ialph == 0) {
                        for (j = jst; j <= jfn; j++) {
                            for (i = ist; i <= ifn; i++) {
                                fprintf(funit, " %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                                    x(j,k,i), y(j,k,i), z(j,k,i),
                                    (float)xmdk(j,i,l+1,nm,iaes),
                                    (float)xmdk(j,i,l+2,nm,iaes),
                                    (float)xmdk(j,i,l+3,nm,iaes));
                            }
                        }
                    } else {
                        for (j = jst; j <= jfn; j++) {
                            for (i = ist; i <= ifn; i++) {
                                fprintf(funit, " %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                                    x(j,k,i), z(j,k,i), -y(j,k,i),
                                    (float)xmdk(j,i,l+1,nm,iaes),
                                    (float)xmdk(j,i,l+3,nm,iaes),
                                    -(float)xmdk(j,i,l+2,nm,iaes));
                            }
                        }
                    }
                }
            }
        }

        // i-direction faces
        for (m = 1; m <= 2; m++) {
            if (m == 1) {
                i    = 1;
                l    = 0;
                nseg = nbci0(nbl1);
            } else {
                i    = idim1;
                l    = 3;
                nseg = nbcidim(nbl1);
            }
            for (ns = 1; ns <= nseg; ns++) {
                lbc  = std::abs(ibcinfo(nbl1, ns, 1, m));
                iflg = (lbc-1006)*(lbc-1005)*(lbc-2004)*(lbc-2014)
                                  *(lbc-2024)*(lbc-2034)*(lbc-2016);
                if (iflg == 0) {
                    ngd++;
                    jst = ibcinfo(nbl1, ns, 2, m);
                    jfn = ibcinfo(nbl1, ns, 3, m);
                    kst = ibcinfo(nbl1, ns, 4, m);
                    kfn = ibcinfo(nbl1, ns, 5, m);
                    // write(iunitw,23111) jfn-jst+1,kfn-kst+1
                    fprintf(funit, " ZONE I = %5d J = %5d\n", jfn-jst+1, kfn-kst+1);
                    if (ialph == 0) {
                        for (k = kst; k <= kfn; k++) {
                            for (j = jst; j <= jfn; j++) {
                                fprintf(funit, " %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                                    x(j,k,i), y(j,k,i), z(j,k,i),
                                    (float)xmdi(j,k,l+1,nm,iaes),
                                    (float)xmdi(j,k,l+2,nm,iaes),
                                    (float)xmdi(j,k,l+3,nm,iaes));
                            }
                        }
                    } else {
                        for (k = kst; k <= kfn; k++) {
                            for (j = jst; j <= jfn; j++) {
                                fprintf(funit, " %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                                    x(j,k,i), z(j,k,i), -y(j,k,i),
                                    (float)xmdi(j,k,l+1,nm,iaes),
                                    (float)xmdi(j,k,l+2,nm,iaes),
                                    -(float)xmdi(j,k,l+3,nm,iaes));
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

} // namespace pltmode_ns
