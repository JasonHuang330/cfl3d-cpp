// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "setslave.h"
#include "lead.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace setslave_ns {

void setslave(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, int& maxbl, int& mxbli, int& maxgr, int& maxseg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nblock, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf, FortranArray3DRef<int> islavept, int& nslave, FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg, int& maxsegdg, FortranArray1DRef<int> iwk, int& iwork, int& nmaster, int& ngrid, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip, FortranArray2DRef<int> iskip, FortranArray1DRef<int> nblg, FortranArray1DRef<int> levelg, int& lfgm, FortranArray2DRef<int> nblk, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray2DRef<int> nblelst, int& nnodes, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax, int& nbli)
{
    // COMMON block references
    int& jdim     = cmn_ginfo.jdim;
    int& kdim     = cmn_ginfo.kdim;
    int& idim     = cmn_ginfo.idim;
    int& mseq     = cmn_info.mseq;
    int* levelt   = cmn_info.levelt;  // levelt(5), 0-based in struct
    int& lglobal  = cmn_mgrd.lglobal;
    int& iunst    = cmn_unst.iunst;
    int& idef_ss  = cmn_elastic_ss.idef_ss;
    int& i2d      = cmn_twod.i2d;
    int& isktyp   = cmn_deformz.isktyp;
    int& meshdef  = cmn_deformz.meshdef;

    // Local variables
    // nsgst(0:maxbl,2) — lower bound 0
    FortranArray2D<int> nsgst;
    nsgst.allocate_range(0, maxbl, 1, 2);

    // n11i(20,nslave), iimax(nslave), iimx1(nslave)
    FortranArray2D<int> n11i(20, nslave);
    FortranArray1D<int> iimax(nslave);
    FortranArray1D<int> iimx1(nslave);

    // Fixed-size local arrays
    // nsgst1(2) declared in Fortran but unused
    int nskp1[4];   // nskp1(3) — 1-based, indices 1..3
    int nskp[4];    // nskp(3)
    int in_arr[4];  // in(3)
    int im_arr[4];  // im(3)
    int ig_arr[4];  // ig(3)
    int id_arr[4];  // id(3)
    int iss[4][3];  // iss(3,2) — Fortran col-major: iss(nn,jj) → iss[jj-1][nn-1]
    int ise[4][3];  // ise(3,2)
    int ijktot[4][3]; // ijktot(2,3) — ijktot(ii,jj) → ijktot[jj-1][ii-1]

    // Local scalars
    int iseq, n, nbl, nn, nbl2, ivert, i, j, k, ll;
    int iskp, jskp, kskp;
    int jm1, km1, im1, jp1, kp1, ip1;
    int nseg, ista, iend, jsta, jend, ksta, kend;
    int i1, j1, k1, iwrap, jwrap, kwrap;
    int itot, jtot, ktot;
    int m, nbl1;
    int nsgst12, nsgst11, nsgst22, nsgst21;
    int idim1, jdim1, kdim1, idim2, jdim2, kdim2;
    int isva21, isva22, isva11, isva12, isva23, isva13;
    int lmb2210, lmb2213, lmb2220, lmb2223, lmb2230, lmb2233;
    int lmb1110, lmb1113, lmb1120, lmb1123, lmb1130, lmb1133;
    int iskp0, jskp0, kskp0, iskp2, jskp2, kskp2;
    int inc1, inc2, jj1, jj2;
    int i2, j2, k2, ll1, ll2;
    int ii, jj, ii3, ii4, ii5;
    int ipt1, ipt2, ni3;
    int n2;
    int issim11, iseim11, issim12, iseim12;
    int issim21, iseim21, issim22, iseim22;
    int issim31, iseim31, issim32, iseim32;
    int neg1_val;
    (void)neg1_val;



    // islavept = 0
    for (int _i3 = islavept.lbound(3); _i3 <= islavept.ubound(3); _i3++)
        for (int _i2 = islavept.lbound(2); _i2 <= islavept.ubound(2); _i2++)
            for (int _i1 = islavept.lbound(1); _i1 <= islavept.ubound(1); _i1++)
                islavept(_i1, _i2, _i3) = 0;

    if (iunst > 1 || idef_ss > 0) {
      for (iseq = 1; iseq <= mseq; iseq++) {
        for (n = 1; n <= nslave; n++) {
            iimax(n) = 1;
        }

        lglobal = lfgm - (mseq - iseq);

        for (nbl = 0; nbl <= nblock; nbl++) {
            for (nn = 1; nn <= 2; nn++) {
                nsgst(nbl, nn) = 0;
            }
        }
        ivert = 0;
        nbl2  = 1;

        for (nbl = 1; nbl <= nblock; nbl++) {
            if ((levelg(nbl) >= lglobal &&
                 levelg(nbl) <= levelt[iseq-1])) {
                lead_ns::lead(nbl, lw, lw2, maxbl);
                iskp = iskip(nbl, 1);
                jskp = jskip(nbl, 1);
                kskp = kskip(nbl, 1);
                if (std::abs(isktyp) == 1) {
                    iskmax(nbl) = (idim - 1) / iskp + 1;
                    jskmax(nbl) = (jdim - 1) / jskp + 1;
                    kskmax(nbl) = (kdim - 1) / kskp + 1;
                    nsgst(nbl2, 1) = nsgst(nbl2 - 1, 1) +
                        (jdim + jskp - 1) * (kdim + kskp - 1) * (idim + iskp - 1)
                        / jskp / kskp / iskp;
                    nsgst(nbl2, 2) = nbl;
                    nbl2 = nbl2 + 1;
                    for (i = 1; i <= idim; i += iskp) {
                        for (k = 1; k <= kdim; k += kskp) {
                            for (j = 1; j <= jdim; j += jskp) {
                                ll = j - 1 + jdim * (k - 1) + jdim * kdim * (i - 1);
                                ivert = ivert + 1;
                                if (ivert > nslave) {
                                    nou(1) = std::min(nou(1) + 1, ibufdim);
                                    std::snprintf(bou(nou(1), 1), 120,
                                        "stopping...parameter nslave too small");
                                    termn8_ns::termn8(myid, neg1_val = -1, ibufdim, nbuf, bou, nou);
                                }
                                islavept(ivert, 1, iseq) = ll;
                                islavept(ivert, 8, iseq) = 1;
                                islavept(ivert, 9, iseq) = nbl;
                                jm1 = -1;
                                km1 = -(jdim + jskp - 1) / jskp;
                                im1 = -(jdim + jskp - 1) * (kdim + kskp - 1) / jskp / kskp;
                                jp1 = 1;
                                kp1 = (jdim + jskp - 1) / jskp;
                                ip1 = (jdim + jskp - 1) * (kdim + kskp - 1) / jskp / kskp;
                                if (j == 1)    jm1 = 0;
                                if (j == jdim) jp1 = 0;
                                if (k == 1)    km1 = 0;
                                if (k == kdim) kp1 = 0;
                                if (i == 1)    im1 = 0;
                                if (i == idim) ip1 = 0;
                                if (i2d != 0) goto label1000_a;
                                for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
                                    ista = 1;
                                    iend = 1;
                                    jsta = ibcinfo(nbl, nseg, 2, 1);
                                    jend = ibcinfo(nbl, nseg, 3, 1);
                                    ksta = ibcinfo(nbl, nseg, 4, 1);
                                    kend = ibcinfo(nbl, nseg, 5, 1);
                                    if (ibcinfo(nbl, nseg, 1, 1) != 2005) {
                                        if (ibcinfo(nbl, nseg, 1, 1) == 1005 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 1006 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2004 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2014 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2024 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2034 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2016 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 1002 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 1000 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 1003 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 1008 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2003 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2006 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2007 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2008 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2102 ||
                                            ibcinfo(nbl, nseg, 1, 1) == 2103) {
                                            if (i >= ista && i <= iend && j >= jsta &&
                                                j <= jend && k >= ksta && k <= kend) {
                                                islavept(ivert, 8, iseq) = 0;
                                            }
                                        }
                                    }
                                }
                                for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
                                    ista = idim;
                                    iend = idim;
                                    jsta = ibcinfo(nbl, nseg, 2, 2);
                                    jend = ibcinfo(nbl, nseg, 3, 2);
                                    ksta = ibcinfo(nbl, nseg, 4, 2);
                                    kend = ibcinfo(nbl, nseg, 5, 2);
                                    if (ibcinfo(nbl, nseg, 1, 2) != 2005) {
                                        if (ibcinfo(nbl, nseg, 1, 2) == 1005 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 1006 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2004 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2014 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2024 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2034 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2016 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 1002 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 1000 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 1003 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 1008 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2003 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2006 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2007 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2008 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2102 ||
                                            ibcinfo(nbl, nseg, 1, 2) == 2103) {
                                            if (i >= ista && i <= iend && j >= jsta &&
                                                j <= jend && k >= ksta && k <= kend) {
                                                islavept(ivert, 8, iseq) = 0;
                                            }
                                        }
                                    }
                                }
                                label1000_a:;

                                // j-face BC checks (nbcj0)
                                for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
                                    ista = jbcinfo(nbl, nseg, 2, 1);
                                    iend = jbcinfo(nbl, nseg, 3, 1);
                                    jsta = 1;
                                    jend = 1;
                                    ksta = jbcinfo(nbl, nseg, 4, 1);
                                    kend = jbcinfo(nbl, nseg, 5, 1);
                                    if (jbcinfo(nbl, nseg, 1, 1) != 2005) {
                                        if (jbcinfo(nbl, nseg, 1, 1) == 1005 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 1006 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2004 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2014 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2024 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2034 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2016 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 1002 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 1000 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 1003 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 1008 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2003 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2006 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2007 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2008 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2102 ||
                                            jbcinfo(nbl, nseg, 1, 1) == 2103 ||
                                            (jbcinfo(nbl, nseg, 1, 1) == 1013 && idim == 2)) {
                                            if (i >= ista && i <= iend && j >= jsta &&
                                                j <= jend && k >= ksta && k <= kend) {
                                                islavept(ivert, 8, iseq) = 0;
                                            }
                                        }
                                    }
                                }
                                for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
                                    ista = jbcinfo(nbl, nseg, 2, 2);
                                    iend = jbcinfo(nbl, nseg, 3, 2);
                                    jsta = jdim;
                                    jend = jdim;
                                    ksta = jbcinfo(nbl, nseg, 4, 2);
                                    kend = jbcinfo(nbl, nseg, 5, 2);
                                    if (jbcinfo(nbl, nseg, 1, 2) != 2005) {
                                        if (jbcinfo(nbl, nseg, 1, 2) == 1005 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 1006 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2004 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2014 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2024 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2034 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2016 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 1002 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 1000 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 1003 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 1008 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2003 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2006 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2007 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2008 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2102 ||
                                            jbcinfo(nbl, nseg, 1, 2) == 2103 ||
                                            (jbcinfo(nbl, nseg, 1, 2) == 1013 && idim == 2)) {
                                            if (i >= ista && i <= iend && j >= jsta &&
                                                j <= jend && k >= ksta && k <= kend) {
                                                islavept(ivert, 8, iseq) = 0;
                                            }
                                        }
                                    }
                                }

                                // k-face BC checks (nbck0)
                                for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
                                    ista = kbcinfo(nbl, nseg, 2, 1);
                                    iend = kbcinfo(nbl, nseg, 3, 1);
                                    jsta = kbcinfo(nbl, nseg, 4, 1);
                                    jend = kbcinfo(nbl, nseg, 5, 1);
                                    ksta = 1;
                                    kend = 1;
                                    if (kbcinfo(nbl, nseg, 1, 1) != 2005) {
                                        if (kbcinfo(nbl, nseg, 1, 1) == 1005 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 1006 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2004 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2014 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2024 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2034 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2016 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 1002 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 1000 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 1003 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 1008 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2003 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2006 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2007 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2008 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2102 ||
                                            kbcinfo(nbl, nseg, 1, 1) == 2103 ||
                                            (kbcinfo(nbl, nseg, 1, 1) == 1013 && idim == 2)) {
                                            if (i >= ista && i <= iend && j >= jsta &&
                                                j <= jend && k >= ksta && k <= kend) {
                                                islavept(ivert, 8, iseq) = 0;
                                            }
                                        }
                                    }
                                }
                                for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
                                    ista = kbcinfo(nbl, nseg, 2, 2);
                                    iend = kbcinfo(nbl, nseg, 3, 2);
                                    jsta = kbcinfo(nbl, nseg, 4, 2);
                                    jend = kbcinfo(nbl, nseg, 5, 2);
                                    ksta = kdim;
                                    kend = kdim;
                                    if (kbcinfo(nbl, nseg, 1, 2) != 2005) {
                                        if (kbcinfo(nbl, nseg, 1, 2) == 1005 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 1006 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2004 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2014 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2024 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2034 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2016 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 1002 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 1000 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 1003 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 1008 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2003 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2006 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2007 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2008 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2102 ||
                                            kbcinfo(nbl, nseg, 1, 2) == 2103 ||
                                            (kbcinfo(nbl, nseg, 1, 2) == 1013 && idim == 2)) {
                                            if (i >= ista && i <= iend && j >= jsta &&
                                                j <= jend && k >= ksta && k <= kend) {
                                                islavept(ivert, 8, iseq) = 0;
                                            }
                                        }
                                    }
                                }
                                islavept(ivert, 2, iseq) = ivert + jm1;
                                islavept(ivert, 3, iseq) = ivert + jp1;
                                islavept(ivert, 4, iseq) = ivert + km1;
                                islavept(ivert, 5, iseq) = ivert + kp1;
                                islavept(ivert, 6, iseq) = ivert + im1;
                                islavept(ivert, 7, iseq) = ivert + ip1;
                            } // j1 loop
                        } // k1 loop
                    } // i1 loop
                } // end else (isktyp != 1)
            } // end if levelg
        } // end nbl loop (first pass)


        // Second pass: shift nsgst and set up block interface connections
        for (nbl = nblock; nbl >= 1; nbl--) {
            nsgst(nbl, 1) = nsgst(nbl - 1, 1);
        }
        ivert = 0;
        for (m = 1; m <= nbli; m++) {
            nbl1 = nblk(1, m);
            if ((levelg(nbl1) >= lglobal &&
                 levelg(nbl1) <= levelt[iseq-1])) {
                for (nbl = 1; nbl <= nblock; nbl++) {
                    if (nsgst(nbl, 2) == nblk(1, m)) {
                        nsgst12 = nsgst(nbl, 2);
                        nsgst11 = nsgst(nbl, 1);
                    }
                    if (nsgst(nbl, 2) == nblk(2, m)) {
                        nsgst22 = nsgst(nbl, 2);
                        nsgst21 = nsgst(nbl, 1);
                    }
                }
                lead_ns::lead(nsgst12, lw, lw2, maxbl);
                idim1 = idim;
                jdim1 = jdim;
                kdim1 = kdim;
                lead_ns::lead(nsgst22, lw, lw2, maxbl);
                idim2 = idim;
                jdim2 = jdim;
                kdim2 = kdim;

                // Initialize iss and ise
                for (nn = 1; nn <= 3; nn++) {
                    iss[1][nn-1] = 0; ise[1][nn-1] = 0;
                    iss[2][nn-1] = 0; ise[2][nn-1] = 0;
                }
                isva21 = isva(2, 1, m);
                isva22 = isva(2, 2, m);
                isva11 = isva(1, 1, m);
                isva12 = isva(1, 2, m);
                if      (isva(2,1,m)==1 && isva(2,2,m)==2) { isva23 = 3; }
                else if (isva(2,2,m)==1 && isva(2,1,m)==2) { isva23 = 3; }
                else if (isva(2,1,m)==2 && isva(2,2,m)==3) { isva23 = 1; }
                else if (isva(2,2,m)==2 && isva(2,1,m)==3) { isva23 = 1; }
                else if (isva(2,1,m)==1 && isva(2,2,m)==3) { isva23 = 2; }
                else                                        { isva23 = 2; }
                if      (isva(1,1,m)==1 && isva(1,2,m)==2) { isva13 = 3; }
                else if (isva(1,2,m)==1 && isva(1,1,m)==2) { isva13 = 3; }
                else if (isva(1,1,m)==2 && isva(1,2,m)==3) { isva13 = 1; }
                else if (isva(1,2,m)==2 && isva(1,1,m)==3) { isva13 = 1; }
                else if (isva(1,1,m)==1 && isva(1,2,m)==3) { isva13 = 2; }
                else                                        { isva13 = 2; }

                lmb2210 = limblk(2, isva21,   m);
                lmb2213 = limblk(2, isva21+3, m);
                lmb2220 = limblk(2, isva22,   m);
                lmb2223 = limblk(2, isva22+3, m);
                lmb2230 = limblk(2, isva23,   m);
                lmb2233 = limblk(2, isva23+3, m);
                lmb1110 = limblk(1, isva11,   m);
                lmb1113 = limblk(1, isva11+3, m);
                lmb1120 = limblk(1, isva12,   m);
                lmb1123 = limblk(1, isva12+3, m);
                lmb1130 = limblk(1, isva13,   m);
                lmb1133 = limblk(1, isva13+3, m);

                // iss/ise indexing: iss(nn,jj) stored as iss[jj-1][nn-1]
                iss[2-1][isva23-1] = lmb2230;
                ise[2-1][isva23-1] = lmb2233;
                if (lmb2210 <= lmb2213) {
                    if (lmb2220 <= lmb2223) {
                        iss[2-1][isva21-1] = lmb2210;
                        ise[2-1][isva21-1] = lmb2213 + 1;
                        iss[2-1][isva22-1] = lmb2220;
                        ise[2-1][isva22-1] = lmb2223 + 1;
                    } else {
                        iss[2-1][isva21-1] = lmb2210;
                        ise[2-1][isva21-1] = lmb2213 + 1;
                        iss[2-1][isva22-1] = lmb2220 + 1;
                        ise[2-1][isva22-1] = lmb2223;
                    }
                } else {
                    if (lmb2220 <= lmb2223) {
                        iss[2-1][isva21-1] = lmb2210 + 1;
                        ise[2-1][isva21-1] = lmb2213;
                        iss[2-1][isva22-1] = lmb2220;
                        ise[2-1][isva22-1] = lmb2223 + 1;
                    } else {
                        iss[2-1][isva21-1] = lmb2210 + 1;
                        ise[2-1][isva21-1] = lmb2213;
                        iss[2-1][isva22-1] = lmb2220 + 1;
                        ise[2-1][isva22-1] = lmb2223;
                    }
                }
                // limblk(1...) data
                iss[1-1][isva13-1] = lmb1130;
                ise[1-1][isva13-1] = lmb1133;
                if (lmb1110 <= lmb1113) {
                    if (lmb1120 <= lmb1123) {
                        iss[1-1][isva11-1] = lmb1110;
                        ise[1-1][isva11-1] = lmb1113 + 1;
                        iss[1-1][isva12-1] = lmb1120;
                        ise[1-1][isva12-1] = lmb1123 + 1;
                    } else {
                        iss[1-1][isva11-1] = lmb1110;
                        ise[1-1][isva11-1] = lmb1113 + 1;
                        iss[1-1][isva12-1] = lmb1120 + 1;
                        ise[1-1][isva12-1] = lmb1123;
                    }
                } else {
                    if (lmb1120 <= lmb1123) {
                        iss[1-1][isva11-1] = lmb1110 + 1;
                        ise[1-1][isva11-1] = lmb1113;
                        iss[1-1][isva12-1] = lmb1120;
                        ise[1-1][isva12-1] = lmb1123 + 1;
                    } else {
                        iss[1-1][isva11-1] = lmb1110 + 1;
                        ise[1-1][isva11-1] = lmb1113;
                        iss[1-1][isva12-1] = lmb1120 + 1;
                        ise[1-1][isva12-1] = lmb1123;
                    }
                }


                if (std::abs(isktyp) == 1) {
                    // Check iss/ise validity
                    for (ii = 1; ii <= 3; ii++) {
                        for (jj = 1; jj <= 2; jj++) {
                            if (iss[jj-1][ii-1] == 0 || ise[jj-1][ii-1] == 0) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    "stopping...deforming mesh control points do not match blocking end points");
                                termn8_ns::termn8(myid, neg1_val = -1, ibufdim, nbuf, bou, nou);
                            }
                        }
                    }
                    iskp0 = iskip(nsgst12, 1);
                    jskp0 = jskip(nsgst12, 1);
                    kskp0 = kskip(nsgst12, 1);
                    nskp[1] = iskip(nsgst12, 1);
                    nskp[2] = jskip(nsgst12, 1);
                    nskp[3] = kskip(nsgst12, 1);
                    iskp2 = iskip(nsgst22, 1);
                    jskp2 = jskip(nsgst22, 1);
                    kskp2 = kskip(nsgst22, 1);
                    nskp1[1] = iskip(nsgst22, 1);
                    nskp1[2] = jskip(nsgst22, 1);
                    nskp1[3] = kskip(nsgst22, 1);
                    for (nn = 1; nn <= 3; nn++) {
                        if (iss[1-1][nn-1] == ise[1-1][nn-1]) nskp[nn] = 0;
                        if (iss[1-1][nn-1] >  ise[1-1][nn-1]) nskp[nn] = -nskp[nn];
                    }
                    for (nn = 1; nn <= 3; nn++) {
                        if (iss[2-1][nn-1] > ise[2-1][nn-1]) nskp1[nn] = -nskp1[nn];
                    }
                    if      (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==1    ) { inc1= 1;                                                                                    jj1=2; }
                    else if (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==jdim1) { inc1=-1;                                                                                    jj1=3; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==1    ) { inc1= (jdim1+jskp0-1)/jskp0;                                                               jj1=4; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==kdim1) { inc1=-(jdim1+jskp0-1)/jskp0;                                                               jj1=5; }
                    else if (iss[2-1][1-1]==ise[2-1][1-1] && iss[2-1][1-1]==1    ) { inc1= (jdim1+jskp0-1)*(kdim1+kskp0-1)/jskp0/kskp0;                                        jj1=6; }
                    else                                                             { inc1=-(jdim1+jskp0-1)*(kdim1+kskp0-1)/jskp0/kskp0;                                        jj1=7; }
                    if      (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==1    ) { inc2= 1;                                                                                    jj2=2; }
                    else if (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==jdim2) { inc2=-1;                                                                                    jj2=3; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==1    ) { inc2= (jdim2+jskp2-1)/jskp2;                                                               jj2=4; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==kdim2) { inc2=-(jdim2+jskp2-1)/jskp2;                                                               jj2=5; }
                    else if (iss[2-1][1-1]==ise[2-1][1-1] && iss[2-1][1-1]==1    ) { inc2= (jdim2+jskp2-1)*(kdim2+kskp2-1)/jskp2/kskp2;                                        jj2=6; }
                    else                                                             { inc2=-(jdim2+jskp2-1)*(kdim2+kskp2-1)/jskp2/kskp2;                                        jj2=7; }
                    for (nn = 1; nn <= 2; nn++) {
                        in_arr[nn] = isva(1, nn, m);
                        im_arr[nn] = isva(2, nn, m);
                    }
                    in_arr[3] = isva13;
                    im_arr[3] = isva23;
                    id_arr[in_arr[1]] = iss[1-1][in_arr[1]-1];
                    for (i2 = iss[2-1][im_arr[1]-1]; nskp1[im_arr[1]] > 0 ? i2 <= ise[2-1][im_arr[1]-1] : i2 >= ise[2-1][im_arr[1]-1]; i2 += nskp1[im_arr[1]]) {
                        ig_arr[im_arr[1]] = i2;
                        id_arr[in_arr[3]] = iss[1-1][in_arr[3]-1];
                        for (k2 = iss[2-1][im_arr[3]-1]; nskp1[im_arr[3]] > 0 ? k2 <= ise[2-1][im_arr[3]-1] : k2 >= ise[2-1][im_arr[3]-1]; k2 += nskp1[im_arr[3]]) {
                            ig_arr[im_arr[3]] = k2;
                            id_arr[in_arr[2]] = iss[1-1][in_arr[2]-1];
                            for (j2 = iss[2-1][im_arr[2]-1]; nskp1[im_arr[2]] > 0 ? j2 <= ise[2-1][im_arr[2]-1] : j2 >= ise[2-1][im_arr[2]-1]; j2 += nskp1[im_arr[2]]) {
                                ig_arr[im_arr[2]] = j2;
                                ll1 = (id_arr[2]+jskp0-1)/jskp0 +
                                      (jdim1+jskp0-1)*(id_arr[3]-1)/jskp0/kskp0 +
                                      (jdim1+jskp0-1)*(kdim1+kskp0-1)*(id_arr[1]-1)/jskp0/kskp0/iskp0;
                                ll2 = (ig_arr[2]+jskp2-1)/jskp2 +
                                      (jdim2+jskp2-1)*(ig_arr[3]-1)/jskp2/kskp2 +
                                      (jdim2+jskp2-1)*(kdim2+kskp2-1)*(ig_arr[1]-1)/jskp2/kskp2/iskp2;
                                islavept(nsgst11+ll1, jj1, iseq) = nsgst21+ll2+inc2;
                                if (nsgst11+ll1 != nsgst21+ll2) {
                                    iimax(nsgst11+ll1) = iimax(nsgst11+ll1) + 1;
                                    ii3 = 10 + iimax(nsgst11+ll1);
                                    islavept(nsgst11+ll1, ii3, iseq) = nsgst21+ll2;
                                    iimax(nsgst21+ll2) = iimax(nsgst21+ll2) + 1;
                                    ii3 = 10 + iimax(nsgst21+ll2);
                                    islavept(nsgst21+ll2, ii3, iseq) = nsgst11+ll1;
                                }
                                islavept(nsgst21+ll2, jj2, iseq) = nsgst11+ll1+inc1;
                                id_arr[in_arr[2]] = id_arr[in_arr[2]] + nskp[in_arr[2]];
                            }
                            id_arr[in_arr[3]] = id_arr[in_arr[3]] + nskp[in_arr[3]];
                        }
                        id_arr[in_arr[1]] = id_arr[in_arr[1]] + nskp[in_arr[1]];
                    }
                } else {

                    // isktyp != 1 branch
                    for (ii = 1; ii <= 3; ii++) {
                        for (jj = 1; jj <= 2; jj++) {
                            if (iss[jj-1][ii-1] == 0 || ise[jj-1][ii-1] == 0) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    "stopping...deforming mesh control points do not match blocking end points");
                                termn8_ns::termn8(myid, neg1_val = -1, ibufdim, nbuf, bou, nou);
                            }
                        }
                    }
                    for (nn = 1; nn <= 3; nn++) {
                        nskp[nn] = 1;
                        if (iss[1-1][nn-1] == ise[1-1][nn-1]) nskp[nn] = 0;
                        if (iss[1-1][nn-1] >  ise[1-1][nn-1]) nskp[nn] = -1;
                    }
                    for (nn = 1; nn <= 3; nn++) {
                        nskp1[nn] = 1;
                        if (iss[2-1][nn-1] > ise[2-1][nn-1]) nskp1[nn] = -1;
                    }
                    // ijktot(ii,jj) stored as ijktot[jj-1][ii-1]
                    ijktot[1-1][1-1] = iskmax(nsgst12);
                    ijktot[1-1][2-1] = jskmax(nsgst12);
                    ijktot[1-1][3-1] = kskmax(nsgst12);
                    ijktot[2-1][1-1] = iskmax(nsgst22);
                    ijktot[2-1][2-1] = jskmax(nsgst22);
                    ijktot[2-1][3-1] = kskmax(nsgst22);
                    if      (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==1    ) { inc1= 1;                                    jj1=2; }
                    else if (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==jdim1) { inc1=-1;                                    jj1=3; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==1    ) { inc1= ijktot[1-1][2-1];                     jj1=4; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==kdim1) { inc1=-ijktot[1-1][2-1];                     jj1=5; }
                    else if (iss[2-1][1-1]==ise[2-1][1-1] && iss[2-1][1-1]==1    ) { inc1= ijktot[1-1][2-1]*ijktot[1-1][3-1];   jj1=6; }
                    else                                                             { inc1=-ijktot[1-1][2-1]*ijktot[1-1][3-1];   jj1=7; }
                    if      (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==1    ) { inc2= 1;                                    jj2=2; }
                    else if (iss[2-1][2-1]==ise[2-1][2-1] && iss[2-1][2-1]==jdim2) { inc2=-1;                                    jj2=3; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==1    ) { inc2= ijktot[2-1][2-1];                     jj2=4; }
                    else if (iss[2-1][3-1]==ise[2-1][3-1] && iss[2-1][3-1]==kdim2) { inc2=-ijktot[2-1][2-1];                     jj2=5; }
                    else if (iss[2-1][1-1]==ise[2-1][1-1] && iss[2-1][1-1]==1    ) { inc2= ijktot[2-1][2-1]*ijktot[2-1][3-1];   jj2=6; }
                    else                                                             { inc2=-ijktot[2-1][2-1]*ijktot[2-1][3-1];   jj2=7; }
                    for (nn = 1; nn <= 2; nn++) {
                        in_arr[nn] = isva(1, nn, m);
                        im_arr[nn] = isva(2, nn, m);
                    }
                    in_arr[3] = isva13;
                    im_arr[3] = isva23;
                    // Find index positions in iskip/jskip/kskip arrays
                    issim11 = 0; iseim11 = 0;
                    issim12 = 0; iseim12 = 0;
                    issim21 = 0; iseim21 = 0;
                    issim22 = 0; iseim22 = 0;
                    issim31 = 0; iseim31 = 0;
                    issim32 = 0; iseim32 = 0;
                    for (ii = 1; ii <= ijktot[1-1][1-1]; ii++) {
                        if (iss[1-1][1-1] == iskip(nsgst12, ii)) issim11 = ii;
                        if (ise[1-1][1-1] == iskip(nsgst12, ii)) iseim11 = ii;
                    }
                    for (ii = 1; ii <= ijktot[2-1][1-1]; ii++) {
                        if (iss[2-1][1-1] == iskip(nsgst22, ii)) issim12 = ii;
                        if (ise[2-1][1-1] == iskip(nsgst22, ii)) iseim12 = ii;
                    }
                    for (ii = 1; ii <= ijktot[1-1][2-1]; ii++) {
                        if (iss[1-1][2-1] == jskip(nsgst12, ii)) issim21 = ii;
                        if (ise[1-1][2-1] == jskip(nsgst12, ii)) iseim21 = ii;
                    }
                    for (ii = 1; ii <= ijktot[2-1][2-1]; ii++) {
                        if (iss[2-1][2-1] == jskip(nsgst22, ii)) issim22 = ii;
                        if (ise[2-1][2-1] == jskip(nsgst22, ii)) iseim22 = ii;
                    }
                    for (ii = 1; ii <= ijktot[1-1][3-1]; ii++) {
                        if (iss[1-1][3-1] == kskip(nsgst12, ii)) issim31 = ii;
                        if (ise[1-1][3-1] == kskip(nsgst12, ii)) iseim31 = ii;
                    }
                    for (ii = 1; ii <= ijktot[2-1][3-1]; ii++) {
                        if (iss[2-1][3-1] == kskip(nsgst22, ii)) issim32 = ii;
                        if (ise[2-1][3-1] == kskip(nsgst22, ii)) iseim32 = ii;
                    }
                    iss[1-1][1-1] = issim11; ise[1-1][1-1] = iseim11;
                    iss[1-1][2-1] = issim21; ise[1-1][2-1] = iseim21;
                    iss[1-1][3-1] = issim31; ise[1-1][3-1] = iseim31;
                    iss[2-1][1-1] = issim12; ise[2-1][1-1] = iseim12;
                    iss[2-1][2-1] = issim22; ise[2-1][2-1] = iseim22;
                    iss[2-1][3-1] = issim32; ise[2-1][3-1] = iseim32;
                    for (ii = 1; ii <= 3; ii++) {
                        for (jj = 1; jj <= 2; jj++) {
                            if (iss[jj-1][ii-1] == 0 || ise[jj-1][ii-1] == 0) {
                                nou(1) = std::min(nou(1) + 1, ibufdim);
                                std::snprintf(bou(nou(1), 1), 120,
                                    "stopping...deforming mesh control points do not match blocking end points");
                                termn8_ns::termn8(myid, neg1_val = -1, ibufdim, nbuf, bou, nou);
                            }
                        }
                    }


                    id_arr[in_arr[1]] = iss[1-1][in_arr[1]-1];
                    for (i2 = iss[2-1][im_arr[1]-1]; nskp1[im_arr[1]] > 0 ? i2 <= ise[2-1][im_arr[1]-1] : i2 >= ise[2-1][im_arr[1]-1]; i2 += nskp1[im_arr[1]]) {
                        ig_arr[im_arr[1]] = i2;
                        id_arr[in_arr[3]] = iss[1-1][in_arr[3]-1];
                        for (k2 = iss[2-1][im_arr[3]-1]; nskp1[im_arr[3]] > 0 ? k2 <= ise[2-1][im_arr[3]-1] : k2 >= ise[2-1][im_arr[3]-1]; k2 += nskp1[im_arr[3]]) {
                            ig_arr[im_arr[3]] = k2;
                            id_arr[in_arr[2]] = iss[1-1][in_arr[2]-1];
                            for (j2 = iss[2-1][im_arr[2]-1]; nskp1[im_arr[2]] > 0 ? j2 <= ise[2-1][im_arr[2]-1] : j2 >= ise[2-1][im_arr[2]-1]; j2 += nskp1[im_arr[2]]) {
                                ig_arr[im_arr[2]] = j2;
                                ll1 = id_arr[2] + (id_arr[3]-1)*ijktot[1-1][2-1] + (id_arr[1]-1)*ijktot[1-1][2-1]*ijktot[1-1][3-1];
                                ll2 = ig_arr[2] + (ig_arr[3]-1)*ijktot[2-1][2-1] + (ig_arr[1]-1)*ijktot[2-1][2-1]*ijktot[2-1][3-1];
                                islavept(nsgst11+ll1, jj1, iseq) = nsgst21+ll2+inc2;
                                if (nsgst11+ll1 != nsgst21+ll2) {
                                    iimax(nsgst11+ll1) = iimax(nsgst11+ll1) + 1;
                                    ii3 = 10 + iimax(nsgst11+ll1);
                                    islavept(nsgst11+ll1, ii3, iseq) = nsgst21+ll2;
                                    iimax(nsgst21+ll2) = iimax(nsgst21+ll2) + 1;
                                    ii3 = 10 + iimax(nsgst21+ll2);
                                    islavept(nsgst21+ll2, ii3, iseq) = nsgst11+ll1;
                                }
                                islavept(nsgst21+ll2, jj2, iseq) = nsgst11+ll1+inc1;
                                id_arr[in_arr[2]] = id_arr[in_arr[2]] + nskp[in_arr[2]];
                            }
                            id_arr[in_arr[3]] = id_arr[in_arr[3]] + nskp[in_arr[3]];
                        }
                        id_arr[in_arr[1]] = id_arr[in_arr[1]] + nskp[in_arr[1]];
                    }
                } // end else (isktyp != 1) in second pass
            } // end if levelg
        } // end m loop (second pass)

        // Determine number of coincident points at block interfaces
        for (n = 1; n <= nslave; n++) {
            if (iimax(n) == 1) iimax(n) = 0;
            islavept(n, 11, iseq) = iimax(n);
            if (islavept(n, 11, iseq) > 0) {
                n11i(1, n) = n;
                for (ii3 = 2; ii3 <= iimax(n); ii3++) {
                    n11i(ii3, n) = islavept(n, 12+ii3-2, iseq);
                }
            }
        }
        for (n = 1; n <= nslave; n++) {
            if (islavept(n, 11, iseq) > 0) {
                iimx1(n) = iimax(n);
                for (n2 = 1; n2 <= nslave; n2++) {
                    if (n2 != n) {
                        for (ii3 = 1; ii3 <= iimax(n); ii3++) {
                            for (ii4 = 1; ii4 <= iimax(n2)-1; ii4++) {
                                if (n11i(ii3, n) == n11i(ii4, n2)) {
                                    for (ii5 = 2; ii5 <= iimax(n); ii5++) {
                                        if (n2 == n11i(ii5, n)) goto label2550;
                                    }
                                    iimx1(n) = iimx1(n) + 1;
                                    n11i(iimx1(n), n) = n2;
                                    goto label2600;
                                }
                                label2550:;
                                if (n11i(ii3, n) == n2) {
                                    for (ii5 = 2; ii5 <= iimax(n); ii5++) {
                                        if (n11i(ii4, n2) == n11i(ii5, n)) goto label2600;
                                    }
                                    iimx1(n) = iimx1(n) + 1;
                                    n11i(iimx1(n), n) = n11i(ii4, n2);
                                    goto label2600;
                                }
                            }
                        }
                        label2600:;
                        islavept(n, 11, iseq) = iimx1(n);
                        if (iimx1(n) > iimax(n) && iimax(n) > 0) {
                            ii4 = iimax(n) + 1;
                            for (ii = 10+iimax(n)+1; ii <= 10+iimx1(n); ii++) {
                                islavept(n, ii, iseq) = n11i(ii4, n);
                                ii4 = ii4 + 1;
                            }
                        }
                    }
                }
                // label2650:;  (unreferenced in Fortran, just continue)
            }
        }


        // Establish the number of points in the ija and sa arrays for each node point
        for (n = 1; n <= nslave; n++) {
            ipt1 = 0;
            ipt2 = 0;
            ii4  = 1;
            if (islavept(n, 8, iseq) == 0) ii4 = 0;
            if (ii4 != 0) {
                iimax(n) = islavept(n, 11, iseq);
                for (int ii2 = 2; ii2 <= iimax(n); ii2++) {
                    ni3 = islavept(n, 12+ii2-2, iseq);
                    if (islavept(ni3, 8, iseq) == 0) ii4 = 0;
                }
            }
            if (ii4 != 0) {
                ipt2 = 2;
                if (islavept(n,2,iseq) != n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq) != n) ipt1 = ipt1 + 1;
                if (islavept(n,4,iseq) != n) ipt1 = ipt1 + 1;
                if (islavept(n,5,iseq) != n) ipt1 = ipt1 + 1;
                if (islavept(n,6,iseq) != n) ipt1 = ipt1 + 1;
                if (islavept(n,7,iseq) != n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,4,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,4,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,4,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,4,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,4,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,4,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,4,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,4,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,5,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,5,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,5,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,5,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,2,iseq)!=n && islavept(n,5,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,5,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,5,iseq)!=n && islavept(n,6,iseq)!=n) ipt1 = ipt1 + 1;
                if (islavept(n,3,iseq)!=n && islavept(n,5,iseq)!=n && islavept(n,7,iseq)!=n) ipt1 = ipt1 + 1;
            }
            islavept(n, 10, iseq) = 3*ipt1 + ipt2;
        }

        // Build nblelst
        for (int _i1 = nblelst.lbound(1); _i1 <= nblelst.ubound(1); _i1++)
            for (int _i2 = nblelst.lbound(2); _i2 <= nblelst.ubound(2); _i2++)
                nblelst(_i1, _i2) = 0;
        ii = 1;
        nblelst(ii, 1) = islavept(1, 9, iseq);
        nblelst(ii, 2) = 1;
        for (n = 2; n <= nslave; n++) {
            if (islavept(n, 9, iseq) != islavept(n-1, 9, iseq)) {
                ii = ii + 1;
                nblelst(ii, 1) = islavept(n, 9, iseq);
                nblelst(ii, 2) = n;
            }
        }
        nblelst(ii+1, 2) = nslave + 1;

      } // end iseq loop
    } // end if (iunst > 1 || idef_ss > 0)

    if (myid == myhost) {
        if (isktyp > 0 && meshdef == 1) {
            FILE* funit = fortran_get_unit(1000 + myid);
            fprintf(funit, " Iteration   Total Resid   N@Maxres  Max Resid\n");
        } else if (isktyp > 0 && meshdef == 0) {
            FILE* funit = fortran_get_unit(1000 + myid);
            fprintf(funit, " Iterations  Total Resid\n");
        }
    }

    return;
}

} // namespace setslave_ns
