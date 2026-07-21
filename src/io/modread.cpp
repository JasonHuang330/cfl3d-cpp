// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "modread.h"
#include "runtime/fortran_io.h"
#include "umalloc.h"
#include <cstdlib>
#include <cmath>

namespace modread_ns {

void modread(int& idim, int& jdim, int& kdim, int& n, int& nbl, int& iunit,
             FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
             FortranArray4DRef<int> ibcinfo,
             FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim,
             FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim,
             FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim,
             int& maxbl, int& maxseg, int& nmds,
             FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk,
             FortranArray5DRef<double> xmdi, int& iaes, int& maxaes)
{
    // COMMON blocks
    int& ip3dgrd = cmn_igrdtyp.ip3dgrd;
    int& ialph   = cmn_igrdtyp.ialph;

    int stats = 0;
    int memuse;

    // Allocate real (float) arrays in which to read modal data
    FortranArray3D<float> xmdir1;
    FortranArray3D<float> xmdjr1;
    FortranArray3D<float> xmdkr1;

    memuse = 0;
    xmdir1.allocate(jdim, kdim, 6);
    stats = 0;
    {
        int sz = jdim * kdim * 6;
        int iflag0 = 0;
        char text[1] = {'x'};
        umalloc_ns::umalloc(sz, iflag0, text, memuse, stats);
    }

    xmdjr1.allocate(kdim, idim, 6);
    stats = 0;
    {
        int sz = kdim * idim * 6;
        int iflag0 = 0;
        char text[1] = {'x'};
        umalloc_ns::umalloc(sz, iflag0, text, memuse, stats);
    }

    xmdkr1.allocate(jdim, idim, 6);
    stats = 0;
    {
        int sz = jdim * idim * 6;
        int iflag0 = 0;
        char text[1] = {'x'};
        umalloc_ns::umalloc(sz, iflag0, text, memuse, stats);
    }

    // -----------------------------------------------------------------------
    // Read in modal deflections for j=const surfaces
    // -----------------------------------------------------------------------
    for (int m = 1; m <= 2; m++) {

        int j_val, l, nseg;
        if (m == 1) {
            j_val = 1;
            l     = 0;
            nseg  = nbcj0(nbl);
        } else {
            j_val = jdim;
            l     = 3;
            nseg  = nbcjdim(nbl);
        }

        for (int ns = 1; ns <= nseg; ns++) {

            int lbc = std::abs(jbcinfo(nbl, ns, 1, m));

            int iflag = (lbc - 1006) * (lbc - 1005) * (lbc - 2004) * (lbc - 2014) *
                        (lbc - 2024) * (lbc - 2034) * (lbc - 2016);

            if (iflag != 0) continue;

            int ist = jbcinfo(nbl, ns, 2, m);
            int ifn = jbcinfo(nbl, ns, 3, m);
            int kst = jbcinfo(nbl, ns, 4, m);
            int kfn = jbcinfo(nbl, ns, 5, m);

            if (ialph == 0) {
                // read(iunit,*)((xmdj(k,i,l+1,n,iaes),xmdj(k,i,l+2,n,iaes),
                //               xmdj(k,i,l+3,n,iaes),k=kst,kfn),i=ist,ifn)
                // Nested implied-DO: outer i=ist..ifn, inner k=kst..kfn, 3 values each
                int nvals = 3 * (kfn - kst + 1) * (ifn - ist + 1);
                FortranArray1D<double> tmp(nvals);
                fortran_read_list_array(iunit, tmp, 1, nvals);
                int idx = 1;
                for (int i = ist; i <= ifn; i++) {
                    for (int k = kst; k <= kfn; k++) {
                        xmdj(k, i, l + 1, n, iaes) = tmp(idx);     idx++;
                        xmdj(k, i, l + 2, n, iaes) = tmp(idx);     idx++;
                        xmdj(k, i, l + 3, n, iaes) = tmp(idx);     idx++;
                    }
                }
            } else {
                // read(iunit,*)((xmdjr1(k,i,l+1),xmdjr1(k,i,l+3),
                //               xmdjr1(k,i,l+2),k=kst,kfn),i=ist,ifn)
                int nvals = 3 * (kfn - kst + 1) * (ifn - ist + 1);
                FortranArray1D<double> tmp(nvals);
                fortran_read_list_array(iunit, tmp, 1, nvals);
                int idx = 1;
                for (int i = ist; i <= ifn; i++) {
                    for (int k = kst; k <= kfn; k++) {
                        xmdjr1(k, i, l + 1) = (float)tmp(idx);     idx++;
                        xmdjr1(k, i, l + 3) = (float)tmp(idx);     idx++;
                        xmdjr1(k, i, l + 2) = (float)tmp(idx);     idx++;
                    }
                }
                for (int i = ist; i <= ifn; i++) {
                    for (int k = kst; k <= kfn; k++) {
                        xmdj(k, i, l + 1, n, iaes) =  (double)xmdjr1(k, i, l + 1);
                        xmdj(k, i, l + 2, n, iaes) = -(double)xmdjr1(k, i, l + 2);
                        xmdj(k, i, l + 3, n, iaes) =  (double)xmdjr1(k, i, l + 3);
                    }
                }
            }

        } // ns loop
    } // m loop

    // -----------------------------------------------------------------------
    // Read in modal deflections for k=const surfaces
    // -----------------------------------------------------------------------
    for (int m = 1; m <= 2; m++) {

        int k_val, l, nseg;
        if (m == 1) {
            k_val = 1;
            l     = 0;
            nseg  = nbck0(nbl);
        } else {
            k_val = kdim;
            l     = 3;
            nseg  = nbckdim(nbl);
        }

        for (int ns = 1; ns <= nseg; ns++) {

            int lbc = std::abs(kbcinfo(nbl, ns, 1, m));

            int iflag = (lbc - 1006) * (lbc - 1005) * (lbc - 2004) * (lbc - 2014) *
                        (lbc - 2024) * (lbc - 2034) * (lbc - 2016);

            if (iflag != 0) continue;

            int ist = kbcinfo(nbl, ns, 2, m);
            int ifn = kbcinfo(nbl, ns, 3, m);
            int jst = kbcinfo(nbl, ns, 4, m);
            int jfn = kbcinfo(nbl, ns, 5, m);

            if (ialph == 0) {
                // read(iunit,*)((xmdk(j,i,l+1,n,iaes),xmdk(j,i,l+2,n,iaes),
                //               xmdk(j,i,l+3,n,iaes),j=jst,jfn),i=ist,ifn)
                int nvals = 3 * (jfn - jst + 1) * (ifn - ist + 1);
                FortranArray1D<double> tmp(nvals);
                fortran_read_list_array(iunit, tmp, 1, nvals);
                int idx = 1;
                for (int i = ist; i <= ifn; i++) {
                    for (int j = jst; j <= jfn; j++) {
                        xmdk(j, i, l + 1, n, iaes) = tmp(idx);     idx++;
                        xmdk(j, i, l + 2, n, iaes) = tmp(idx);     idx++;
                        xmdk(j, i, l + 3, n, iaes) = tmp(idx);     idx++;
                    }
                }
            } else {
                // read(iunit,*)((xmdkr1(j,i,l+1),xmdkr1(j,i,l+3),
                //               xmdkr1(j,i,l+2),j=jst,jfn),i=ist,ifn)
                int nvals = 3 * (jfn - jst + 1) * (ifn - ist + 1);
                FortranArray1D<double> tmp(nvals);
                fortran_read_list_array(iunit, tmp, 1, nvals);
                int idx = 1;
                for (int i = ist; i <= ifn; i++) {
                    for (int j = jst; j <= jfn; j++) {
                        xmdkr1(j, i, l + 1) = (float)tmp(idx);     idx++;
                        xmdkr1(j, i, l + 3) = (float)tmp(idx);     idx++;
                        xmdkr1(j, i, l + 2) = (float)tmp(idx);     idx++;
                    }
                }
                for (int i = ist; i <= ifn; i++) {
                    for (int j = jst; j <= jfn; j++) {
                        xmdk(j, i, l + 1, n, iaes) =  (double)xmdkr1(j, i, l + 1);
                        xmdk(j, i, l + 2, n, iaes) = -(double)xmdkr1(j, i, l + 2);
                        xmdk(j, i, l + 3, n, iaes) =  (double)xmdkr1(j, i, l + 3);
                    }
                }
            }

        } // ns loop
    } // m loop

    // -----------------------------------------------------------------------
    // Read in modal deflections for i=const surfaces
    // -----------------------------------------------------------------------
    for (int m = 1; m <= 2; m++) {

        int i_val, l, nseg;
        if (m == 1) {
            i_val = 1;
            l     = 0;
            nseg  = nbci0(nbl);
        } else {
            i_val = idim;
            l     = 3;
            nseg  = nbcidim(nbl);
        }

        for (int ns = 1; ns <= nseg; ns++) {

            int lbc = std::abs(ibcinfo(nbl, ns, 1, m));

            int iflag = (lbc - 1006) * (lbc - 1005) * (lbc - 2004) * (lbc - 2014) *
                        (lbc - 2024) * (lbc - 2034) * (lbc - 2016);

            if (iflag != 0) continue;

            int jst = ibcinfo(nbl, ns, 2, m);
            int jfn = ibcinfo(nbl, ns, 3, m);
            int kst = ibcinfo(nbl, ns, 4, m);
            int kfn = ibcinfo(nbl, ns, 5, m);

            if (ialph == 0) {
                // read(iunit,*)((xmdi(j,k,l+1,n,iaes),xmdi(j,k,l+2,n,iaes),
                //               xmdi(j,k,l+3,n,iaes),j=jst,jfn),k=kst,kfn)
                int nvals = 3 * (jfn - jst + 1) * (kfn - kst + 1);
                FortranArray1D<double> tmp(nvals);
                fortran_read_list_array(iunit, tmp, 1, nvals);
                int idx = 1;
                for (int k = kst; k <= kfn; k++) {
                    for (int j = jst; j <= jfn; j++) {
                        xmdi(j, k, l + 1, n, iaes) = tmp(idx);     idx++;
                        xmdi(j, k, l + 2, n, iaes) = tmp(idx);     idx++;
                        xmdi(j, k, l + 3, n, iaes) = tmp(idx);     idx++;
                    }
                }
            } else {
                // read(iunit,*)((xmdir1(j,k,l+1),xmdir1(j,k,l+3),
                //               xmdir1(j,k,l+2),j=jst,jfn),k=kst,kfn)
                int nvals = 3 * (jfn - jst + 1) * (kfn - kst + 1);
                FortranArray1D<double> tmp(nvals);
                fortran_read_list_array(iunit, tmp, 1, nvals);
                int idx = 1;
                for (int k = kst; k <= kfn; k++) {
                    for (int j = jst; j <= jfn; j++) {
                        xmdir1(j, k, l + 1) = (float)tmp(idx);     idx++;
                        xmdir1(j, k, l + 3) = (float)tmp(idx);     idx++;
                        xmdir1(j, k, l + 2) = (float)tmp(idx);     idx++;
                    }
                }
                for (int k = kst; k <= kfn; k++) {
                    for (int j = jst; j <= jfn; j++) {
                        xmdi(j, k, l + 1, n, iaes) =  (double)xmdir1(j, k, l + 1);
                        xmdi(j, k, l + 2, n, iaes) = -(double)xmdir1(j, k, l + 2);
                        xmdi(j, k, l + 3, n, iaes) =  (double)xmdir1(j, k, l + 3);
                    }
                }
            }

        } // ns loop
    } // m loop

    // Deallocate real arrays in which modal data was read
    xmdjr1.deallocate();
    xmdkr1.deallocate();
    xmdir1.deallocate();

    return;
}

} // namespace modread_ns
