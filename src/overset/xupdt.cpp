// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "xupdt.h"
#include <cstdio>
#include <algorithm>

namespace xupdt_ns {

void xupdt(FortranArray4DRef<double> q, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, FortranArray3DRef<double> bcj0, FortranArray3DRef<double> bck0, FortranArray3DRef<double> bci0, int& maxbl, int& iitot, FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray3DRef<double> qb, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& int_updt)
{
    // temporarily set isklton = 0 to suppress output messages...
    // the messages are now output from bc_info
    int isklt_sav = cmn_sklton.isklton;
    cmn_sklton.isklton = 0;

    // Set index for loading into appropriate qb array
    //   (1 for q, 2 for vk0,vj0,vi0, 3 for tj0,tk0,ti0)
    int iset;
    if (ldim == 5) {
        iset = 1;
    } else if (ldim == 1) {
        iset = 2;
    } else {
        iset = 3;
    }

    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    int lsta = lbg(nbl);
    int lend = lsta - 1;

    // Block 1: update interior (fringe) q values
    if (ibpntsg(nbl, 1) > 0 && int_updt != 0) {
        lend = lsta + ibpntsg(nbl, 1) - 1;
        for (int l = lsta; l <= lend; l++) {
            for (int ll = 1; ll <= ldim; ll++) {
                q(jjbg(l), kkbg(l), iibg(l), ll) = qb(ibcg(l), ll, iset);
            }
        }
        if (cmn_sklton.isklton == 1) {
            // write(15,*)' interior values updated = ',ibpntsg(nbl,1)
            if (iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   hole    chimera grid interpolation     type    0");
            }
        }
    }

    // Block 2: update qj0 and bcj0
    if (ibpntsg(nbl, 2) > 0) {
        int mp0   = 0;
        int mpdim = 0;
        lsta = lend + 1;
        lend = lsta + ibpntsg(nbl, 2) - 1;
        for (int l = lsta; l <= lend; l++) {
            int mp  = std::max(0, 1 + jjbg(l) - jdim1) + 1 + std::max(0, -jjbg(l));
            int mpp = std::max(0, 3 - mp);
            mp0   = mp0   + std::min(1, mpp);
            mpp   = std::max(0, mp - 2);
            mpdim = mpdim + std::min(1, mpp);
            int mm = mp / 3 + 1;
            for (int ll = 1; ll <= ldim; ll++) {
                qj0(kkbg(l), iibg(l), ll, mp) = qb(ibcg(l), ll, iset);
                bcj0(kkbg(l), iibg(l), mm) = 0.0;
            }
        }
        if (cmn_sklton.isklton == 1) {
            // write(15,*)' QJ0 values updated = ',ibpntsg(nbl,2)
            if (mp0 > 0 && iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   j=   1  chimera grid interpolation     type    0");
            }
            if (mpdim > 0 && iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   j=jdim  chimera grid interpolation     type    0");
            }
        }
    }

    // Block 3: update qk0 and bck0
    if (ibpntsg(nbl, 3) > 0) {
        int mp0   = 0;
        int mpdim = 0;
        lsta = lend + 1;
        lend = lsta + ibpntsg(nbl, 3) - 1;
        for (int l = lsta; l <= lend; l++) {
            int mp  = std::max(0, 1 + kkbg(l) - kdim1) + 1 + std::max(0, -kkbg(l));
            int mpp = std::max(0, 3 - mp);
            mp0   = mp0   + std::min(1, mpp);
            mpp   = std::max(0, mp - 2);
            mpdim = mpdim + std::min(1, mpp);
            int mm = mp / 3 + 1;
            for (int ll = 1; ll <= ldim; ll++) {
                qk0(jjbg(l), iibg(l), ll, mp) = qb(ibcg(l), ll, iset);
                bck0(jjbg(l), iibg(l), mm) = 0.0;
            }
        }
        if (cmn_sklton.isklton == 1) {
            // write(15,*)' QK0 values updated = ',ibpntsg(nbl,3)
            if (mp0 > 0 && iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   k=   1  chimera grid interpolation     type    0");
            }
            if (mpdim > 0 && iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   k=kdim  chimera grid interpolation     type    0");
            }
        }
    }

    // Block 4: update qi0 and bci0
    if (ibpntsg(nbl, 4) > 0) {
        int mp0   = 0;
        int mpdim = 0;
        lsta = lend + 1;
        lend = lsta + ibpntsg(nbl, 4) - 1;
        for (int l = lsta; l <= lend; l++) {
            int mp  = std::max(0, 1 + iibg(l) - idim1) + 1 + std::max(0, -iibg(l));
            int mpp = std::max(0, 3 - mp);
            mp0   = mp0   + std::min(1, mpp);
            mpp   = std::max(0, mp - 2);
            mpdim = mpdim + std::min(1, mpp);
            int mm = mp / 3 + 1;
            for (int ll = 1; ll <= ldim; ll++) {
                qi0(jjbg(l), kkbg(l), ll, mp) = qb(ibcg(l), ll, iset);
                bci0(jjbg(l), kkbg(l), mm) = 0.0;
            }
        }
        if (cmn_sklton.isklton == 1) {
            // write(15,*)' QI0 values updated = ',ibpntsg(nbl,4)
            if (mp0 > 0 && iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   i=   1  chimera grid interpolation     type    0");
            }
            if (mpdim > 0 && iset == 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "   i=idim  chimera grid interpolation     type    0");
            }
        }
    }

    cmn_sklton.isklton = isklt_sav;
}

} // namespace xupdt_ns
