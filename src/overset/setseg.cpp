// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "setseg.h"
#include <cmath>

namespace setseg_ns {

void setseg(int& maxgr, int& maxbl, int& maxseg,
            FortranArray1DRef<int> nblg, FortranArray1DRef<int> ncgg,
            FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
            FortranArray1DRef<int> kdimg,
            FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
            FortranArray4DRef<int> kbcinfo,
            FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim,
            FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim,
            FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim,
            int& ndefrm0,
            FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nsegdfrm,
            FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
            int& maxsegdg, int& ngrid)
{
    int ng, nbl, ncg, idim, jdim, kdim, nseg, iseg, n;

    for (ng = 1; ng <= ngrid; ng++) {
        nbl  = nblg(ng);
        ncg  = ncgg(ng);
        idim = idimg(nbl);
        jdim = jdimg(nbl);
        kdim = kdimg(nbl);

        // i=0 face
        for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
            if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2004 ||
                std::abs(ibcinfo(nbl,nseg,1,1)) == 2014 ||
                std::abs(ibcinfo(nbl,nseg,1,1)) == 2024 ||
                std::abs(ibcinfo(nbl,nseg,1,1)) == 2034 ||
                std::abs(ibcinfo(nbl,nseg,1,1)) == 2016 ||
                std::abs(ibcinfo(nbl,nseg,1,1)) == 1005 ||
                std::abs(ibcinfo(nbl,nseg,1,1)) == 1006) {
                ndefrm0           = ndefrm0 + 1;
                idefrm(nbl)       = 1;
                iseg              = nsegdfrm(nbl) + 1;
                nsegdfrm(nbl)     = iseg;
                icsi(nbl,iseg)    = 1;
                icsf(nbl,iseg)    = 1;
                jcsi(nbl,iseg)    = ibcinfo(nbl,nseg,2,1);
                jcsf(nbl,iseg)    = ibcinfo(nbl,nseg,3,1);
                kcsi(nbl,iseg)    = ibcinfo(nbl,nseg,4,1);
                kcsf(nbl,iseg)    = ibcinfo(nbl,nseg,5,1);
            }
        }

        // i=idim face
        for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
            if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2004 ||
                std::abs(ibcinfo(nbl,nseg,1,2)) == 2014 ||
                std::abs(ibcinfo(nbl,nseg,1,2)) == 2024 ||
                std::abs(ibcinfo(nbl,nseg,1,2)) == 2034 ||
                std::abs(ibcinfo(nbl,nseg,1,2)) == 2016 ||
                std::abs(ibcinfo(nbl,nseg,1,2)) == 1005 ||
                std::abs(ibcinfo(nbl,nseg,1,2)) == 1006) {
                ndefrm0           = ndefrm0 + 1;
                idefrm(nbl)       = 1;
                iseg              = nsegdfrm(nbl) + 1;
                nsegdfrm(nbl)     = iseg;
                icsi(nbl,iseg)    = idimg(nbl);
                icsf(nbl,iseg)    = idimg(nbl);
                jcsi(nbl,iseg)    = ibcinfo(nbl,nseg,2,2);
                jcsf(nbl,iseg)    = ibcinfo(nbl,nseg,3,2);
                kcsi(nbl,iseg)    = ibcinfo(nbl,nseg,4,2);
                kcsf(nbl,iseg)    = ibcinfo(nbl,nseg,5,2);
            }
        }

        // j=0 face
        for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
            if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2004 ||
                std::abs(jbcinfo(nbl,nseg,1,1)) == 2014 ||
                std::abs(jbcinfo(nbl,nseg,1,1)) == 2024 ||
                std::abs(jbcinfo(nbl,nseg,1,1)) == 2034 ||
                std::abs(jbcinfo(nbl,nseg,1,1)) == 2016 ||
                std::abs(jbcinfo(nbl,nseg,1,1)) == 1005 ||
                std::abs(jbcinfo(nbl,nseg,1,1)) == 1006) {
                ndefrm0           = ndefrm0 + 1;
                idefrm(nbl)       = 1;
                iseg              = nsegdfrm(nbl) + 1;
                nsegdfrm(nbl)     = iseg;
                icsi(nbl,iseg)    = jbcinfo(nbl,nseg,2,1);
                icsf(nbl,iseg)    = jbcinfo(nbl,nseg,3,1);
                jcsi(nbl,iseg)    = 1;
                jcsf(nbl,iseg)    = 1;
                kcsi(nbl,iseg)    = jbcinfo(nbl,nseg,4,1);
                kcsf(nbl,iseg)    = jbcinfo(nbl,nseg,5,1);
            }
        }

        // j=jdim face
        for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
            if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2004 ||
                std::abs(jbcinfo(nbl,nseg,1,2)) == 2014 ||
                std::abs(jbcinfo(nbl,nseg,1,2)) == 2024 ||
                std::abs(jbcinfo(nbl,nseg,1,2)) == 2034 ||
                std::abs(jbcinfo(nbl,nseg,1,2)) == 2016 ||
                std::abs(jbcinfo(nbl,nseg,1,2)) == 1005 ||
                std::abs(jbcinfo(nbl,nseg,1,2)) == 1006) {
                ndefrm0           = ndefrm0 + 1;
                idefrm(nbl)       = 1;
                iseg              = nsegdfrm(nbl) + 1;
                nsegdfrm(nbl)     = iseg;
                icsi(nbl,iseg)    = jbcinfo(nbl,nseg,2,2);
                icsf(nbl,iseg)    = jbcinfo(nbl,nseg,3,2);
                jcsi(nbl,iseg)    = jdimg(nbl);
                jcsf(nbl,iseg)    = jdimg(nbl);
                kcsi(nbl,iseg)    = jbcinfo(nbl,nseg,4,2);
                kcsf(nbl,iseg)    = jbcinfo(nbl,nseg,5,2);
            }
        }

        // k=0 face
        for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
            if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2004 ||
                std::abs(kbcinfo(nbl,nseg,1,1)) == 2014 ||
                std::abs(kbcinfo(nbl,nseg,1,1)) == 2024 ||
                std::abs(kbcinfo(nbl,nseg,1,1)) == 2034 ||
                std::abs(kbcinfo(nbl,nseg,1,1)) == 2016 ||
                std::abs(kbcinfo(nbl,nseg,1,1)) == 1005 ||
                std::abs(kbcinfo(nbl,nseg,1,1)) == 1006) {
                ndefrm0           = ndefrm0 + 1;
                idefrm(nbl)       = 1;
                iseg              = nsegdfrm(nbl) + 1;
                nsegdfrm(nbl)     = iseg;
                icsi(nbl,iseg)    = kbcinfo(nbl,nseg,2,1);
                icsf(nbl,iseg)    = kbcinfo(nbl,nseg,3,1);
                jcsi(nbl,iseg)    = kbcinfo(nbl,nseg,4,1);
                jcsf(nbl,iseg)    = kbcinfo(nbl,nseg,5,1);
                kcsi(nbl,iseg)    = 1;
                kcsf(nbl,iseg)    = 1;
            }
        }

        // k=kdim face
        for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
            if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2004 ||
                std::abs(kbcinfo(nbl,nseg,1,2)) == 2014 ||
                std::abs(kbcinfo(nbl,nseg,1,2)) == 2024 ||
                std::abs(kbcinfo(nbl,nseg,1,2)) == 2034 ||
                std::abs(kbcinfo(nbl,nseg,1,2)) == 2016 ||
                std::abs(kbcinfo(nbl,nseg,1,2)) == 1005 ||
                std::abs(kbcinfo(nbl,nseg,1,2)) == 1006) {
                ndefrm0           = ndefrm0 + 1;
                idefrm(nbl)       = 1;
                iseg              = nsegdfrm(nbl) + 1;
                nsegdfrm(nbl)     = iseg;
                icsi(nbl,iseg)    = kbcinfo(nbl,nseg,2,2);
                icsf(nbl,iseg)    = kbcinfo(nbl,nseg,3,2);
                jcsi(nbl,iseg)    = kbcinfo(nbl,nseg,4,2);
                jcsf(nbl,iseg)    = kbcinfo(nbl,nseg,5,2);
                kcsi(nbl,iseg)    = kdimg(nbl);
                kcsf(nbl,iseg)    = kdimg(nbl);
            }
        }

        // Coarser grid levels
        if (ncg > 0) {
            for (n = 1; n <= ncg; n++) {
                nbl = nbl + 1;
                idefrm(nbl)   = idefrm(nbl-1);
                nsegdfrm(nbl) = nsegdfrm(nbl-1);
                for (iseg = 1; iseg <= nsegdfrm(nbl); iseg++) {
                    icsi(nbl,iseg) = icsi(nbl-1,iseg)/2 + 1;
                    icsf(nbl,iseg) = icsf(nbl-1,iseg)/2 + 1;
                    jcsi(nbl,iseg) = jcsi(nbl-1,iseg)/2 + 1;
                    jcsf(nbl,iseg) = jcsf(nbl-1,iseg)/2 + 1;
                    kcsi(nbl,iseg) = kcsi(nbl-1,iseg)/2 + 1;
                    kcsf(nbl,iseg) = kcsf(nbl-1,iseg)/2 + 1;
                }
            }
        }
    }

    return;
}

} // namespace setseg_ns
