// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bcchk.h"
#include "bc.h"
#include "termn8.h"
#include <algorithm>
#include <cstdio>
#include <cstring>

namespace bcchk_ns {

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
        double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat,
        int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim,
              maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg,
              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
              lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bcchk(int& idim, int& jdim, int& kdim,
           FortranArray4DRef<double> q, FortranArray4DRef<double> qi0,
           FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
           FortranArray3DRef<double> blank, int& ibcflg, int& nbl,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, int& myid, int& istop,
           FortranArray1DRef<int> igridg, int& maxbl)
{
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    if (ibcflg == 0) {
        //
        // don't reset boundary values from initial values for checking if
        // adjacent interior point is blanked out...otherwise this routine (below)
        // will think there is a problem when in fact there is not
        //
        for (int i = 1; i <= idim1; i++) {
            for (int k = 1; k <= kdim; k++) {
                qj0(k,i,1,1) = (1.0 - blank(1,k,i)) * qj0(k,i,1,1)
                              - 1.0 * blank(1,k,i);
                qj0(k,i,5,1) = (1.0 - blank(1,k,i)) * qj0(k,i,5,1)
                              - 1.0 * blank(1,k,i);
                qj0(k,i,1,3) = (1.0 - blank(jdim1,k,i)) * qj0(k,i,1,3)
                              - 1.0 * blank(jdim1,k,i);
                qj0(k,i,5,3) = (1.0 - blank(jdim1,k,i)) * qj0(k,i,5,3)
                              - 1.0 * blank(jdim1,k,i);
            }
            for (int j = 1; j <= jdim; j++) {
                qk0(j,i,1,1) = (1.0 - blank(j,1,i)) * qk0(j,i,1,1)
                              - 1.0 * blank(j,1,i);
                qk0(j,i,5,1) = (1.0 - blank(j,1,i)) * qk0(j,i,5,1)
                              - 1.0 * blank(j,1,i);
                qk0(j,i,1,3) = (1.0 - blank(j,kdim1,i)) * qk0(j,i,1,3)
                              - 1.0 * blank(j,kdim1,i);
                qk0(j,i,5,3) = (1.0 - blank(j,kdim1,i)) * qk0(j,i,5,3)
                              - 1.0 * blank(j,kdim1,i);
            }
        }

        for (int k = 1; k <= kdim; k++) {
            for (int j = 1; j <= jdim; j++) {
                qi0(j,k,1,1) = (1.0 - blank(j,k,1)) * qi0(j,k,1,1)
                              - 1.0 * blank(j,k,1);
                qi0(j,k,5,1) = (1.0 - blank(j,k,1)) * qi0(j,k,5,1)
                              - 1.0 * blank(j,k,1);
                qi0(j,k,1,3) = (1.0 - blank(j,k,idim1)) * qi0(j,k,1,3)
                              - 1.0 * blank(j,k,idim1);
                qi0(j,k,5,3) = (1.0 - blank(j,k,idim1)) * qi0(j,k,5,3)
                              - 1.0 * blank(j,k,idim1);
            }
        }

    } else if (ibcflg == 1) {
        //
        // fill endpoints for safety with multi-plane vectorization
        //
        for (int l = 1; l <= 5; l++) {
            for (int i = 1; i <= idim1; i++) {
                q(jdim,kdim,i,l) = q(jdim1,kdim1,i,l);

                qj0(kdim,i,l,1) = q(1,kdim1,i,l);
                qj0(kdim,i,l,2) = q(1,kdim1,i,l);
                qj0(kdim,i,l,3) = q(jdim1,kdim1,i,l);
                qj0(kdim,i,l,4) = q(jdim1,kdim1,i,l);

                qk0(jdim,i,l,1) = qk0(jdim1,i,l,1);
                qk0(jdim,i,l,2) = qk0(jdim1,i,l,2);
                qk0(jdim,i,l,3) = qk0(jdim1,i,l,3);
                qk0(jdim,i,l,4) = qk0(jdim1,i,l,4);
            }
            for (int izz = 1; izz <= jdim1; izz++) {
                q(izz,kdim,idim,l) = q(izz,kdim1,idim1,l);

                qi0(izz,kdim,l,1) = q(izz,kdim1,1,l);
                qi0(izz,kdim,l,2) = q(izz,kdim1,1,l);
                qi0(izz,kdim,l,3) = q(izz,kdim1,idim1,l);
                qi0(izz,kdim,l,4) = q(izz,kdim1,idim1,l);
            }

            q(jdim,kdim,idim,l) = q(jdim1,kdim1,idim1,l);

            qi0(jdim,kdim,l,1) = q(jdim1,kdim1,1,l);
            qi0(jdim,kdim,l,2) = q(jdim1,kdim1,1,l);
            qi0(jdim,kdim,l,3) = q(jdim1,kdim1,idim1,l);
            qi0(jdim,kdim,l,4) = q(jdim1,kdim1,idim1,l);

            for (int k = 1; k <= kdim1; k++) {
                qi0(jdim,k,l,1) = qi0(jdim1,k,l,1);
                qi0(jdim,k,l,2) = qi0(jdim1,k,l,2);
                qi0(jdim,k,l,3) = qi0(jdim1,k,l,3);
                qi0(jdim,k,l,4) = qi0(jdim1,k,l,4);

                q(jdim,k,idim,l) = q(jdim1,k,idim,l);
            }
        }

    } else if (ibcflg == 2) {
        //
        // check for negative (or large) densities and pressures
        //
        double epsz = 1.0e-10;
        double epss = 1.0e+05;
        int m1 = -1;

        for (int i = 1; i <= idim1; i++) {
            for (int k = 1; k <= kdim; k++) {

                int ib = (int)blank(1,k,i);

                if (((float)qj0(k,i,1,1) <= (float)epsz ||
                     (float)qj0(k,i,5,1) <= (float)epsz ||
                     (float)qj0(k,i,1,1) >= (float)epss ||
                     (float)qj0(k,i,5,1) >= (float)epss) &&
                    ib > 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "boundary conditions resulted in negative (or large)"
                        " density and/or pressure");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " block%4d (grid %4d)  on j=1 boundary at k,i=%5d%5d"
                        " qj0(1),qj0(5)=%12.5e%12.5e",
                        nbl, igridg(nbl), k, i,
                        (float)qj0(k,i,1,1), (float)qj0(k,i,5,1));
                    if (qj0(k,i,1,1) == -1.0 && qj0(k,i,5,1) == -1.0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "After analysis, most probable cause is:"
                            " boundary data not set.");
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    } else {
                        if (istop > 0) {
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }

                ib = (int)blank(jdim1,k,i);

                if (((float)qj0(k,i,1,3) <= (float)epsz ||
                     (float)qj0(k,i,5,3) <= (float)epsz ||
                     (float)qj0(k,i,1,3) >= (float)epss ||
                     (float)qj0(k,i,5,3) >= (float)epss) &&
                    ib > 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "boundary conditions resulted in negative (or large)"
                        " density and/or pressure");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " block%4d (grid %4d)  on j=jdim boundary at"
                        " k,i=%5d%5d qj0(1),qj0(5)=%12.5e%12.5e",
                        nbl, igridg(nbl), k, i,
                        (float)qj0(k,i,1,3), (float)qj0(k,i,5,3));
                    if (qj0(k,i,1,3) == -1.0 && qj0(k,i,5,3) == -1.0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "After analysis, most probable cause is:"
                            " boundary data not set.");
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    } else {
                        if (istop > 0) {
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }
            }

            for (int j = 1; j <= jdim; j++) {

                int ib = (int)blank(j,1,i);

                if (((float)qk0(j,i,1,1) <= (float)epsz ||
                     (float)qk0(j,i,5,1) <= (float)epsz ||
                     (float)qk0(j,i,1,1) >= (float)epss ||
                     (float)qk0(j,i,5,1) >= (float)epss) &&
                    ib > 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "boundary conditions resulted in negative (or large)"
                        " density and/or pressure");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " block%4d (grid %4d)  on k=1 boundary at j,i=%5d%5d"
                        " qk0(1),qk0(5)=%12.5e%12.5e",
                        nbl, igridg(nbl), j, i,
                        (float)qk0(j,i,1,1), (float)qk0(j,i,5,1));
                    if (qk0(j,i,1,1) == -1.0 && qk0(j,i,5,1) == -1.0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "After analysis, most probable cause is:"
                            " boundary data not set.");
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    } else {
                        if (istop > 0) {
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }

                ib = (int)blank(j,kdim1,i);

                if (((float)qk0(j,i,1,3) <= (float)epsz ||
                     (float)qk0(j,i,5,3) <= (float)epsz ||
                     (float)qk0(j,i,1,3) >= (float)epss ||
                     (float)qk0(j,i,5,3) >= (float)epss) &&
                    ib > 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "boundary conditions resulted in negative (or large)"
                        " density and/or pressure");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " block%4d (grid %4d)  on k=kdim boundary at"
                        " j,i=%5d%5d qk0(1),qk0(5)=%12.5e%12.5e",
                        nbl, igridg(nbl), j, i,
                        (float)qk0(j,i,1,3), (float)qk0(j,i,5,3));
                    if ((float)qk0(j,i,1,3) == -1.0f &&
                        (float)qk0(j,i,5,3) == -1.0f) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "After analysis, most probable cause is:"
                            " boundary data not set.");
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    } else {
                        if (istop > 0) {
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }
            }
        }

        for (int k = 1; k <= kdim; k++) {
            for (int j = 1; j <= jdim; j++) {

                int ib = (int)blank(j,k,1);

                if (((float)qi0(j,k,1,1) <= (float)epsz ||
                     (float)qi0(j,k,5,1) <= (float)epsz ||
                     (float)qi0(j,k,1,1) >= (float)epss ||
                     (float)qi0(j,k,5,1) >= (float)epss) &&
                    ib > 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "boundary conditions resulted in negative (or large)"
                        " density and/or pressure");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " block%4d (grid %4d)  on i=1 boundary at j,k=%5d%5d"
                        " qi0(1),qi0(5)=%12.5e%12.5e",
                        nbl, igridg(nbl), j, k,
                        (float)qi0(j,k,1,1), (float)qi0(j,k,5,1));
                    if ((float)qi0(j,k,1,1) == -1.0f &&
                        (float)qi0(j,k,5,1) == -1.0f) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "After analysis, most probable cause is:"
                            " boundary data not set.");
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    } else {
                        if (istop > 0) {
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }

                ib = (int)blank(j,k,idim1);

                if (((float)qi0(j,k,1,3) <= (float)epsz ||
                     (float)qi0(j,k,5,3) <= (float)epsz ||
                     (float)qi0(j,k,1,3) >= (float)epss ||
                     (float)qi0(j,k,5,3) >= (float)epss) &&
                    ib > 0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        "boundary conditions resulted in negative (or large)"
                        " density and/or pressure");
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120,
                        " block%4d (grid %4d)  on i=idim boundary at"
                        " j,k=%5d%5d qi0(1),qi0(5)=%12.5e%12.5e",
                        nbl, igridg(nbl), j, k,
                        (float)qi0(j,k,1,3), (float)qi0(j,k,5,3));
                    if (qi0(j,k,1,3) == -1.0 && qi0(j,k,5,3) == -1.0) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            "After analysis, most probable cause is:"
                            " boundary data not set.");
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    } else {
                        if (istop > 0) {
                            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                        }
                    }
                }
            }
        }
    }
}

} // namespace bcchk_ns
