// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "dird.h"
#include "q8smin.h"
#include "ccomplex.h"
#include "runtime/fortran_array.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace dird_ns {

void dird(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> snj0, FortranArray3DRef<double> snk0, FortranArray3DRef<double> sni0, FortranArray3DRef<double> snjm, FortranArray3DRef<double> snkm, FortranArray3DRef<double> snim, FortranArray4DRef<double> w, FortranArray1DRef<int> ivisc, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    using namespace q8smin_ns;
    using namespace ccomplex_ns;

    int jdim1, kdim1, idim1;
    int i, j, k, km, jm, im;
    int icount;
    double snk0min, snk0mg, snj0min, snj0mg, sni0min, sni0mg;

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;


    // cell-center locations
    for (i = 1; i <= idim1; i++) {
    for (k = 1; k <= kdim1; k++) {
    for (j = 1; j <= jdim1; j++) {
        w(j,k,i,1) = (  x(j  ,k  ,i) + x(j  ,k  ,i+1) +
                        x(j  ,k+1,i) + x(j  ,k+1,i+1) +
                        x(j+1,k  ,i) + x(j+1,k  ,i+1) +
                        x(j+1,k+1,i) + x(j+1,k+1,i+1) );
        w(j,k,i,2) = (  y(j  ,k  ,i) + y(j  ,k  ,i+1) +
                        y(j  ,k+1,i) + y(j  ,k+1,i+1) +
                        y(j+1,k  ,i) + y(j+1,k  ,i+1) +
                        y(j+1,k+1,i) + y(j+1,k+1,i+1) );
        w(j,k,i,3) = (  z(j  ,k  ,i) + z(j  ,k  ,i+1) +
                        z(j  ,k+1,i) + z(j  ,k+1,i+1) +
                        z(j+1,k  ,i) + z(j+1,k  ,i+1) +
                        z(j+1,k+1,i) + z(j+1,k+1,i+1) );
    }}}

    // k=0 wall
    if (ivisc(3) > 1) {

        k = 1;
        for (i = 1; i <= idim1; i++) {
        for (j = 1; j <= jdim1; j++) {
            w(j,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) +
                                   x(j+1,k,i) + x(j+1,k,i+1) ) *sk(j,k,i,1)
                              + (  y(j  ,k,i) + y(j  ,k,i+1) +
                                   y(j+1,k,i) + y(j+1,k,i+1) ) *sk(j,k,i,2)
                              + (  z(j  ,k,i) + z(j  ,k,i+1) +
                                   z(j+1,k,i) + z(j+1,k,i+1) ) *sk(j,k,i,3) );
        }}

        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            snk0(j,k,i) = ( w(j,k,i,1)*sk(j,1,i,1)
                          + w(j,k,i,2)*sk(j,1,i,2)
                          + w(j,k,i,3)*sk(j,1,i,3) )*.125e0 - w(j,1,i,4);
        }}}

        icount  = 0;
        snk0mg  = 1.0e+06;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            // q8smin expects FortranArray1DRef<double>; snk0(1,k,i)..snk0(jdim1,k,i) are contiguous
            FortranArray1DRef<double> snk0_col(&snk0(1,k,i), jdim1);
            snk0min = q8smin(jdim1, snk0_col);

            if ((double)snk0min <= 0.) {
                // write(15,8001) - commented out in original
                if (icount <= 100) {
                    for (j = 1; j <= jdim1; j++) {
                        if ((double)snk0(j,k,i) > 0.e0) continue;
                        icount = icount + 1;
                        // if (icount.le.100) write(15,8903) j,k,i,real(snk0(j,k,i))
                    }
                }
            }

            snk0mg = ccmin(snk0min, snk0mg);
        }}

        // ensure grid wrapping back in toward body doesn't mess up snk0
        if (icount > 0) {
            // write(15,8111) - commented out in original
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   negative directed distances have been detected...");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   this is not usually a problem, but check output on");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   unit 15 to be sure (may need to check grid too)");
            for (i = 1; i <= idim1; i++) {
            for (j = 1; j <= jdim1; j++) {
            for (k = 2; k <= kdim1; k++) {
                snk0(j,k,i) = ccmax(snk0(j,k,i), snk0(j,k-1,i));
            }}}
        }
        // write(15,9901) real(snk0mg) - commented out in original

    } // end if ivisc(3)>1 for k=0 wall


    // k=kdim wall
    if (ivisc(3) > 1) {

        k  = kdim;
        km = k - 1;
        for (i = 1; i <= idim1; i++) {
        for (j = 1; j <= jdim1; j++) {
            w(j,km,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) +
                                    x(j+1,k,i) + x(j+1,k,i+1) ) *sk(j,k,i,1)
                               + (  y(j  ,k,i) + y(j  ,k,i+1) +
                                    y(j+1,k,i) + y(j+1,k,i+1) ) *sk(j,k,i,2)
                               + (  z(j  ,k,i) + z(j  ,k,i+1) +
                                    z(j+1,k,i) + z(j+1,k,i+1) ) *sk(j,k,i,3) );
        }}

        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            snkm(j,k,i) =  w(j,kdim-1,i,4) - ( w(j,k,i,1)*sk(j,kdim,i,1)
                          +                     w(j,k,i,2)*sk(j,kdim,i,2)
                          +                     w(j,k,i,3)*sk(j,kdim,i,3) )
                                             * .125e0;
        }}}

        // write(15,11) - commented out in original
        // write(15,1601) loops - commented out in original

        icount  = 0;
        snk0mg  = 1.0e+06;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            FortranArray1DRef<double> snkm_col(&snkm(1,k,i), jdim1);
            snk0min = q8smin(jdim1, snkm_col);

            if ((double)snk0min <= 0.) {
                // write(15,801) - commented out in original
                if (icount <= 100) {
                    for (j = 1; j <= jdim1; j++) {
                        if ((double)snkm(j,k,i) > 0.e0) continue;
                        icount = icount + 1;
                        // if (icount.le.100) write(15,893) j,k,i,real(snkm(j,k,i))
                    }
                }
            }

            snk0mg = ccmin(snk0min, snk0mg);
        }}

        // ensure grid wrapping back in toward body doesn't mess up snkm
        if (icount > 0) {
            // write(15,811) - commented out in original
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   negative directed distances have been detected...");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   this is not usually a problem, but check output on");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   unit 15 to be sure (may need to check grid too)");
            for (i = 1; i <= idim1; i++) {
            for (j = 1; j <= jdim1; j++) {
            for (k = 1; k <= kdim1-1; k++) {
                snkm(j,k,i) = ccmax(snkm(j,k,i), snkm(j,k+1,i));
            }}}
        }
        // write(15,1991) real(snk0mg) - commented out in original

    } // end if ivisc(3)>1 for k=kdim wall


    // j=0 wall
    if (ivisc(2) > 1) {

        j = 1;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            w(j,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) +
                                   x(j,k+1,i) + x(j,k+1,i+1) ) *sj(j,k,i,1)
                              + (  y(j  ,k,i) + y(j  ,k,i+1) +
                                   y(j,k+1,i) + y(j,k+1,i+1) ) *sj(j,k,i,2)
                              + (  z(j  ,k,i) + z(j  ,k,i+1) +
                                   z(j,k+1,i) + z(j,k+1,i+1) ) *sj(j,k,i,3) );
        }}

        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            snj0(j,k,i) = ( w(j,k,i,1)*sj(1,k,i,1)
                          + w(j,k,i,2)*sj(1,k,i,2)
                          + w(j,k,i,3)*sj(1,k,i,3) )*.125e0 - w(1,k,i,4);
        }}}

        icount  = 0;
        snj0mg  = 1.0e+06;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            FortranArray1DRef<double> snj0_col(&snj0(1,k,i), jdim1);
            snj0min = q8smin(jdim1, snj0_col);

            if ((double)snj0min <= 0.) {
                // write(15,8001) - commented out in original
                if (icount <= 100) {
                    for (j = 1; j <= jdim1; j++) {
                        if ((double)snj0(j,k,i) > 0.e0) continue;
                        icount = icount + 1;
                        // if (icount.le.100) write(15,8703) j,k,i,real(snj0(j,k,i))
                    }
                }
            }

            snj0mg = ccmin(snj0min, snj0mg);
        }}

        // ensure grid wrapping back in toward body doesn't mess up snj0
        if (icount > 0) {
            // write(15,8111) - commented out in original
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   negative directed distances have been detected...");
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {
            for (j = 2; j <= jdim1; j++) {
                snj0(j,k,i) = ccmax(snj0(j,k,i), snj0(j-1,k,i));
            }}}
        }
        // write(15,9701) real(snj0mg) - commented out in original

    } // end if ivisc(2)>1 for j=0 wall


    // j=jdim wall
    if (ivisc(2) > 1) {

        j  = jdim;
        jm = j - 1;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            w(jm,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k,i+1) +
                                    x(j,k+1,i) + x(j,k+1,i+1) ) *sj(j,k,i,1)
                               + (  y(j  ,k,i) + y(j  ,k,i+1) +
                                    y(j,k+1,i) + y(j,k+1,i+1) ) *sj(j,k,i,2)
                               + (  z(j  ,k,i) + z(j  ,k,i+1) +
                                    z(j,k+1,i) + z(j,k+1,i+1) ) *sj(j,k,i,3) );
        }}

        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            snjm(j,k,i) =  w(jdim-1,k,i,4) - ( w(j,k,i,1)*sj(jdim,k,i,1)
                          +                     w(j,k,i,2)*sj(jdim,k,i,2)
                          +                     w(j,k,i,3)*sj(jdim,k,i,3) )
                                             * .125e0;
        }}}

        // write(15,116) - commented out in original
        // write(15,701) loops - commented out in original

        icount  = 0;
        snj0mg  = 1.0e+06;
        {
            FortranArray1D<double> snjm_tmp(kdim1);
            for (i = 1; i <= idim1; i++) {
            for (j = 1; j <= jdim1; j++) {
                // snjm(j,1,i)..snjm(j,kdim1,i) are not contiguous (j varies fastest)
                // copy to temp array for q8smin
                for (k = 1; k <= kdim1; k++) {
                    snjm_tmp(k) = snjm(j,k,i);
                }
                FortranArray1DRef<double> snjm_ref = snjm_tmp.ref();
                snj0min = q8smin(kdim1, snjm_ref);

                if ((double)snj0min <= 0.) {
                    // write(15,1811) - commented out in original
                    if (icount <= 100) {
                        for (k = 1; k <= kdim1; k++) {
                            if ((double)snjm(j,k,i) > 0.e0) continue;
                            icount = icount + 1;
                            // if (icount.le.100) write(15,895) j,k,i,real(snjm(j,k,i))
                        }
                    }
                }

                snj0mg = ccmin(snj0min, snj0mg);
            }}
        }

        // ensure grid wrapping back in toward body doesn't mess up snjm
        if (icount > 0) {
            // write(15,6811) - commented out in original
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   negative directed distances have been detected...");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   this is not usually a problem, but check output on");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   unit 15 to be sure (may need to check grid too)");
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {
            for (j = 1; j <= jdim1-1; j++) {
                snjm(j,k,i) = ccmax(snjm(j,k,i), snjm(j+1,k,i));
            }}}
        }
        // write(15,6991) real(snj0mg) - commented out in original

    } // end if ivisc(2)>1 for j=jdim wall


    // i=0 wall
    if (ivisc(1) > 1) {

        i = 1;
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            w(j,k,i,4) = .25e0*((  x(j  ,k,i) + x(j  ,k+1,i) +
                                   x(j+1,k,i) + x(j+1,k+1,i) ) *si(j,k,i,1)
                              + (  y(j  ,k,i) + y(j  ,k+1,i) +
                                   y(j+1,k,i) + y(j+1,k+1,i) ) *si(j,k,i,2)
                              + (  z(j  ,k,i) + z(j  ,k+1,i) +
                                   z(j+1,k,i) + z(j+1,k+1,i) ) *si(j,k,i,3) );
        }}

        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            sni0(j,k,i) = ( w(j,k,i,1)*si(j,k,1,1)
                          + w(j,k,i,2)*si(j,k,1,2)
                          + w(j,k,i,3)*si(j,k,1,3) )*.125e0 - w(j,k,1,4);
        }}}

        icount  = 0;
        sni0mg  = 1.0e+06;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            FortranArray1DRef<double> sni0_col(&sni0(1,k,i), jdim1);
            sni0min = q8smin(jdim1, sni0_col);

            if ((double)sni0min <= 0.) {
                // write(15,8001) - commented out in original
                if (icount <= 100) {
                    for (j = 1; j <= jdim1; j++) {
                        if ((double)sni0(j,k,i) > 0.e0) continue;
                        icount = icount + 1;
                        // if (icount.le.100) write(15,8003) j,k,i,real(sni0(j,k,i))
                    }
                }
            }

            sni0mg = ccmin(sni0min, sni0mg);
        }}

        // ensure grid wrapping back in toward body doesn't mess up sni0
        if (icount > 0) {
            // write(15,8111) - commented out in original
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   negative directed distances have been detected...");
            for (k = 1; k <= kdim1; k++) {
            for (j = 1; j <= jdim1; j++) {
            for (i = 2; i <= idim1; i++) {
                sni0(j,k,i) = ccmax(sni0(j,k,i), sni0(j,k,i-1));
            }}}
        }
        // write(15,9001) real(sni0mg) - commented out in original

    } // end if ivisc(1)>1 for i=0 wall


    // i=idim wall
    if (ivisc(1) > 1) {

        i  = idim;
        im = i - 1;
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            w(j,k,im,4) = .25e0*((  x(j  ,k,i) + x(j  ,k+1,i) +
                                    x(j+1,k,i) + x(j+1,k+1,i) ) *si(j,k,i,1)
                               + (  y(j  ,k,i) + y(j  ,k+1,i) +
                                    y(j+1,k,i) + y(j+1,k+1,i) ) *si(j,k,i,2)
                               + (  z(j  ,k,i) + z(j  ,k+1,i) +
                                    z(j+1,k,i) + z(j+1,k+1,i) ) *si(j,k,i,3) );
        }}

        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            snim(j,k,i) =  w(j,k,idim-1,4) - ( w(j,k,i,1)*si(j,k,idim,1)
                          +                     w(j,k,i,2)*si(j,k,idim,2)
                          +                     w(j,k,i,3)*si(j,k,idim,3) )
                                             * .125e0;
        }}}

        // write(15,126) - commented out in original
        // write(15,704) loops - commented out in original

        icount  = 0;
        sni0mg  = 1.0e+06;
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
            FortranArray1DRef<double> snim_col(&snim(1,k,i), jdim1);
            sni0min = q8smin(jdim1, snim_col);

            if ((double)sni0min <= 0.) {
                // write(15,817) - commented out in original
                if (icount <= 100) {
                    for (j = 1; j <= jdim1; j++) {
                        if ((double)snim(j,k,i) > 0.e0) continue;
                        icount = icount + 1;
                        // if (icount.le.100) write(15,899) j,k,i,real(snim(j,k,i))
                    }
                }
            }

            sni0mg = ccmin(sni0min, sni0mg);
        }}

        // ensure grid wrapping back in toward body doesn't mess up snim
        if (icount > 0) {
            // write(15,6711) - commented out in original
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   negative directed distances have been detected..."
                "\n   this is not usually a problem, but check output on"
                "\n   unit 15 to be sure (may need to check grid too)");
            for (k = 1; k <= kdim1; k++) {
            for (j = 1; j <= jdim1; j++) {
            for (i = 1; i <= idim1-1; i++) {
                snim(j,k,i) = ccmax(snim(j,k,i), snim(j,k,i+1));
            }}}
        }
        // write(15,6997) real(sni0mg) - commented out in original

    } // end if ivisc(1)>1 for i=idim wall

} // end dird

} // namespace dird_ns
