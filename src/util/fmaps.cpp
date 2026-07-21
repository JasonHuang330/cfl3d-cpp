// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fmaps.h"
#include "ccomplex.h"
#include <cmath>

namespace fmaps_ns {

void fmaps(FortranArray1DRef<double> nx, FortranArray1DRef<double> ny, FortranArray1DRef<double> nz,
           FortranArray1DRef<double> area, FortranArray1DRef<double> at,
           FortranArray2DRef<double> f, FortranArray2DRef<double> qr, FortranArray2DRef<double> ql,
           int& n, int& nvtq)
{
    // local real scalars (translated as double to match array element type)
    double m_r, m_l, m_av, m_up, m_lo, m_max1, m_min1, m_dif, m_upd, m_max;
    double m_0, m_m1;

    // COMMON /fluid/
    double gamma = (double)cmn_fluid.gamma;
    double gm1   = (double)cmn_fluid.gm1;

    double x1   = gamma / gm1;
    double c1   = 1.0e0 / gm1;
    double eps4 = 0.0e0;

    for (int i = 1; i <= n; i++) {

        // ------------------------------------------------------------------
        // left/right density, velocity, pressure, enthalpy and sound speed
        // ------------------------------------------------------------------

        double r_l  = ql(i,1);
        double vx_l = ql(i,2);
        double vy_l = ql(i,3);
        double vz_l = ql(i,4);
        double p_l  = ql(i,5);
        double h_l  = x1*ql(i,5)/ql(i,1) + .5e0*(ql(i,2)*ql(i,2) + ql(i,3)*ql(i,3)
                                                  + ql(i,4)*ql(i,4));
        double c_l  = std::sqrt(gamma*ql(i,5)/ql(i,1));

        double r_r  = qr(i,1);
        double vx_r = qr(i,2);
        double vy_r = qr(i,3);
        double vz_r = qr(i,4);
        double p_r  = qr(i,5);
        double h_r  = x1*qr(i,5)/qr(i,1) + .5e0*(qr(i,2)*qr(i,2) + qr(i,3)*qr(i,3)
                                                  + qr(i,4)*qr(i,4));
        double c_r  = std::sqrt(gamma*qr(i,5)/qr(i,1));

        double c_av = 0.5*(c_l + c_r);

        // nx, ny, nz are the components of the normal vector at a cell face
        // anx, any, anz are the face vector components { (nx, ny, nz) * area }

        double anx = nx(i)*area(i);
        double any = ny(i)*area(i);
        double anz = nz(i)*area(i);

        double q_l   = (vx_l * nx(i) + vy_l * ny(i) + vz_l * nz(i));
        double q_r   = (vx_r * nx(i) + vy_r * ny(i) + vz_r * nz(i));

        m_l    = q_l / c_av;
        m_r    = q_r / c_av;

        // ------------------------------------------------------------------
        // compute conditions at cell face:
        // ------------------------------------------------------------------

        double am_l   = ccomplex_ns::ccabs(m_l);
        double am_r   = ccomplex_ns::ccabs(m_r);

        double one_d  = 1.0;
        double zero_d = 0.0;
        double sm_l   = ccomplex_ns::ccsignrc(one_d, m_l);
        double sm_r   = ccomplex_ns::ccsignrc(one_d, m_r);

        m_av   = 0.5*(m_l + m_r);
        m_up   = ccomplex_ns::ccmax(am_l, am_r);
        m_lo   = ccomplex_ns::ccmin(am_l, am_r);
        m_max1 = ccomplex_ns::ccmincr(m_up, one_d);
        m_min1 = ccomplex_ns::ccmincr(m_lo, one_d);
        m_dif  = 0.5*(am_r - am_l);

        double two_m_max1_m1 = 2.0*m_max1 - 1.0;
        double two_m_min1_m1 = 2.0*m_min1 - 1.0;
        double b_m    = ccomplex_ns::ccmaxrc(zero_d, two_m_max1_m1);
        double b_p    = ccomplex_ns::ccmaxrc(zero_d, two_m_min1_m1);

        // no lower bound on dissipation

        m_upd  = m_up;
        double dfak   = 1.0;

        // take aspect ratio into account
        //
        //     fasp   = scale**zeta
        //     dfak   = ccmaxcr(fasp, 1.0)
        //     dcut   = ccmincr(fasp, 1.0)
        //     dcut   = ccmax(hfl1*dcut, hfl2)
        //     m_upd  = ccmax(m_up, dcut):

        double q_cd   = 0.5 * area(i) * c_av * ( m_av - b_m*m_dif );
        double q_up   = 0.5 * area(i) * c_av *   m_upd * dfak;

        double p_av  = 0.5*(p_r + p_l);
        double p_ds  = 0.5*(p_r*sm_r - p_l*sm_l);
        double p_12  = p_av - b_p*p_ds;

        // ------------------------------------------------------------------
        // compute flux through cell face (maps scheme)
        // ------------------------------------------------------------------

        f(i,1)       = q_cd * (r_r      + r_l) -
                       q_up * (r_r      - r_l);
        f(i,2)       = q_cd * (r_r*vx_r + r_l*vx_l) -
                       q_up * (r_r*vx_r - r_l*vx_l) + anx * p_12;
        f(i,3)       = q_cd * (r_r*vy_r + r_l*vy_l) -
                       q_up * (r_r*vy_r - r_l*vy_l) + any * p_12;
        f(i,4)       = q_cd * (r_r*vz_r + r_l*vz_l) -
                       q_up * (r_r*vz_r - r_l*vz_l) + anz * p_12;
        f(i,5)       = q_cd * (r_r*h_r  + r_l*h_l)  -
                       q_up * (r_r*h_r  - r_l*h_l);

        // need to check for moving grids (gs related to at(i))
        //
        //     f(i,5)       = q_cd * (r_r*h_r  + r_l*h_l)  -
        //    .               q_up * (r_r*h_r  - r_l*h_l)  + gs * area(i) * p_12

        // ------------------------------------------------------------------
        // additional terms (maps+ scheme)
        // ------------------------------------------------------------------
        //
        //     rho   = ccmin(r_l, r_r)
        //     h     = ccmin(h_l, h_r)

        double p_d   = p_r - p_l;
        double q_d   = q_r - q_l;
        double rho   = 0.5*( r_l + r_r);
        double u     = 0.5*(vx_l + vx_r);
        double v     = 0.5*(vy_l + vy_r);
        double w     = 0.5*(vz_l + vz_r);
        double h     = 0.5*( h_l + h_r);
        double c_min = ccomplex_ns::ccmin(c_l, c_r);
        double c_max = ccomplex_ns::ccmax(c_l, c_r);
        double q_min = m_lo*c_min;
        double q_max = m_up*c_max;
        double qq    = u*u + v*v + w*w;
        double cc    = c_max*c_max;

        // ------------------------------------------------------------------
        // rescale speed of sound for preconditioning
        // ------------------------------------------------------------------
        //
        //     qs     = q_max
        //     qref2  = ccmin( ccmax( upc*qq, epslocm2*cc), cc)
        //     alp    = (1.-qref2/cc)/2.0
        //     q_maxp = q_max*(1.-alp)
        //     q_minp = q_maxp
        //     c_maxp = sqrt(qs*qs * alp*alp + qref2)
        //     c_minp = c_maxp
        //
        //     if (no preconditioning) then (var, ...) else (varp, ...) */
        //     (globals.precond.preconditioning == 0) ?
        //     (q_max = q_max , q_min = q_min , c_max = c_max , c_min = c_min ) :
        //     (q_max = q_maxp, q_min = q_minp, c_max = c_maxp, c_min = c_minp)

        double rc    = 1.0/c_max;
        m_max = q_max/c_min;
        m_0   = ccomplex_ns::ccmincr(m_max, one_d);
        //     m_m1  = 2.0*ccmaxcr(0.5-m_0, 0.0)
        double one_minus_m0 = 1.0 - m_0;
        m_m1  = ccomplex_ns::ccmaxcr(one_minus_m0, zero_d);

        double p_scal = 0.5 * rc  * m_m1         * p_d;
        double q_scal = 0.5 * rho * m_m1 * c_min * q_d;

        f(i,1) = f(i,1) - area(i) *     p_scal;
        f(i,2) = f(i,2) - area(i) * u * p_scal - anx * q_scal;
        f(i,3) = f(i,3) - area(i) * v * p_scal - any * q_scal;
        f(i,4) = f(i,4) - area(i) * w * p_scal - anz * q_scal;
        f(i,5) = f(i,5) - area(i) * h * p_scal;

    } // end do 1000

    return;
}

} // namespace fmaps_ns
