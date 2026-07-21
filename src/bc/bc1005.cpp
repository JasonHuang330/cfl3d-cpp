// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc1005.h"
#include "ccomplex.h"
#include <cmath>
#include <algorithm>

namespace bc1005_ns {

using namespace ccomplex_ns;

// bc dispatcher: loops over all blocks, dispatches to bc1005 for bctype 1005 or 1006
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
    int iuns = 0;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        FortranArray4DRef<double> q    (&w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0  (&w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0  (&w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0  (&w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj   (&w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk   (&w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si   (&w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj  (&w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck  (&w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci  (&w(lw(10,ibl)), jdim,   kdim,   2);
        FortranArray4DRef<double> xtbj (&w(lw(11,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbk (&w(lw(12,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbi (&w(lw(13,ibl)), jdim,   kdim,   3,      2);
        FortranArray4DRef<double> atbj (&w(lw(14,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbk (&w(lw(15,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbi (&w(lw(16,ibl)), jdim,   kdim,   3,      2);

        FortranArray4DRef<double> tursav(&w(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&w(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&w(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&w(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&w(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&w(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&w(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&w(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // j-face segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            int bctype = jbcinfo(ibl, iseg, 1, 1);
            if (bctype == 1005 || bctype == 1006) {
                int ista      = jbcinfo(ibl, iseg, 1, 2);
                int iend      = jbcinfo(ibl, iseg, 1, 3);
                int jsta      = jbcinfo(ibl, iseg, 1, 4);
                int jend      = jbcinfo(ibl, iseg, 1, 5);
                int ksta      = jbcinfo(ibl, iseg, 1, 6);
                int kend      = jbcinfo(ibl, iseg, 1, 7);
                int nface     = jbcinfo(ibl, iseg, 1, 8);
                int inormmom  = jbcinfo(ibl, iseg, 1, 9);
                bc1005(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, inormmom,
                       nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // k-face segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            int bctype = kbcinfo(ibl, iseg, 1, 1);
            if (bctype == 1005 || bctype == 1006) {
                int ista      = kbcinfo(ibl, iseg, 1, 2);
                int iend      = kbcinfo(ibl, iseg, 1, 3);
                int jsta      = kbcinfo(ibl, iseg, 1, 4);
                int jend      = kbcinfo(ibl, iseg, 1, 5);
                int ksta      = kbcinfo(ibl, iseg, 1, 6);
                int kend      = kbcinfo(ibl, iseg, 1, 7);
                int nface     = kbcinfo(ibl, iseg, 1, 8);
                int inormmom  = kbcinfo(ibl, iseg, 1, 9);
                bc1005(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, inormmom,
                       nou, bou, nbuf, ibufdim, nummem);
            }
        }

        // i-face segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            int bctype = ibcinfo(ibl, iseg, 1, 1);
            if (bctype == 1005 || bctype == 1006) {
                int ista      = ibcinfo(ibl, iseg, 1, 2);
                int iend      = ibcinfo(ibl, iseg, 1, 3);
                int jsta      = ibcinfo(ibl, iseg, 1, 4);
                int jend      = ibcinfo(ibl, iseg, 1, 5);
                int ksta      = ibcinfo(ibl, iseg, 1, 6);
                int kend      = ibcinfo(ibl, iseg, 1, 7);
                int nface     = ibcinfo(ibl, iseg, 1, 8);
                int inormmom  = ibcinfo(ibl, iseg, 1, 9);
                bc1005(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, inormmom,
                       nou, bou, nbuf, ibufdim, nummem);
            }
        }
    }
}

// bc1005: Set inviscid surface boundary conditions (bctype 1005/1006)
void bc1005(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj, FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& iuns, int& inormmom,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem)
{
    // COMMON block references
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;
    int*   ivisc   = cmn_reyue.ivisc;
    int&   iexp    = cmn_zero.iexp;

    // Note: (10.**(-iexp) is machine zero)
    double xminn = std::pow(10.0, (double)(-iexp+1));
    double one_d = 1.0;  // for ccmaxcr calls requiring lvalue reference

    int jdim1 = jdim-1;
    int kdim1 = kdim-1;
    int idim1 = idim-1;

    int jend1 = jend-1;
    int kend1 = kend-1;
    int iend1 = iend-1;

// ******************************************************************************
//      j=1 boundary             inviscid surface                    bctype 1005
// ******************************************************************************
    if (nface == 3) {

    for (int i = ista; i <= iend1; i++) {
    for (int k = ksta; k <= kend1; k++) {
        double pb            = q(1,k,i,5);
        qj0(k,i,1,1)        = q(1,k,i,1);

        double contra        = q(1,k,i,2)*sj(1,k,i,1)
                             + q(1,k,i,3)*sj(1,k,i,2)
                             + q(1,k,i,4)*sj(1,k,i,3) + sj(1,k,i,5);

        qj0(k,i,2,1) = q(1,k,i,2) - contra*sj(1,k,i,1);
        qj0(k,i,3,1) = q(1,k,i,3) - contra*sj(1,k,i,2);
        qj0(k,i,4,1) = q(1,k,i,4) - contra*sj(1,k,i,3);

        if (inormmom == 1) {
            double sixa = si(1,k,i,1)*si(1,k,i,4) + si(1,k,i+1,1)*si(1,k,i+1,4);
            double siya = si(1,k,i,2)*si(1,k,i,4) + si(1,k,i+1,2)*si(1,k,i+1,4);
            double siza = si(1,k,i,3)*si(1,k,i,4) + si(1,k,i+1,3)*si(1,k,i+1,4);

            double skxa = sk(1,k,i,1)*sk(1,k,i,4) + sk(1,k+1,i,1)*sk(1,k+1,i,4);
            double skya = sk(1,k,i,2)*sk(1,k,i,4) + sk(1,k+1,i,2)*sk(1,k+1,i,4);
            double skza = sk(1,k,i,3)*sk(1,k,i,4) + sk(1,k+1,i,3)*sk(1,k+1,i,4);

            double sjxa = 2.0*sj(1,k,i,1)*sj(1,k,i,4);
            double sjya = 2.0*sj(1,k,i,2)*sj(1,k,i,4);
            double sjza = 2.0*sj(1,k,i,3)*sj(1,k,i,4);

            int ip = std::min(i+1, iend1);
            int im = std::max(i-1, ista);
            double factor = (double)(ip-im);
            factor = ccmaxcr(factor, one_d);
            double rxi = (sj(1,k,ip,1) - sj(1,k,im,1))/factor;
            double ryi = (sj(1,k,ip,2) - sj(1,k,im,2))/factor;
            double rzi = (sj(1,k,ip,3) - sj(1,k,im,3))/factor;
            double pi  = (q(1,k,ip,5)  - q(1,k,im,5)) /factor;

            int kp = std::min(k+1, kend1);
            int km = std::max(k-1, ksta);
            factor = (double)(kp-km);
            factor = ccmaxcr(factor, one_d);
            double rxk = (sj(1,kp,i,1) - sj(1,km,i,1))/factor;
            double ryk = (sj(1,kp,i,2) - sj(1,km,i,2))/factor;
            double rzk = (sj(1,kp,i,3) - sj(1,km,i,3))/factor;
            double pk  = (q(1,kp,i,5)  - q(1,km,i,5)) /factor;

            double sii = sj(1,k,i,1)*sixa + sj(1,k,i,2)*siya + sj(1,k,i,3)*siza;
            double sjj = sj(1,k,i,1)*sjxa + sj(1,k,i,2)*sjya + sj(1,k,i,3)*sjza;
            double skk = sj(1,k,i,1)*skxa + sj(1,k,i,2)*skya + sj(1,k,i,3)*skza;

            double qi2 = qj0(k,i,2,1)*sixa + qj0(k,i,3,1)*siya + qj0(k,i,4,1)*siza;
            double qk2 = qj0(k,i,2,1)*skxa + qj0(k,i,3,1)*skya + qj0(k,i,4,1)*skza;

            double dp = ( (qi2*(qj0(k,i,2,1)*rxi + qj0(k,i,3,1)*ryi +
                               qj0(k,i,4,1)*rzi)
                          +qk2*(qj0(k,i,2,1)*rxk + qj0(k,i,3,1)*ryk +
                               qj0(k,i,4,1)*rzk))*qj0(k,i,1,1)
                         -sii*pi - skk*pk) / sjj;
            int jp2 = std::min(2, jdim1);
            pb = 1.125*pb - 0.125*q(jp2,k,i,5) - dp*0.375;
            pb = ccmax(pb, xminn);
            // now correct density for entropy and velocity for enthalpy
            double sent = q(1,k,i,5) / std::pow(q(1,k,i,1), (double)gamma);
            double hent = ( (double)gamma*q(1,k,i,5)/q(1,k,i,1) )/(double)gm1 +
                          0.5*( q(1,k,i,2)*q(1,k,i,2) + q(1,k,i,3)*q(1,k,i,3)
                               +q(1,k,i,4)*q(1,k,i,4) );
            double rhob  = std::pow(pb/sent, 1.0/(double)gamma);
            double term2 = rhob/qj0(k,i,1,1);
            double vmag  = 2.0*( hent - ((double)gamma*pb/rhob)/(double)gm1 );
            double term1 = std::sqrt( ccabs(vmag) /
                           (qj0(k,i,2,1)*qj0(k,i,2,1)+qj0(k,i,3,1)*qj0(k,i,3,1)
                           +qj0(k,i,4,1)*qj0(k,i,4,1)));
            double uvelb = qj0(k,i,2,1)*term1;
            double vvelb = qj0(k,i,3,1)*term1;
            double wvelb = qj0(k,i,4,1)*term1;

            qj0(k,i,1,1) = rhob;
            qj0(k,i,2,1) = uvelb;
            qj0(k,i,3,1) = vvelb;
            qj0(k,i,4,1) = wvelb;
        }

        qj0(k,i,5,1) = pb;
        bcj(k,i,1)   = 1.0;

        // f23 = 0.0  -  2-point extrapolation
        //       1.0  -  3-point extrapolation
        double f23 = 0.0;

        int j2 = std::min(2, jdim1);
        if (j2 == 1) f23 = 0.0;

        double z1 =   2.0 + 1.5*f23;
        double z2 =       - 0.5*f23;
        double z3 = -(2.0 +     f23);

        qj0(k,i,1,2) = z1*q(1,k,i,1) + z2*q(j2,k,i,1) + z3*qj0(k,i,1,1);
        qj0(k,i,2,2) = z1*q(1,k,i,2) + z2*q(j2,k,i,2) + z3*qj0(k,i,2,1);
        qj0(k,i,3,2) = z1*q(1,k,i,3) + z2*q(j2,k,i,3) + z3*qj0(k,i,3,1);
        qj0(k,i,4,2) = z1*q(1,k,i,4) + z2*q(j2,k,i,4) + z3*qj0(k,i,4,1);
        qj0(k,i,5,2) = z1*q(1,k,i,5) + z2*q(j2,k,i,5) + z3*qj0(k,i,5,1);
    } // k
    } // i

    if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {
            vj0(k,i,1,1) = vist3d(1,k,i);
            vj0(k,i,1,2) = vist3d(1,k,i);
        }
        }
    }
    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal) {
    if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
        for (int l = 1; l <= nummem; l++) {
        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {
            tj0(k,i,l,1) = tursav(1,k,i,l);
            tj0(k,i,l,2) = tursav(1,k,i,l);
        }
        }
        }
    }
    }

    } // end nface==3


// ******************************************************************************
//      j=jdim boundary          inviscid surface                    bctype 1005
// ******************************************************************************
    if (nface == 4) {

    for (int i = ista; i <= iend1; i++) {
    for (int k = ksta; k <= kend1; k++) {
        double pb            = q(jdim1,k,i,5);
        qj0(k,i,1,3)        = q(jdim1,k,i,1);

        double contra        = q(jdim1,k,i,2)*sj(jdim,k,i,1)
                             + q(jdim1,k,i,3)*sj(jdim,k,i,2)
                             + q(jdim1,k,i,4)*sj(jdim,k,i,3) + sj(jdim,k,i,5);

        qj0(k,i,2,3) = q(jdim1,k,i,2) - contra*sj(jdim,k,i,1);
        qj0(k,i,3,3) = q(jdim1,k,i,3) - contra*sj(jdim,k,i,2);
        qj0(k,i,4,3) = q(jdim1,k,i,4) - contra*sj(jdim,k,i,3);

        if (inormmom == 1) {
            double sixa = si(jdim1,k,i,1)*si(jdim1,k,i,4) +
                          si(jdim1,k,i+1,1)*si(jdim1,k,i+1,4);
            double siya = si(jdim1,k,i,2)*si(jdim1,k,i,4) +
                          si(jdim1,k,i+1,2)*si(jdim1,k,i+1,4);
            double siza = si(jdim1,k,i,3)*si(jdim1,k,i,4) +
                          si(jdim1,k,i+1,3)*si(jdim1,k,i+1,4);

            double skxa = sk(jdim1,k,i,1)*sk(jdim1,k,i,4) +
                          sk(jdim1,k+1,i,1)*sk(jdim1,k+1,i,4);
            double skya = sk(jdim1,k,i,2)*sk(jdim1,k,i,4) +
                          sk(jdim1,k+1,i,2)*sk(jdim1,k+1,i,4);
            double skza = sk(jdim1,k,i,3)*sk(jdim1,k,i,4) +
                          sk(jdim1,k+1,i,3)*sk(jdim1,k+1,i,4);

            double sjxa = 2.0*sj(jdim,k,i,1)*sj(jdim,k,i,4);
            double sjya = 2.0*sj(jdim,k,i,2)*sj(jdim,k,i,4);
            double sjza = 2.0*sj(jdim,k,i,3)*sj(jdim,k,i,4);

            int ip = std::min(i+1, iend1);
            int im = std::max(i-1, ista);
            double factor = (double)(ip-im);
            factor = ccmaxcr(factor, one_d);
            double rxi = (sj(jdim,k,ip,1) - sj(jdim,k,im,1))/factor;
            double ryi = (sj(jdim,k,ip,2) - sj(jdim,k,im,2))/factor;
            double rzi = (sj(jdim,k,ip,3) - sj(jdim,k,im,3))/factor;
            double pi  = (q(jdim1,k,ip,5) - q(jdim1,k,im,5))/factor;

            int kp = std::min(k+1, kend1);
            int km = std::max(k-1, ksta);
            factor = (double)(kp-km);
            factor = ccmaxcr(factor, one_d);
            double rxk = (sj(jdim,kp,i,1) - sj(jdim,km,i,1))/factor;
            double ryk = (sj(jdim,kp,i,2) - sj(jdim,km,i,2))/factor;
            double rzk = (sj(jdim,kp,i,3) - sj(jdim,km,i,3))/factor;
            double pk  = (q(jdim1,kp,i,5) - q(jdim1,km,i,5))/factor;

            double sii = sj(jdim,k,i,1)*sixa + sj(jdim,k,i,2)*siya +
                         sj(jdim,k,i,3)*siza;
            double sjj = sj(jdim,k,i,1)*sjxa + sj(jdim,k,i,2)*sjya +
                         sj(jdim,k,i,3)*sjza;
            double skk = sj(jdim,k,i,1)*skxa + sj(jdim,k,i,2)*skya +
                         sj(jdim,k,i,3)*skza;

            double qi2 = qj0(k,i,2,3)*sixa + qj0(k,i,3,3)*siya + qj0(k,i,4,3)*siza;
            double qk2 = qj0(k,i,2,3)*skxa + qj0(k,i,3,3)*skya + qj0(k,i,4,3)*skza;

            double dp = -( (qi2*(qj0(k,i,2,3)*rxi + qj0(k,i,3,3)*ryi +
                                qj0(k,i,4,3)*rzi)
                           +qk2*(qj0(k,i,2,3)*rxk + qj0(k,i,3,3)*ryk +
                                qj0(k,i,4,3)*rzk))*qj0(k,i,1,3)
                          -sii*pi - skk*pk) / sjj;
            int jm2 = std::max(jdim-2, 1);
            pb = 1.125*pb - 0.125*q(jm2,k,i,5) - dp*0.375;
            pb = ccmax(pb, xminn);
            // now correct density for entropy and velocity for enthalpy
            double sent = q(jdim1,k,i,5) / std::pow(q(jdim1,k,i,1), (double)gamma);
            double hent = ( (double)gamma*q(jdim1,k,i,5)/q(jdim1,k,i,1) )/(double)gm1 +
                          0.5*( q(jdim1,k,i,2)*q(jdim1,k,i,2)
                               +q(jdim1,k,i,3)*q(jdim1,k,i,3)
                               +q(jdim1,k,i,4)*q(jdim1,k,i,4) );
            double rhob  = std::pow(pb/sent, 1.0/(double)gamma);
            double term2 = rhob/qj0(k,i,1,3);
            double vmag  = 2.0*( hent - ((double)gamma*pb/rhob)/(double)gm1 );
            double term1 = std::sqrt( ccabs(vmag) /
                           (qj0(k,i,2,3)*qj0(k,i,2,3)+qj0(k,i,3,3)*qj0(k,i,3,3)
                           +qj0(k,i,4,3)*qj0(k,i,4,3)));
            double uvelb = qj0(k,i,2,3)*term1;
            double vvelb = qj0(k,i,3,3)*term1;
            double wvelb = qj0(k,i,4,3)*term1;

            qj0(k,i,1,3) = rhob;
            qj0(k,i,2,3) = uvelb;
            qj0(k,i,3,3) = vvelb;
            qj0(k,i,4,3) = wvelb;
        }

        qj0(k,i,5,3) = pb;
        bcj(k,i,2)   = 1.0;

        // f23 = 0.0  -  2-point extrapolation
        //       1.0  -  3-point extrapolation
        double f23 = 0.0;

        int j2 = std::max(1, jdim-2);
        if (j2 == 1) f23 = 0.0;

        double z1 =  -2.0 - 1.5*f23;
        double z2 =       + 0.5*f23;
        double z3 = +(2.0 +     f23);

        qj0(k,i,1,4) = z1*q(jdim1,k,i,1)+z2*q(j2,k,i,1)+z3*qj0(k,i,1,3);
        qj0(k,i,2,4) = z1*q(jdim1,k,i,2)+z2*q(j2,k,i,2)+z3*qj0(k,i,2,3);
        qj0(k,i,3,4) = z1*q(jdim1,k,i,3)+z2*q(j2,k,i,3)+z3*qj0(k,i,3,3);
        qj0(k,i,4,4) = z1*q(jdim1,k,i,4)+z2*q(j2,k,i,4)+z3*qj0(k,i,4,3);
        qj0(k,i,5,4) = z1*q(jdim1,k,i,5)+z2*q(j2,k,i,5)+z3*qj0(k,i,5,3);
    } // k
    } // i

    if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {
            vj0(k,i,1,3) = vist3d(jdim1,k,i);
            vj0(k,i,1,4) = vist3d(jdim1,k,i);
        }
        }
    }
    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal) {
    if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
        for (int l = 1; l <= nummem; l++) {
        for (int i = ista; i <= iend1; i++) {
        for (int k = ksta; k <= kend1; k++) {
            tj0(k,i,l,3) = tursav(jdim1,k,i,l);
            tj0(k,i,l,4) = tursav(jdim1,k,i,l);
        }
        }
        }
    }
    }

    } // end nface==4



// ******************************************************************************
//      k=1 boundary             inviscid surface                    bctype 1005
// ******************************************************************************
    if (nface == 5) {

    for (int i = ista; i <= iend1; i++) {
    for (int j = jsta; j <= jend1; j++) {
        double pb            = q(j,1,i,5);
        qk0(j,i,1,1)        = q(j,1,i,1);

        double contra        = q(j,1,i,2)*sk(j,1,i,1)
                             + q(j,1,i,3)*sk(j,1,i,2)
                             + q(j,1,i,4)*sk(j,1,i,3) + sk(j,1,i,5);

        qk0(j,i,2,1) = q(j,1,i,2) - contra*sk(j,1,i,1);
        qk0(j,i,3,1) = q(j,1,i,3) - contra*sk(j,1,i,2);
        qk0(j,i,4,1) = q(j,1,i,4) - contra*sk(j,1,i,3);

        if (inormmom == 1) {
            double sixa = si(j,1,i,1)*si(j,1,i,4) + si(j,1,i+1,1)*si(j,1,i+1,4);
            double siya = si(j,1,i,2)*si(j,1,i,4) + si(j,1,i+1,2)*si(j,1,i+1,4);
            double siza = si(j,1,i,3)*si(j,1,i,4) + si(j,1,i+1,3)*si(j,1,i+1,4);

            double sjxa = sj(j,1,i,1)*sj(j,1,i,4) + sj(j+1,1,i,1)*sj(j+1,1,i,4);
            double sjya = sj(j,1,i,2)*sj(j,1,i,4) + sj(j+1,1,i,2)*sj(j+1,1,i,4);
            double sjza = sj(j,1,i,3)*sj(j,1,i,4) + sj(j+1,1,i,3)*sj(j+1,1,i,4);

            double skxa = 2.0*sk(j,1,i,1)*sk(j,1,i,4);
            double skya = 2.0*sk(j,1,i,2)*sk(j,1,i,4);
            double skza = 2.0*sk(j,1,i,3)*sk(j,1,i,4);

            int ip = std::min(i+1, iend1);
            int im = std::max(i-1, ista);
            double factor = (double)(ip-im);
            factor = ccmaxcr(factor, one_d);
            double rxi = (sk(j,1,ip,1) - sk(j,1,im,1))/factor;
            double ryi = (sk(j,1,ip,2) - sk(j,1,im,2))/factor;
            double rzi = (sk(j,1,ip,3) - sk(j,1,im,3))/factor;
            double pi  = (q(j,1,ip,5)  - q(j,1,im,5)) /factor;

            int jp = std::min(j+1, jend1);
            int jm = std::max(j-1, jsta);
            factor = (double)(jp-jm);
            factor = ccmaxcr(factor, one_d);
            double rxj = (sk(jp,1,i,1) - sk(jm,1,i,1))/factor;
            double ryj = (sk(jp,1,i,2) - sk(jm,1,i,2))/factor;
            double rzj = (sk(jp,1,i,3) - sk(jm,1,i,3))/factor;
            double pj  = (q(jp,1,i,5)  - q(jm,1,i,5)) /factor;

            double sii = sk(j,1,i,1)*sixa + sk(j,1,i,2)*siya + sk(j,1,i,3)*siza;
            double sjj = sk(j,1,i,1)*sjxa + sk(j,1,i,2)*sjya + sk(j,1,i,3)*sjza;
            double skk = sk(j,1,i,1)*skxa + sk(j,1,i,2)*skya + sk(j,1,i,3)*skza;

            double qi2 = qk0(j,i,2,1)*sixa + qk0(j,i,3,1)*siya + qk0(j,i,4,1)*siza;
            double qj2 = qk0(j,i,2,1)*sjxa + qk0(j,i,3,1)*sjya + qk0(j,i,4,1)*sjza;

            double dp = ( ( qi2*(qk0(j,i,2,1)*rxi + qk0(j,i,3,1)*ryi +
                                 qk0(j,i,4,1)*rzi)
                           +qj2*(qk0(j,i,2,1)*rxj + qk0(j,i,3,1)*ryj +
                                 qk0(j,i,4,1)*rzj))*qk0(j,i,1,1)
                         -sii*pi - sjj*pj ) / skk;
            int kp2 = std::min(2, kdim1);
            pb = 1.125*pb - 0.125*q(j,kp2,i,5) - dp*0.375;
            pb = ccmax(pb, xminn);
            // now correct density for entropy and velocity for enthalpy
            double sent = q(j,1,i,5) / std::pow(q(j,1,i,1), (double)gamma);
            double hent = ( (double)gamma*q(j,1,i,5)/q(j,1,i,1) )/(double)gm1 +
                          0.5*( q(j,1,i,2)*q(j,1,i,2) + q(j,1,i,3)*q(j,1,i,3)
                               +q(j,1,i,4)*q(j,1,i,4) );
            double rhob  = std::pow(pb/sent, 1.0/(double)gamma);
            double term2 = rhob/qk0(j,i,1,1);
            double vmag  = 2.0*( hent - ((double)gamma*pb/rhob)/(double)gm1 );
            double term1 = std::sqrt( ccabs(vmag) /
                           (qk0(j,i,2,1)*qk0(j,i,2,1)+qk0(j,i,3,1)*qk0(j,i,3,1)
                           +qk0(j,i,4,1)*qk0(j,i,4,1)));
            double uvelb = qk0(j,i,2,1)*term1;
            double vvelb = qk0(j,i,3,1)*term1;
            double wvelb = qk0(j,i,4,1)*term1;

            qk0(j,i,1,1) = rhob;
            qk0(j,i,2,1) = uvelb;
            qk0(j,i,3,1) = vvelb;
            qk0(j,i,4,1) = wvelb;
        }

        qk0(j,i,5,1) = pb;
        bck(j,i,1)   = 1.0;

        // f23 = 0.0  -  2-point extrapolation
        //       1.0  -  3-point extrapolation
        double f23 = 0.0;

        int k2 = std::min(2, kdim1);
        if (k2 == 1) f23 = 0.0;

        double z1 =   2.0 + 1.5*f23;
        double z2 =       - 0.5*f23;
        double z3 = -(2.0 +     f23);

        qk0(j,i,1,2) = z1*q(j,1,i,1) + z2*q(j,k2,i,1) + z3*qk0(j,i,1,1);
        qk0(j,i,2,2) = z1*q(j,1,i,2) + z2*q(j,k2,i,2) + z3*qk0(j,i,2,1);
        qk0(j,i,3,2) = z1*q(j,1,i,3) + z2*q(j,k2,i,3) + z3*qk0(j,i,3,1);
        qk0(j,i,4,2) = z1*q(j,1,i,4) + z2*q(j,k2,i,4) + z3*qk0(j,i,4,1);
        qk0(j,i,5,2) = z1*q(j,1,i,5) + z2*q(j,k2,i,5) + z3*qk0(j,i,5,1);
    } // j
    } // i

    if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {
            vk0(j,i,1,1) = vist3d(j,1,i);
            vk0(j,i,1,2) = vist3d(j,1,i);
        }
        }
    }
    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal) {
    if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
        for (int l = 1; l <= nummem; l++) {
        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {
            tk0(j,i,l,1) = tursav(j,1,i,l);
            tk0(j,i,l,2) = tursav(j,1,i,l);
        }
        }
        }
    }
    }

    } // end nface==5



// ******************************************************************************
//      k=kdim boundary          inviscid surface                    bctype 1005
// ******************************************************************************
    if (nface == 6) {

    for (int i = ista; i <= iend1; i++) {
    for (int j = jsta; j <= jend1; j++) {
        double pb            = q(j,kdim1,i,5);
        qk0(j,i,1,3)        = q(j,kdim1,i,1);

        double contra        = q(j,kdim1,i,2)*sk(j,kdim,i,1)
                             + q(j,kdim1,i,3)*sk(j,kdim,i,2)
                             + q(j,kdim1,i,4)*sk(j,kdim,i,3) + sk(j,kdim,i,5);

        qk0(j,i,2,3) = q(j,kdim1,i,2) - contra*sk(j,kdim,i,1);
        qk0(j,i,3,3) = q(j,kdim1,i,3) - contra*sk(j,kdim,i,2);
        qk0(j,i,4,3) = q(j,kdim1,i,4) - contra*sk(j,kdim,i,3);

        if (inormmom == 1) {
            double sixa = si(j,kdim1,i,1)*si(j,kdim1,i,4) +
                          si(j,kdim1,i+1,1)*si(j,kdim1,i+1,4);
            double siya = si(j,kdim1,i,2)*si(j,kdim1,i,4) +
                          si(j,kdim1,i+1,2)*si(j,kdim1,i+1,4);
            double siza = si(j,kdim1,i,3)*si(j,kdim1,i,4) +
                          si(j,kdim1,i+1,3)*si(j,kdim1,i+1,4);

            double sjxa = sj(j,kdim1,i,1)*sj(j,kdim1,i,4) +
                          sj(j+1,kdim1,i,1)*sj(j+1,kdim1,i,4);
            double sjya = sj(j,kdim1,i,2)*sj(j,kdim1,i,4) +
                          sj(j+1,kdim1,i,2)*sj(j+1,kdim1,i,4);
            double sjza = sj(j,kdim1,i,3)*sj(j,kdim1,i,4) +
                          sj(j+1,kdim1,i,3)*sj(j+1,kdim1,i,4);

            double skxa = 2.0*sk(j,kdim,i,1)*sk(j,kdim,i,4);
            double skya = 2.0*sk(j,kdim,i,2)*sk(j,kdim,i,4);
            double skza = 2.0*sk(j,kdim,i,3)*sk(j,kdim,i,4);

            int ip = std::min(i+1, iend1);
            int im = std::max(i-1, ista);
            double factor = (double)(ip-im);
            factor = ccmaxcr(factor, one_d);
            double rxi = (sk(j,kdim,ip,1) - sk(j,kdim,im,1))/factor;
            double ryi = (sk(j,kdim,ip,2) - sk(j,kdim,im,2))/factor;
            double rzi = (sk(j,kdim,ip,3) - sk(j,kdim,im,3))/factor;
            double pi  = (q(j,kdim1,ip,5) - q(j,kdim1,im,5))/factor;

            int jp = std::min(j+1, jend1);
            int jm = std::max(j-1, jsta);
            factor = (double)(jp-jm);
            factor = ccmaxcr(factor, one_d);
            double rxj = (sk(jp,kdim,i,1) - sk(jm,kdim,i,1))/factor;
            double ryj = (sk(jp,kdim,i,2) - sk(jm,kdim,i,2))/factor;
            double rzj = (sk(jp,kdim,i,3) - sk(jm,kdim,i,3))/factor;
            double pj  = (q(jp,kdim1,i,5) - q(jm,kdim1,i,5))/factor;

            double sii = sk(j,kdim,i,1)*sixa + sk(j,kdim,i,2)*siya +
                         sk(j,kdim,i,3)*siza;
            double sjj = sk(j,kdim,i,1)*sjxa + sk(j,kdim,i,2)*sjya +
                         sk(j,kdim,i,3)*sjza;
            double skk = sk(j,kdim,i,1)*skxa + sk(j,kdim,i,2)*skya +
                         sk(j,kdim,i,3)*skza;

            double qi2 = qk0(j,i,2,3)*sixa + qk0(j,i,3,3)*siya + qk0(j,i,4,3)*siza;
            double qj2 = qk0(j,i,2,3)*sjxa + qk0(j,i,3,3)*sjya + qk0(j,i,4,3)*sjza;

            double dp = -( ( qi2*(qk0(j,i,2,3)*rxi + qk0(j,i,3,3)*ryi +
                                  qk0(j,i,4,3)*rzi)
                            +qj2*(qk0(j,i,2,3)*rxj + qk0(j,i,3,3)*ryj +
                                  qk0(j,i,4,3)*rzj))*qk0(j,i,1,3)
                          -sii*pi - sjj*pj ) / skk;
            int km2 = std::max(kdim-2, 1);
            pb = 1.125*pb - 0.125*q(j,km2,i,5) - dp*0.375;
            pb = ccmax(pb, xminn);
            // now correct density for entropy and velocity for enthalpy
            double sent = q(j,kdim1,i,5) / std::pow(q(j,kdim1,i,1), (double)gamma);
            double hent = ( (double)gamma*q(j,kdim1,i,5)/q(j,kdim1,i,1) )/(double)gm1 +
                          0.5*( q(j,kdim1,i,2)*q(j,kdim1,i,2)
                               +q(j,kdim1,i,3)*q(j,kdim1,i,3)
                               +q(j,kdim1,i,4)*q(j,kdim1,i,4) );
            double rhob  = std::pow(pb/sent, 1.0/(double)gamma);
            double term2 = rhob/qk0(j,i,1,3);
            double vmag  = 2.0*( hent - ((double)gamma*pb/rhob)/(double)gm1 );
            double term1 = std::sqrt( ccabs(vmag) /
                           (qk0(j,i,2,3)*qk0(j,i,2,3)+qk0(j,i,3,3)*qk0(j,i,3,3)
                           +qk0(j,i,4,3)*qk0(j,i,4,3)));
            double uvelb = qk0(j,i,2,3)*term1;
            double vvelb = qk0(j,i,3,3)*term1;
            double wvelb = qk0(j,i,4,3)*term1;

            qk0(j,i,1,3) = rhob;
            qk0(j,i,2,3) = uvelb;
            qk0(j,i,3,3) = vvelb;
            qk0(j,i,4,3) = wvelb;
        }

        qk0(j,i,5,3) = pb;
        bck(j,i,2)   = 1.0;

        // f23 = 0.0  -  2-point extrapolation
        //       1.0  -  3-point extrapolation
        double f23 = 0.0;

        int k2 = std::max(1, kdim-2);
        if (k2 == 1) f23 = 0.0;

        double z1 =  -2.0 - 1.5*f23;
        double z2 =       + 0.5*f23;
        double z3 = +(2.0 +     f23);

        qk0(j,i,1,4) = z1*q(j,kdim1,i,1)+z2*q(j,k2,i,1)+z3*qk0(j,i,1,3);
        qk0(j,i,2,4) = z1*q(j,kdim1,i,2)+z2*q(j,k2,i,2)+z3*qk0(j,i,2,3);
        qk0(j,i,3,4) = z1*q(j,kdim1,i,3)+z2*q(j,k2,i,3)+z3*qk0(j,i,3,3);
        qk0(j,i,4,4) = z1*q(j,kdim1,i,4)+z2*q(j,k2,i,4)+z3*qk0(j,i,4,3);
        qk0(j,i,5,4) = z1*q(j,kdim1,i,5)+z2*q(j,k2,i,5)+z3*qk0(j,i,5,3);
    } // j
    } // i

    if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {
            vk0(j,i,1,3) = vist3d(j,kdim1,i);
            vk0(j,i,1,4) = vist3d(j,kdim1,i);
        }
        }
    }
    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal) {
    if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
        for (int l = 1; l <= nummem; l++) {
        for (int i = ista; i <= iend1; i++) {
        for (int j = jsta; j <= jend1; j++) {
            tk0(j,i,l,3) = tursav(j,kdim1,i,l);
            tk0(j,i,l,4) = tursav(j,kdim1,i,l);
        }
        }
        }
    }
    }

    } // end nface==6



// ******************************************************************************
//      i=1 boundary             inviscid surface                    bctype 1005
// ******************************************************************************
    if (nface == 1) {

    for (int k = ksta; k <= kend1; k++) {
    for (int j = jsta; j <= jend1; j++) {
        double pb            = q(j,k,1,5);
        qi0(j,k,1,1)        = q(j,k,1,1);

        double contra        = q(j,k,1,2)*si(j,k,1,1)
                             + q(j,k,1,3)*si(j,k,1,2)
                             + q(j,k,1,4)*si(j,k,1,3) + si(j,k,1,5);

        qi0(j,k,2,1) = q(j,k,1,2) - contra*si(j,k,1,1);
        qi0(j,k,3,1) = q(j,k,1,3) - contra*si(j,k,1,2);
        qi0(j,k,4,1) = q(j,k,1,4) - contra*si(j,k,1,3);

        if (inormmom == 1) {
            double skxa = sk(j,k,1,1)*sk(j,k,1,4) + sk(j,k+1,1,1)*sk(j,k+1,1,4);
            double skya = sk(j,k,1,2)*sk(j,k,1,4) + sk(j,k+1,1,2)*sk(j,k+1,1,4);
            double skza = sk(j,k,1,3)*sk(j,k,1,4) + sk(j,k+1,1,3)*sk(j,k+1,1,4);

            double sjxa = sj(j,k,1,1)*sj(j,k,1,4) + sj(j+1,k,1,1)*sj(j+1,k,1,4);
            double sjya = sj(j,k,1,2)*sj(j,k,1,4) + sj(j+1,k,1,2)*sj(j+1,k,1,4);
            double sjza = sj(j,k,1,3)*sj(j,k,1,4) + sj(j+1,k,1,3)*sj(j+1,k,1,4);

            double sixa = 2.0*si(j,k,1,1)*si(j,k,1,4);
            double siya = 2.0*si(j,k,1,2)*si(j,k,1,4);
            double siza = 2.0*si(j,k,1,3)*si(j,k,1,4);

            int kp = std::min(k+1, kend1);
            int km = std::max(k-1, ksta);
            double factor = (double)(kp-km);
            factor = ccmaxcr(factor, one_d);
            double rxk = (si(j,kp,1,1) - si(j,km,1,1))/factor;
            double ryk = (si(j,kp,1,2) - si(j,km,1,2))/factor;
            double rzk = (si(j,kp,1,3) - si(j,km,1,3))/factor;
            double pk  = (q(j,kp,1,5)  - q(j,km,1,5)) /factor;

            int jp = std::min(j+1, jend1);
            int jm = std::max(j-1, jsta);
            factor = (double)(jp-jm);
            factor = ccmaxcr(factor, one_d);
            double rxj = (si(jp,k,1,1) - si(jm,k,1,1))/factor;
            double ryj = (si(jp,k,1,2) - si(jm,k,1,2))/factor;
            double rzj = (si(jp,k,1,3) - si(jm,k,1,3))/factor;
            double pj  = (q(jp,k,1,5)  - q(jm,k,1,5)) /factor;

            double sii = si(j,k,1,1)*sixa + si(j,k,1,2)*siya + si(j,k,1,3)*siza;
            double sjj = si(j,k,1,1)*sjxa + si(j,k,1,2)*sjya + si(j,k,1,3)*sjza;
            double skk = si(j,k,1,1)*skxa + si(j,k,1,2)*skya + si(j,k,1,3)*skza;

            double qk2 = qi0(j,k,2,1)*skxa + qi0(j,k,3,1)*skya + qi0(j,k,4,1)*skza;
            double qj2 = qi0(j,k,2,1)*sjxa + qi0(j,k,3,1)*sjya + qi0(j,k,4,1)*sjza;

            double dp = ( (qk2*(qi0(j,k,2,1)*rxk + qi0(j,k,3,1)*ryk +
                               qi0(j,k,4,1)*rzk)
                          +qj2*(qi0(j,k,2,1)*rxj + qi0(j,k,3,1)*ryj +
                               qi0(j,k,4,1)*rzj))*qi0(j,k,1,1)
                         -skk*pk - sjj*pj) / sii;
            int ip2 = std::min(2, idim1);
            pb = 1.125*pb - 0.125*q(j,k,ip2,5) - dp*0.375;
            pb = ccmax(pb, xminn);
            // now correct density for entropy and velocity for enthalpy
            double sent = q(j,k,1,5) / std::pow(q(j,k,1,1), (double)gamma);
            double hent = ( (double)gamma*q(j,k,1,5)/q(j,k,1,1) )/(double)gm1 +
                          0.5*( q(j,k,1,2)*q(j,k,1,2) + q(j,k,1,3)*q(j,k,1,3)
                               +q(j,k,1,4)*q(j,k,1,4) );
            double rhob  = std::pow(pb/sent, 1.0/(double)gamma);
            double term2 = rhob/qi0(j,k,1,1);
            double vmag  = 2.0*( hent - ((double)gamma*pb/rhob)/(double)gm1 );
            double term1 = std::sqrt( ccabs(vmag) /
                           (qi0(j,k,2,1)*qi0(j,k,2,1)+qi0(j,k,3,1)*qi0(j,k,3,1)
                           +qi0(j,k,4,1)*qi0(j,k,4,1)));
            double uvelb = qi0(j,k,2,1)*term1;
            double vvelb = qi0(j,k,3,1)*term1;
            double wvelb = qi0(j,k,4,1)*term1;

            qi0(j,k,1,1) = rhob;
            qi0(j,k,2,1) = uvelb;
            qi0(j,k,3,1) = vvelb;
            qi0(j,k,4,1) = wvelb;
        }

        qi0(j,k,5,1) = pb;
        bci(j,k,1)   = 1.0;

        // f23 = 0.0  -  2-point extrapolation
        //       1.0  -  3-point extrapolation
        double f23 = 0.0;

        int i2 = std::min(2, idim1);
        if (i2 == 1) f23 = 0.0;

        double z1 =   2.0 + 1.5*f23;
        double z2 =       - 0.5*f23;
        double z3 = -(2.0 +     f23);

        qi0(j,k,1,2) = z1*q(j,k,1,1) + z2*q(j,k,i2,1) + z3*qi0(j,k,1,1);
        qi0(j,k,2,2) = z1*q(j,k,1,2) + z2*q(j,k,i2,2) + z3*qi0(j,k,2,1);
        qi0(j,k,3,2) = z1*q(j,k,1,3) + z2*q(j,k,i2,3) + z3*qi0(j,k,3,1);
        qi0(j,k,4,2) = z1*q(j,k,1,4) + z2*q(j,k,i2,4) + z3*qi0(j,k,4,1);
        qi0(j,k,5,2) = z1*q(j,k,1,5) + z2*q(j,k,i2,5) + z3*qi0(j,k,5,1);
    } // j
    } // k

    if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
        for (int k = ksta; k <= kend1; k++) {
        for (int j = jsta; j <= jend1; j++) {
            vi0(j,k,1,1) = vist3d(j,k,1);
            vi0(j,k,1,2) = vist3d(j,k,1);
        }
        }
    }
    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal) {
    if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
        for (int l = 1; l <= nummem; l++) {
        for (int k = ksta; k <= kend1; k++) {
        for (int j = jsta; j <= jend1; j++) {
            ti0(j,k,l,1) = tursav(j,k,1,l);
            ti0(j,k,l,2) = tursav(j,k,1,l);
        }
        }
        }
    }
    }

    } // end nface==1



// ******************************************************************************
//      i=idim boundary          inviscid surface                    bctype 1005
// ******************************************************************************
    if (nface == 2) {

    for (int k = ksta; k <= kend1; k++) {
    for (int j = jsta; j <= jend1; j++) {
        double pb            = q(j,k,idim1,5);
        qi0(j,k,1,3)        = q(j,k,idim1,1);

        double contra        = q(j,k,idim1,2)*si(j,k,idim,1)
                             + q(j,k,idim1,3)*si(j,k,idim,2)
                             + q(j,k,idim1,4)*si(j,k,idim,3) + si(j,k,idim,5);

        qi0(j,k,2,3) = q(j,k,idim1,2) - contra*si(j,k,idim,1);
        qi0(j,k,3,3) = q(j,k,idim1,3) - contra*si(j,k,idim,2);
        qi0(j,k,4,3) = q(j,k,idim1,4) - contra*si(j,k,idim,3);

        if (inormmom == 1) {
            double skxa = sk(j,k,idim1,1)*sk(j,k,idim1,4) +
                          sk(j,k+1,idim1,1)*sk(j,k+1,idim1,4);
            double skya = sk(j,k,idim1,2)*sk(j,k,idim1,4) +
                          sk(j,k+1,idim1,2)*sk(j,k+1,idim1,4);
            double skza = sk(j,k,idim1,3)*sk(j,k,idim1,4) +
                          sk(j,k+1,idim1,3)*sk(j,k+1,idim1,4);

            double sjxa = sj(j,k,idim1,1)*sj(j,k,idim1,4) +
                          sj(j+1,k,idim1,1)*sj(j+1,k,idim1,4);
            double sjya = sj(j,k,idim1,2)*sj(j,k,idim1,4) +
                          sj(j+1,k,idim1,2)*sj(j+1,k,idim1,4);
            double sjza = sj(j,k,idim1,3)*sj(j,k,idim1,4) +
                          sj(j+1,k,idim1,3)*sj(j+1,k,idim1,4);

            double sixa = 2.0*si(j,k,idim,1)*si(j,k,idim,4);
            double siya = 2.0*si(j,k,idim,2)*si(j,k,idim,4);
            double siza = 2.0*si(j,k,idim,3)*si(j,k,idim,4);

            int kp = std::min(k+1, kend1);
            int km = std::max(k-1, ksta);
            double factor = (double)(kp-km);
            factor = ccmaxcr(factor, one_d);
            double rxk = (si(j,kp,idim,1) - si(j,km,idim,1))/factor;
            double ryk = (si(j,kp,idim,2) - si(j,km,idim,2))/factor;
            double rzk = (si(j,kp,idim,3) - si(j,km,idim,3))/factor;
            double pk  = (q(j,kp,idim1,5) - q(j,km,idim1,5))/factor;

            int jp = std::min(j+1, jend1);
            int jm = std::max(j-1, jsta);
            factor = (double)(jp-jm);
            factor = ccmaxcr(factor, one_d);
            double rxj = (si(jp,k,idim,1) - si(jm,k,idim,1))/factor;
            double ryj = (si(jp,k,idim,2) - si(jm,k,idim,2))/factor;
            double rzj = (si(jp,k,idim,3) - si(jm,k,idim,3))/factor;
            double pj  = (q(jp,k,idim1,5) - q(jm,k,idim1,5))/factor;

            double sii = si(j,k,idim,1)*sixa + si(j,k,idim,2)*siya +
                         si(j,k,idim,3)*siza;
            double sjj = si(j,k,idim,1)*sjxa + si(j,k,idim,2)*sjya +
                         si(j,k,idim,3)*sjza;
            double skk = si(j,k,idim,1)*skxa + si(j,k,idim,2)*skya +
                         si(j,k,idim,3)*skza;

            double qk2 = qi0(j,k,2,3)*skxa + qi0(j,k,3,3)*skya + qi0(j,k,4,3)*skza;
            double qj2 = qi0(j,k,2,3)*sjxa + qi0(j,k,3,3)*sjya + qi0(j,k,4,3)*sjza;

            double dp = -( (qk2*(qi0(j,k,2,3)*rxk + qi0(j,k,3,3)*ryk +
                                qi0(j,k,4,3)*rzk)
                           +qj2*(qi0(j,k,2,3)*rxj + qi0(j,k,3,3)*ryj +
                                qi0(j,k,4,3)*rzj))*qi0(j,k,1,3)
                          -skk*pk - sjj*pj) / sii;
            int im2 = std::max(idim-2, 1);
            pb = 1.125*pb - 0.125*q(j,k,im2,5) - dp*0.375;
            pb = ccmax(pb, xminn);
            // now correct density for entropy and velocity for enthalpy
            double sent = q(j,k,idim1,5) / std::pow(q(j,k,idim1,1), (double)gamma);
            double hent = ( (double)gamma*q(j,k,idim1,5)/q(j,k,idim1,1) )/(double)gm1 +
                          0.5*( q(j,k,idim1,2)*q(j,k,idim1,2)
                               +q(j,k,idim1,3)*q(j,k,idim1,3)
                               +q(j,k,idim1,4)*q(j,k,idim1,4) );
            double rhob  = std::pow(pb/sent, 1.0/(double)gamma);
            double term2 = rhob/qi0(j,k,1,3);
            double vmag  = 2.0*( hent - ((double)gamma*pb/rhob)/(double)gm1 );
            double term1 = std::sqrt( ccabs(vmag) /
                           (qi0(j,k,2,3)*qi0(j,k,2,3)+qi0(j,k,3,3)*qi0(j,k,3,3)
                           +qi0(j,k,4,3)*qi0(j,k,4,3)));
            double uvelb = qi0(j,k,2,3)*term1;
            double vvelb = qi0(j,k,3,3)*term1;
            double wvelb = qi0(j,k,4,3)*term1;

            qi0(j,k,1,3) = rhob;
            qi0(j,k,2,3) = uvelb;
            qi0(j,k,3,3) = vvelb;
            qi0(j,k,4,3) = wvelb;
        }

        qi0(j,k,5,3) = pb;
        bci(j,k,2)   = 1.0;

        // f23 = 0.0  -  2-point extrapolation
        //       1.0  -  3-point extrapolation
        double f23 = 0.0;

        int i2 = std::max(1, idim-2);
        if (i2 == 1) f23 = 0.0;

        double z1 =  -2.0 - 1.5*f23;
        double z2 =       + 0.5*f23;
        double z3 = +(2.0 +     f23);

        qi0(j,k,1,4) = z1*q(j,k,idim1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,3);
        qi0(j,k,2,4) = z1*q(j,k,idim1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,3);
        qi0(j,k,3,4) = z1*q(j,k,idim1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,3);
        qi0(j,k,4,4) = z1*q(j,k,idim1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,3);
        qi0(j,k,5,4) = z1*q(j,k,idim1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,3);
    } // j
    } // k

    if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
        for (int k = ksta; k <= kend1; k++) {
        for (int j = jsta; j <= jend1; j++) {
            vi0(j,k,1,3) = vist3d(j,k,idim1);
            vi0(j,k,1,4) = vist3d(j,k,idim1);
        }
        }
    }
    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal) {
    if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
        for (int l = 1; l <= nummem; l++) {
        for (int k = ksta; k <= kend1; k++) {
        for (int j = jsta; j <= jend1; j++) {
            ti0(j,k,l,3) = tursav(j,k,idim1,l);
            ti0(j,k,l,4) = tursav(j,k,idim1,l);
        }
        }
        }
    }
    }
    } // end nface==2

} // end bc1005


} // namespace bc1005_ns
