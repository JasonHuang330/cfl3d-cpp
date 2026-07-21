// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "mvdat.h"

namespace mvdat_ns {

void mvdat(int& nbl, int& idir, int& maxbl,
           FortranArray1DRef<double> utrans, FortranArray1DRef<double> vtrans,
           FortranArray1DRef<double> wtrans, FortranArray1DRef<double> omegax,
           FortranArray1DRef<double> omegay, FortranArray1DRef<double> omegaz,
           FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig,
           FortranArray1DRef<double> zorig, FortranArray1DRef<double> dxmx,
           FortranArray1DRef<double> dymx, FortranArray1DRef<double> dzmx,
           FortranArray1DRef<double> dthxmx, FortranArray1DRef<double> dthymx,
           FortranArray1DRef<double> dthzmx, FortranArray1DRef<double> thetax,
           FortranArray1DRef<double> thetay, FortranArray1DRef<double> thetaz,
           FortranArray1DRef<double> rfreqt, FortranArray1DRef<double> rfreqr,
           FortranArray1DRef<double> xorig0, FortranArray1DRef<double> yorig0,
           FortranArray1DRef<double> zorig0, FortranArray1DRef<double> time2,
           FortranArray1DRef<double> thetaxl, FortranArray1DRef<double> thetayl,
           FortranArray1DRef<double> thetazl, FortranArray1DRef<int> itrans,
           FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm,
           FortranArray2DRef<double> utrnsae, FortranArray2DRef<double> vtrnsae,
           FortranArray2DRef<double> wtrnsae, FortranArray2DRef<double> omgxae,
           FortranArray2DRef<double> omgyae, FortranArray2DRef<double> omgzae,
           FortranArray2DRef<double> xorgae, FortranArray2DRef<double> yorgae,
           FortranArray2DRef<double> zorgae, FortranArray2DRef<double> thtxae,
           FortranArray2DRef<double> thtyae, FortranArray2DRef<double> thtzae,
           FortranArray2DRef<double> rfrqtae, FortranArray2DRef<double> rfrqrae,
           FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
           FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
           FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
           FortranArray2DRef<double> freq, FortranArray2DRef<double> gmass,
           FortranArray2DRef<double> damp, FortranArray2DRef<double> x0,
           FortranArray2DRef<double> gf0, int& nmds, int& maxaes,
           FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb,
           FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg,
           FortranArray2DRef<int> iaesurf, int& maxsegdg,
           FortranArray1DRef<double> wk, int& nwork)
{
    // COMMON block aliases
    float& sref      = cmn_fsum.sref;
    float& cref      = cmn_fsum.cref;
    float& bref      = cmn_fsum.bref;
    float& xmc       = cmn_fsum.xmc;
    float& ymc       = cmn_fsum.ymc;
    float& zmc       = cmn_fsum.zmc;

    float& dt        = cmn_info.dt;

    float& xmc0      = cmn_motionmc.xmc0;
    float& ymc0      = cmn_motionmc.ymc0;
    float& zmc0      = cmn_motionmc.zmc0;
    float& utransmc  = cmn_motionmc.utransmc;
    float& vtransmc  = cmn_motionmc.vtransmc;
    float& wtransmc  = cmn_motionmc.wtransmc;
    float& omegaxmc  = cmn_motionmc.omegaxmc;
    float& omegaymc  = cmn_motionmc.omegaymc;
    float& omegazmc  = cmn_motionmc.omegazmc;
    float& xorigmc   = cmn_motionmc.xorigmc;
    float& yorigmc   = cmn_motionmc.yorigmc;
    float& zorigmc   = cmn_motionmc.zorigmc;
    float& xorig0mc  = cmn_motionmc.xorig0mc;
    float& yorig0mc  = cmn_motionmc.yorig0mc;
    float& zorig0mc  = cmn_motionmc.zorig0mc;
    float& thetaxmc  = cmn_motionmc.thetaxmc;
    float& thetaymc  = cmn_motionmc.thetaymc;
    float& thetazmc  = cmn_motionmc.thetazmc;
    float& dxmxmc    = cmn_motionmc.dxmxmc;
    float& dymxmc    = cmn_motionmc.dymxmc;
    float& dzmxmc    = cmn_motionmc.dzmxmc;
    float& dthxmxmc  = cmn_motionmc.dthxmxmc;
    float& dthymxmc  = cmn_motionmc.dthymxmc;
    float& dthzmxmc  = cmn_motionmc.dthzmxmc;
    float& rfreqtmc  = cmn_motionmc.rfreqtmc;
    float& rfreqrmc  = cmn_motionmc.rfreqrmc;
    int32_t& itransmc  = cmn_motionmc.itransmc;
    int32_t& irotatmc  = cmn_motionmc.irotatmc;
    float& time2mc   = cmn_motionmc.time2mc;

    if (idir == 0) {
        wk( 1) = (double)itrans(nbl);
        wk( 2) = rfreqt(nbl);
        wk( 3) = xorig(nbl);
        wk( 4) = yorig(nbl);
        wk( 5) = zorig(nbl);
        wk( 6) = xorig0(nbl);
        wk( 7) = yorig0(nbl);
        wk( 8) = zorig0(nbl);
        wk( 9) = utrans(nbl);
        wk(10) = vtrans(nbl);
        wk(11) = wtrans(nbl);
        wk(12) = dxmx(nbl);
        wk(13) = dymx(nbl);
        wk(14) = dzmx(nbl);
        wk(15) = (double)itransmc;
        wk(16) = (double)rfreqtmc;
        wk(17) = (double)xorigmc;
        wk(18) = (double)yorigmc;
        wk(19) = (double)zorigmc;
        wk(20) = (double)xorig0mc;
        wk(21) = (double)yorig0mc;
        wk(22) = (double)zorig0mc;
        wk(23) = (double)utransmc;
        wk(24) = (double)vtransmc;
        wk(25) = (double)wtransmc;
        wk(26) = (double)xmc;
        wk(27) = (double)ymc;
        wk(28) = (double)zmc;
        wk(29) = (double)dxmxmc;
        wk(30) = (double)dymxmc;
        wk(31) = (double)dzmxmc;
        wk(32) = (double)irotat(nbl);
        wk(33) = rfreqr(nbl);
        wk(34) = thetax(nbl);
        wk(35) = thetay(nbl);
        wk(36) = thetaz(nbl);
        wk(37) = omegax(nbl);
        wk(38) = omegay(nbl);
        wk(39) = omegaz(nbl);
        wk(40) = dthxmx(nbl);
        wk(41) = dthymx(nbl);
        wk(42) = dthzmx(nbl);
        wk(43) = (double)irotatmc;
        wk(44) = (double)rfreqrmc;
        wk(45) = (double)thetaxmc;
        wk(46) = (double)thetaymc;
        wk(47) = (double)thetazmc;
        wk(48) = (double)omegaxmc;
        wk(49) = (double)omegaymc;
        wk(50) = (double)omegazmc;
        wk(51) = (double)dthxmxmc;
        wk(52) = (double)dthymxmc;
        wk(53) = (double)dthzmxmc;
        wk(54) = time2(nbl);
        wk(55) = (double)time2mc;
        wk(56) = (double)dt;
        wk(57) = (double)idefrm(nbl);
        wk(58) = (double)nsegdfrm(nbl);
        int nlast = 58;
        int nseg  = nsegdfrm(nbl);
        for (int is = 1; is <= nseg; is++) {
            wk(is+nlast)          = utrnsae(nbl,is);
            wk(is+nlast+   nseg)  = vtrnsae(nbl,is);
            wk(is+nlast+ 2*nseg)  = wtrnsae(nbl,is);
            wk(is+nlast+ 3*nseg)  = omgxae(nbl,is);
            wk(is+nlast+ 4*nseg)  = omgyae(nbl,is);
            wk(is+nlast+ 5*nseg)  = omgzae(nbl,is);
            wk(is+nlast+ 6*nseg)  = xorgae(nbl,is);
            wk(is+nlast+ 7*nseg)  = yorgae(nbl,is);
            wk(is+nlast+ 8*nseg)  = zorgae(nbl,is);
            wk(is+nlast+ 9*nseg)  = thtxae(nbl,is);
            wk(is+nlast+10*nseg)  = thtyae(nbl,is);
            wk(is+nlast+11*nseg)  = thtzae(nbl,is);
            wk(is+nlast+12*nseg)  = rfrqtae(nbl,is);
            wk(is+nlast+13*nseg)  = rfrqrae(nbl,is);
            wk(is+nlast+14*nseg)  = (double)icsi(nbl,is);
            wk(is+nlast+15*nseg)  = (double)icsf(nbl,is);
            wk(is+nlast+16*nseg)  = (double)jcsi(nbl,is);
            wk(is+nlast+17*nseg)  = (double)jcsf(nbl,is);
            wk(is+nlast+18*nseg)  = (double)kcsi(nbl,is);
            wk(is+nlast+19*nseg)  = (double)kcsf(nbl,is);
        }
    } else {
        itrans(nbl)   = (int)wk( 1);
        rfreqt(nbl)   = wk( 2);
        xorig(nbl)    = wk( 3);
        yorig(nbl)    = wk( 4);
        zorig(nbl)    = wk( 5);
        xorig0(nbl)   = wk( 6);
        yorig0(nbl)   = wk( 7);
        zorig0(nbl)   = wk( 8);
        utrans(nbl)   = wk( 9);
        vtrans(nbl)   = wk(10);
        wtrans(nbl)   = wk(11);
        dxmx(nbl)     = wk(12);
        dymx(nbl)     = wk(13);
        dzmx(nbl)     = wk(14);
        itransmc      = (int32_t)(int)wk(15);
        rfreqtmc      = (float)wk(16);
        xorigmc       = (float)wk(17);
        yorigmc       = (float)wk(18);
        zorigmc       = (float)wk(19);
        xorig0mc      = (float)wk(20);
        yorig0mc      = (float)wk(21);
        zorig0mc      = (float)wk(22);
        utransmc      = (float)wk(23);
        vtransmc      = (float)wk(24);
        wtransmc      = (float)wk(25);
        xmc           = (float)wk(26);
        ymc           = (float)wk(27);
        zmc           = (float)wk(28);
        dxmxmc        = (float)wk(29);
        dymxmc        = (float)wk(30);
        dzmxmc        = (float)wk(31);
        irotat(nbl)   = (int)wk(32);
        rfreqr(nbl)   = wk(33);
        thetax(nbl)   = wk(34);
        thetay(nbl)   = wk(35);
        thetaz(nbl)   = wk(36);
        omegax(nbl)   = wk(37);
        omegay(nbl)   = wk(38);
        omegaz(nbl)   = wk(39);
        dthxmx(nbl)   = wk(40);
        dthymx(nbl)   = wk(41);
        dthzmx(nbl)   = wk(42);
        irotatmc      = (int32_t)(int)wk(43);
        rfreqrmc      = (float)wk(44);
        thetaxmc      = (float)wk(45);
        thetaymc      = (float)wk(46);
        thetazmc      = (float)wk(47);
        omegaxmc      = (float)wk(48);
        omegaymc      = (float)wk(49);
        omegazmc      = (float)wk(50);
        dthxmxmc      = (float)wk(51);
        dthymxmc      = (float)wk(52);
        dthzmxmc      = (float)wk(53);
        time2(nbl)    = wk(54);
        time2mc       = (float)wk(55);
        dt            = (float)wk(56);
        idefrm(nbl)   = (int)wk(57);
        nsegdfrm(nbl) = (int)wk(58);
        int nlast = 58;
        int nseg  = nsegdfrm(nbl);
        for (int is = 1; is <= nseg; is++) {
            utrnsae(nbl,is) = wk(is+nlast);
            vtrnsae(nbl,is) = wk(is+nlast+   nseg);
            wtrnsae(nbl,is) = wk(is+nlast+ 2*nseg);
            omgxae(nbl,is)  = wk(is+nlast+ 3*nseg);
            omgyae(nbl,is)  = wk(is+nlast+ 4*nseg);
            omgzae(nbl,is)  = wk(is+nlast+ 5*nseg);
            xorgae(nbl,is)  = wk(is+nlast+ 6*nseg);
            yorgae(nbl,is)  = wk(is+nlast+ 7*nseg);
            zorgae(nbl,is)  = wk(is+nlast+ 8*nseg);
            thtxae(nbl,is)  = wk(is+nlast+ 9*nseg);
            thtyae(nbl,is)  = wk(is+nlast+10*nseg);
            thtzae(nbl,is)  = wk(is+nlast+11*nseg);
            rfrqtae(nbl,is) = wk(is+nlast+12*nseg);
            rfrqrae(nbl,is) = wk(is+nlast+13*nseg);
            icsi(nbl,is)    = (int)wk(is+nlast+14*nseg);
            icsf(nbl,is)    = (int)wk(is+nlast+15*nseg);
            jcsi(nbl,is)    = (int)wk(is+nlast+16*nseg);
            jcsf(nbl,is)    = (int)wk(is+nlast+17*nseg);
            kcsi(nbl,is)    = (int)wk(is+nlast+18*nseg);
            kcsf(nbl,is)    = (int)wk(is+nlast+19*nseg);
        }
    }

    return;
}

} // namespace mvdat_ns
