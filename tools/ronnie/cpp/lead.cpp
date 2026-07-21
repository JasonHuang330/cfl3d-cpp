// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// lead.cpp — install block attributes into common blocks given a block number.
// lw is dimensioned (65,maxbl), lw2 is (43,maxbl); both column-major, 1-based.
#include "ron_common.h"

void lead(int nbl,int* lw,int* lw2,int /*maxbl*/){
    // lw2(i,nbl) at flat (i-1)+(nbl-1)*43 ; lw(i,nbl) at (i-1)+(nbl-1)*65
    auto W2=[&](int i)->int{ return lw2[(i-1)+(nbl-1)*43]; };
    auto W =[&](int i)->int{ return lw [(i-1)+(nbl-1)*65]; };

    ginfo_.jdim = W2(1);
    ginfo_.kdim = W2(2);
    ginfo_.idim = W2(3);
    ginfo_.nblc = W2(4);
    ginfo_.jj2  = W2(5);
    ginfo_.kk2  = W2(6);
    ginfo_.ii2  = W2(7);
    lam_.ilamlo = W2(8);
    lam_.ilamhi = W2(9);
    lam_.jlamlo = W2(10);
    lam_.jlamhi = W2(11);
    lam_.klamlo = W2(12);
    lam_.klamhi = W2(13);
    reyue_.ivisc[0] = W2(14);
    reyue_.ivisc[1] = W2(15);
    reyue_.ivisc[2] = W2(16);
    degshf_.ideg[0] = W2(17);
    degshf_.ideg[1] = W2(18);
    degshf_.ideg[2] = W2(19);
    wallfun_.iwf[0] = W2(20);
    wallfun_.iwf[1] = W2(21);
    wallfun_.iwf[2] = W2(22);
    info_.idiag[0]  = W2(23);
    info_.idiag[1]  = W2(24);
    info_.idiag[2]  = W2(25);
    info_.iflim[0]  = W2(26);
    info_.iflim[1]  = W2(27);
    info_.iflim[2]  = W2(28);
    fvfds_.ifds[0]  = W2(29);
    fvfds_.ifds[1]  = W2(30);
    fvfds_.ifds[2]  = W2(31);
    // rkap0 divides by 1e6 (must be consistent with pointers)
    fvfds_.rkap0[0] = (double)W2(32)/1e6;
    fvfds_.rkap0[1] = (double)W2(33)/1e6;
    fvfds_.rkap0[2] = (double)W2(34)/1e6;
    ginfo_.js = W2(35);
    ginfo_.ks = W2(36);
    ginfo_.is = W2(37);
    ginfo_.je = W2(38);
    ginfo_.ke = W2(39);
    ginfo_.ie = W2(40);

    ginfo_.lq    = W(1);
    ginfo_.lqj0  = W(2);
    ginfo_.lqk0  = W(3);
    ginfo_.lqi0  = W(4);
    ginfo_.lsj   = W(5);
    ginfo_.lsk   = W(6);
    ginfo_.lsi   = W(7);
    ginfo_.lvol  = W(8);
    ginfo_.ldtj  = W(9);
    ginfo_.lx    = W(10);
    ginfo_.ly    = W(11);
    ginfo_.lz    = W(12);
    ginfo_.lvis  = W(13);
    ginfo_.lsnk0 = W(14);
    ginfo_.lsni0 = W(15);
    ginfo_.lq1   = W(16);
    ginfo_.lqr   = W(17);
    ginfo_.lblk  = W(18);
    ginfo_.lxib  = W(19);
    ginfo_.lsig  = W(20);
    ginfo_.lsqtq = W(21);
    ginfo_.lg    = W(22);
    ginfo_.ltj0  = W(23);
    ginfo_.ltk0  = W(24);
    ginfo_.lti0  = W(25);
    ginfo_.lxkb  = W(26);
    ginfo_.lnbl  = W(27);
    ginfo_.lvj0  = W(28);
    ginfo_.lvk0  = W(29);
    ginfo_.lvi0  = W(30);
    ginfo_.lbcj  = W(31);
    ginfo_.lbck  = W(32);
    ginfo_.lbci  = W(33);
    ginfo_.lqc0  = W(34);
    ginfo_.ldqc0 = W(35);
    ginfo_.lxtbj = W(36);
    ginfo_.lxtbk = W(37);
    ginfo_.lxtbi = W(38);
    ginfo_.latbj = W(39);
    ginfo_.latbk = W(40);
    ginfo_.latbi = W(41);
    ginfo_.lbcdj = W(42);
    ginfo_.lbcdk = W(43);
    ginfo_.lbcdi = W(44);
    ginfo_.lxib2 = W(45);
    ginfo_.lqavg = W(46);
    ginfo_.lux   = W(47);
    ginfo_.lcmuv = W(48);
    ginfo_.lvolj0= W(49);
    ginfo_.lvolk0= W(50);
    ginfo_.lvoli0= W(51);
    ginfo_.lxmdj = W(52);
    ginfo_.lxmdk = W(53);
    ginfo_.lxmdi = W(54);
    ginfo_.lvelg = W(55);
    ginfo_.lxnm2 = W(56);
    ginfo_.lynm2 = W(57);
    ginfo_.lznm2 = W(58);
    ginfo_.ldeltj= W(59);
    ginfo_.ldeltk= W(60);
    ginfo_.ldelti= W(61);
    ginfo_.lxnm1 = W(62);
    ginfo_.lynm1 = W(63);
    ginfo_.lznm1 = W(64);
    ginfo2_.lq2avg = W(65);
}
