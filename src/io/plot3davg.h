// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace plot3davg_ns {

void plot3d(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray4DRef<double> xv, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, int& nset);
void plot3davg(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray4DRef<double> xv, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, int& nset, FortranArray4DRef<double> qavg, FortranArray4DRef<double> q2avg, int& nt, int& movabs);

} // namespace plot3davg_ns
