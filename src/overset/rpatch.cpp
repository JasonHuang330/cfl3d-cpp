// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rpatch.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cstdint>

namespace rpatch_ns {

// Helper: read a Fortran unformatted sequential record from fp.
// Reads the leading 4-byte record marker, then 'nbytes' bytes into buf,
// then the trailing 4-byte record marker.
// Returns true on success, false on EOF/error.
static bool read_fortran_record(FILE* fp, void* buf, int32_t nbytes) {
    int32_t rec_start, rec_end;
    if (fread(&rec_start, sizeof(int32_t), 1, fp) != 1) return false;
    if (fread(buf, 1, (size_t)nbytes, fp) != (size_t)nbytes) return false;
    if (fread(&rec_end, sizeof(int32_t), 1, fp) != 1) return false;
    return true;
}

// Helper: read a single int from a Fortran unformatted record.
static bool read_fortran_int(FILE* fp, int& val) {
    int32_t ival;
    if (!read_fortran_record(fp, &ival, sizeof(int32_t))) return false;
    val = (int)ival;
    return true;
}

// Helper: read a single double from a Fortran unformatted record.
static bool read_fortran_double(FILE* fp, double& val) {
    return read_fortran_record(fp, &val, sizeof(double));
}

// Helper: read an array of ints from a Fortran unformatted record.
static bool read_fortran_int_array(FILE* fp, int* arr, int count) {
    int32_t rec_start, rec_end;
    if (fread(&rec_start, sizeof(int32_t), 1, fp) != 1) return false;
    for (int i = 0; i < count; i++) {
        int32_t ival;
        if (fread(&ival, sizeof(int32_t), 1, fp) != 1) return false;
        arr[i] = (int)ival;
    }
    if (fread(&rec_end, sizeof(int32_t), 1, fp) != 1) return false;
    return true;
}

// Helper: read an array of doubles from a Fortran unformatted record.
static bool read_fortran_double_array(FILE* fp, double* arr, int count) {
    int32_t rec_start, rec_end;
    if (fread(&rec_start, sizeof(int32_t), 1, fp) != 1) return false;
    for (int i = 0; i < count; i++) {
        if (fread(&arr[i], sizeof(double), 1, fp) != 1) return false;
    }
    if (fread(&rec_end, sizeof(int32_t), 1, fp) != 1) return false;
    return true;
}

void rpatch(int& maxbl, int& maxxe, int& intmax, int& nsub1,
            FortranArray2DRef<double> windex, int& ninter,
            FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
            FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy,
            FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> jdimg,
            FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    FILE* fp22 = fortran_get_unit(22);

    // read patched-grid connectivity file
    // read(22) ninter
    read_fortran_int(fp22, ninter);

    if (ninter > intmax) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            " too many blocks to interpolate; intmax should be at least %4d",
            ninter);
        int m1 = -1;
        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
    }

    for (int n = 1; n <= std::abs(ninter); n++) {
        // read(22) iindex(n,1)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 1) = val;
        }
        int nfb = iindex(n, 1);

        // read(22) (iindex(n,1+ll),ll=1,nfb)
        {
            int32_t rec_start, rec_end;
            fread(&rec_start, sizeof(int32_t), 1, fp22);
            for (int ll = 1; ll <= nfb; ll++) {
                int32_t ival;
                fread(&ival, sizeof(int32_t), 1, fp22);
                iindex(n, 1 + ll) = (int)ival;
            }
            fread(&rec_end, sizeof(int32_t), 1, fp22);
        }

        // read(22) iindex(n,nfb+2)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, nfb + 2) = val;
        }

        // read(22) (iindex(n,nfb+2+ll),ll=1,nfb)
        {
            int32_t rec_start, rec_end;
            fread(&rec_start, sizeof(int32_t), 1, fp22);
            for (int ll = 1; ll <= nfb; ll++) {
                int32_t ival;
                fread(&ival, sizeof(int32_t), 1, fp22);
                iindex(n, nfb + 2 + ll) = (int)ival;
            }
            fread(&rec_end, sizeof(int32_t), 1, fp22);
        }

        // read(22) iindex(n,2*nfb+3)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 3) = val;
        }

        // read(22) iindex(n,2*nfb+4)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 4) = val;
        }

        // read(22) iindex(n,2*nfb+5)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 5) = val;
        }

        // read(22) iindex(n,2*nfb+6)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 6) = val;
        }

        // read(22) iindex(n,2*nfb+7)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 7) = val;
        }

        // read(22) iindex(n,2*nfb+8)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 8) = val;
        }

        // read(22) iindex(n,2*nfb+9)
        {
            int val;
            read_fortran_int(fp22, val);
            iindex(n, 2 * nfb + 9) = val;
        }

        int nbl = iindex(n, nfb + 2);
        int lst  = iindex(n, 2 * nfb + 5);
        int len  = lst + iindex(n, 2 * nfb + 4) - 1;

        if (n == ninter && len > maxxe) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "\n stopping...parameter maxxe is too small...must be at least%6d",
                len);
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }

        // read(22) (nblkpt(nnn),nnn=lst,len)
        {
            int count = len - lst + 1;
            int32_t rec_start, rec_end;
            fread(&rec_start, sizeof(int32_t), 1, fp22);
            for (int nnn = lst; nnn <= len; nnn++) {
                int32_t ival;
                fread(&ival, sizeof(int32_t), 1, fp22);
                nblkpt(nnn) = (int)ival;
            }
            fread(&rec_end, sizeof(int32_t), 1, fp22);
        }

        // read(22) ((windex(nnn,ll),nnn=lst,len),ll=1,2)
        {
            int32_t rec_start, rec_end;
            fread(&rec_start, sizeof(int32_t), 1, fp22);
            for (int ll = 1; ll <= 2; ll++) {
                for (int nnn = lst; nnn <= len; nnn++) {
                    double dval;
                    fread(&dval, sizeof(double), 1, fp22);
                    windex(nnn, ll) = dval;
                }
            }
            fread(&rec_end, sizeof(int32_t), 1, fp22);
        }
    } // end do 1500

    ninter = std::abs(ninter);

    //
    // if angular displacements are in the patch file, read them
    // (obtained with ioflag=2 in ronnie input); if not (ioflag=0 or 1)
    // set angular patch displacements to zero
    //
    int itest = 0;
    bool hit_eof = false;
    for (int n = 1; n <= ninter && !hit_eof; n++) {
        int nfb = iindex(n, 1);
        for (int ll = 1; ll <= nfb && !hit_eof; ll++) {
            // read(22,end=999) dthetxx(n,ll),dthetyy(n,ll),dthetzz(n,ll)
            // This reads 3 doubles from one unformatted record
            int32_t rec_start;
            if (fread(&rec_start, sizeof(int32_t), 1, fp22) != 1) {
                hit_eof = true;
                break;
            }
            double vx, vy, vz;
            if (fread(&vx, sizeof(double), 1, fp22) != 1) { hit_eof = true; break; }
            if (fread(&vy, sizeof(double), 1, fp22) != 1) { hit_eof = true; break; }
            if (fread(&vz, sizeof(double), 1, fp22) != 1) { hit_eof = true; break; }
            int32_t rec_end;
            if (fread(&rec_end, sizeof(int32_t), 1, fp22) != 1) { hit_eof = true; break; }
            dthetxx(n, ll) = vx;
            dthetyy(n, ll) = vy;
            dthetzz(n, ll) = vz;
        }
    }
    if (!hit_eof) {
        itest = 1;
    }

    // 999 continue
    if (itest == 0) {
        for (int n = 1; n <= ninter; n++) {
            int nfb = iindex(n, 1);
            for (int ll = 1; ll <= nfb; ll++) {
                dthetxx(n, ll) = 0.;
                dthetyy(n, ll) = 0.;
                dthetzz(n, ll) = 0.;
            }
        }
    }

    //
    // set search range on from side to full index range
    // (added to be compatible with dynamic patch storage
    // - not needed in flow code for static patching)
    //
    for (int n = 1; n <= ninter; n++) {
        int nfb = iindex(n, 1);
        for (int ll = 1; ll <= nfb; ll++) {
            int nbl = iindex(n, ll + 1);
            int itop = (iindex(n, ll + nfb + 2) - iindex(n, ll + nfb + 2) / 100 * 100) / 10;
            int jd1 = 0, kd1 = 0;
            if (itop == 1) {
                jd1 = jdimg(nbl);
                kd1 = kdimg(nbl);
            }
            if (itop == 2) {
                jd1 = kdimg(nbl);
                kd1 = idimg(nbl);
            }
            if (itop == 3) {
                jd1 = jdimg(nbl);
                kd1 = idimg(nbl);
            }
            iindex(n, 2 * nfb + 9 + ll) = 1;
            iindex(n, 3 * nfb + 9 + ll) = jd1;
            iindex(n, 4 * nfb + 9 + ll) = 1;
            iindex(n, 5 * nfb + 9 + ll) = kd1;
        }
    } // end do 1700

    return;
}

} // namespace rpatch_ns
