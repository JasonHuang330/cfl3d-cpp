// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "module_stm_2005.h"
#include "runtime/fortran_array.h"
#include <cmath>

namespace module_stm_2005_ns {

// Private helper: rotate the tensor.
// Tout = R * Tin * R_t
// (rotates the tensor itself, not the reference frame)
static void tensor_rotate(FortranArray1DRef<float> tin,
                          FortranArray2DRef<float> rn,
                          FortranArray2DRef<float> rnt,
                          FortranArray1DRef<float> tout)
{
    // Build symmetric 3x3 tensor from 6-component Voigt notation
    // Fortran: t33(row,col), stored column-major
    // t33(1,1)=tin(1), t33(2,2)=tin(2), t33(3,3)=tin(3)
    // t33(1,2)=tin(4), t33(2,3)=tin(5), t33(1,3)=tin(6)
    // t33(2,1)=t33(1,2), t33(3,2)=t33(2,3), t33(3,1)=t33(1,3)
    float t33[3][3], s[3][3], v[3][3];

    t33[0][0] = tin(1);
    t33[1][1] = tin(2);
    t33[2][2] = tin(3);

    t33[0][1] = tin(4);  // t33(1,2) -> row=1,col=2 -> [col-1][row-1] = [1][0]
    t33[1][2] = tin(5);  // t33(2,3) -> [2][1]
    t33[0][2] = tin(6);  // t33(1,3) -> [2][0]

    t33[1][0] = t33[0][1];  // t33(2,1) = t33(1,2)
    t33[2][1] = t33[1][2];  // t33(3,2) = t33(2,3)
    t33[2][0] = t33[0][2];  // t33(3,1) = t33(1,3)

    // s(:,1) = t33(:,1)*rnt(1,1) + t33(:,2)*rnt(2,1) + t33(:,3)*rnt(3,1)
    // s(:,2) = t33(:,1)*rnt(1,2) + t33(:,2)*rnt(2,2) + t33(:,3)*rnt(3,2)
    // s(:,3) = t33(:,1)*rnt(1,3) + t33(:,2)*rnt(2,3) + t33(:,3)*rnt(3,3)
    // In C: s[col][row] = sum_m t33[m][row] * rnt(m+1, col+1)
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            s[col][row] = t33[0][row] * rnt(1, col+1)
                        + t33[1][row] * rnt(2, col+1)
                        + t33[2][row] * rnt(3, col+1);
        }
    }

    // v(:,1) = rn(:,1)*s(1,1) + rn(:,2)*s(2,1) + rn(:,3)*s(3,1)
    // v(:,2) = rn(:,1)*s(1,2) + rn(:,2)*s(2,2) + rn(:,3)*s(3,2)
    // v(:,3) = rn(:,1)*s(1,3) + rn(:,2)*s(2,3) + rn(:,3)*s(3,3)
    // In C: v[col][row] = sum_m rn(row+1, m+1) * s[col][m]
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            v[col][row] = rn(row+1, 1) * s[col][0]
                        + rn(row+1, 2) * s[col][1]
                        + rn(row+1, 3) * s[col][2];
        }
    }

    // tout(1)=v(1,1), tout(2)=v(2,2), tout(3)=v(3,3)
    // tout(4)=v(1,2), tout(5)=v(2,3), tout(6)=v(1,3)
    tout(1) = v[0][0];  // v(1,1)
    tout(2) = v[1][1];  // v(2,2)
    tout(3) = v[2][2];  // v(3,3)
    tout(4) = v[1][0];  // v(1,2) -> col=2,row=1 -> [1][0]
    tout(5) = v[2][1];  // v(2,3) -> col=3,row=2 -> [2][1]
    tout(6) = v[2][0];  // v(1,3) -> col=3,row=1 -> [2][0]
}

// Build the rotation matrix and its transpose
// based on the rotation angle specified
void stm2k5_get_rotmat(float& thx, float& thy, float& thz,
                       FortranArray2DRef<float> rn,
                       FortranArray2DRef<float> rnt)
{
    if (std::abs(thx) > 0) {
        // rn(:,1) = (/1., 0., 0./)
        rn(1,1) = 1.f;  rn(2,1) = 0.f;  rn(3,1) = 0.f;
        // rn(:,2) = (/0., cos(thx), sin(thx)/)
        rn(1,2) = 0.f;  rn(2,2) = std::cos(thx);  rn(3,2) = std::sin(thx);
        // rn(:,3) = (/0.,-sin(thx), cos(thx)/)
        rn(1,3) = 0.f;  rn(2,3) = -std::sin(thx); rn(3,3) = std::cos(thx);
    } else if (std::abs(thy) > 0) {
        // rn(:,1) = (/cos(thy), 0. , -sin(thy)/)
        rn(1,1) = std::cos(thy);  rn(2,1) = 0.f;  rn(3,1) = -std::sin(thy);
        // rn(:,2) = (/0.,  1.0,  0./)
        rn(1,2) = 0.f;            rn(2,2) = 1.f;  rn(3,2) = 0.f;
        // rn(:,3) = (/sin(thy), 0. ,  cos(thy)/)
        rn(1,3) = std::sin(thy);  rn(2,3) = 0.f;  rn(3,3) = std::cos(thy);
    } else {
        // rn(:,1) = (/ cos(thz), sin(thz), 0./)
        rn(1,1) = std::cos(thz);  rn(2,1) = std::sin(thz);  rn(3,1) = 0.f;
        // rn(:,2) = (/-sin(thz), cos(thz), 0./)
        rn(1,2) = -std::sin(thz); rn(2,2) = std::cos(thz);  rn(3,2) = 0.f;
        // rn(:,3) = (/0., 0., 1.0/)
        rn(1,3) = 0.f;            rn(2,3) = 0.f;             rn(3,3) = 1.f;
    }

    // get the transpose of the rotation matrix
    // rnt(:,1) = rn(1,:)
    rnt(1,1) = rn(1,1);  rnt(2,1) = rn(1,2);  rnt(3,1) = rn(1,3);
    // rnt(:,2) = rn(2,:)
    rnt(1,2) = rn(2,1);  rnt(2,2) = rn(2,2);  rnt(3,2) = rn(2,3);
    // rnt(:,3) = rn(3,:)
    rnt(1,3) = rn(3,1);  rnt(2,3) = rn(3,2);  rnt(3,3) = rn(3,3);
}

// Apply the boundary condition
void stm2k5_bc(FortranArray1DRef<float> vin,
               FortranArray2DRef<float> rn,
               FortranArray2DRef<float> rnt,
               FortranArray1DRef<float> vout)
{
    tensor_rotate(vin, rn, rnt, vout);
    vout(7) = vin(7);
}

} // namespace module_stm_2005_ns
