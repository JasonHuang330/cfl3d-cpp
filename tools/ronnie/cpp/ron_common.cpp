// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ron_common.cpp — definitions of the global COMMON-block structs + unit table.
#include "ron_common.h"

Ginfo ginfo_{};
Info info_{};
Params params_{};
Filenam filenam_{};
Mydist2 mydist2_{};
Igrdtyp igrdtyp_{};
Conversion conversion_{};
Unit5 unit5_{};
Tracer tracer_{};
Zero zero_{};
Sklt1 sklt1_{};
IsDpatch is_dpatch_{};
Save save_{};
Tacos tacos_{};
Areas areas_{};
Hist hist_{};
Degshf degshf_{};
Fvfds fvfds_{};
Ginfo2 ginfo2_{};
Lam lam_{};
Reyue reyue_{};
Wallfun wallfun_{};
Tol tol_{};
Time1 time1_{};
Time2 time2_{};
Rjbdbgi rjbdbgi_{};

FILE* g_units[128] = {nullptr};
