// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "module_profileout.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

namespace module_profileout_ns {

// Internal storage for output_parm struct arrays
// (The header declares oparm/cfparm as FortranArray1D<double> for module
//  compatibility, but the actual struct data is stored here.)
static std::vector<output_parm> s_oparm;
static std::vector<output_parm> s_cfparm;

// -----------------------------------------------------------------------
// Private helper: set_default
// -----------------------------------------------------------------------
static void set_default(int& istr, int& iend, int imin, int imax) {
    if (istr == 0) {
        istr = imin;
    }
    if (iend == 0) {
        iend = imax;
    }
}


// -----------------------------------------------------------------------
// profile_read_input
// -----------------------------------------------------------------------
void profile_read_input(char (&filename)[1]) {
    // filename is a Fortran CHARACTER(len=*) — treat as space-terminated C string
    const char* fname = reinterpret_cast<const char*>(&filename[0]);

    // Build trimmed filename string
    std::string fname_str(fname);
    size_t end = fname_str.find_last_not_of(' ');
    if (end == std::string::npos) return; // empty filename
    fname_str = fname_str.substr(0, end + 1);

    // INQUIRE(file=filename, exist=lexist)
    bool lexist = fortran_file_exists(fname_str.c_str());
    if (!lexist) return;

    int ncp_profs = 0;
    int nprofiles = 0;

    // Inner toupper lambda (mirrors contained SUBROUTINE toupper)
    auto toupper_str = [](std::string& s) {
        for (char& c : s) {
            if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
        }
    };

    // First pass: open file and count PROFILE / CP_PROFILE keywords
    fortran_open_unit(1235, fname_str.c_str(), "r");
    {
        FILE* fp = fortran_get_unit(1235);
        char cline[171];
        while (fgets(cline, sizeof(cline), fp) != nullptr) {
            int len = (int)strlen(cline);
            if (len > 0 && cline[len-1] == '\n') { cline[len-1] = '\0'; len--; }

            std::string cline_str(cline);
            // IF(LEN_TRIM(cline)==0) CYCLE
            size_t trim_end = cline_str.find_last_not_of(" \t\r");
            if (trim_end == std::string::npos) continue;

            // str = ADJUSTL(cline); IF(str(1:1)=='!') CYCLE
            size_t trim_start = cline_str.find_first_not_of(" \t");
            std::string str = (trim_start == std::string::npos) ? "" : cline_str.substr(trim_start);
            if (!str.empty() && str[0] == '!') continue;

            // READ(cline,*,iostat=iost1)key,rval — parse first token as key
            std::string key;
            {
                size_t sp = str.find_first_of(" \t");
                key = (sp == std::string::npos) ? str : str.substr(0, sp);
            }
            toupper_str(key);
            size_t ke = key.find_last_not_of(" \t");
            if (ke != std::string::npos) key = key.substr(0, ke + 1);

            if (key == "PROFILE")    nprofiles++;
            if (key == "CP_PROFILE") ncp_profs++;
        }
    }

    if (nprofiles == 0 && ncp_profs == 0) {
        fortran_close_unit(1235);
        return;
    }

    // IF(nprofiles>0) ALLOCATE(oparm(nprofiles))
    if (nprofiles > 0) {
        s_oparm.resize(nprofiles);
        oparm.allocate(nprofiles);
    }
    // IF(ncp_profs>0) ALLOCATE(cfparm(ncp_profs))
    if (ncp_profs > 0) {
        s_cfparm.resize(ncp_profs);
        cfparm.allocate(ncp_profs);
    }

    // REWIND(1235)
    {
        FILE* fp = fortran_get_unit(1235);
        rewind(fp);
    }

    nprofiles = 0;
    ncp_profs = 0;

    // Second pass: read parameters
    {
        FILE* fp = fortran_get_unit(1235);
        char cline[171];
        while (fgets(cline, sizeof(cline), fp) != nullptr) {
            int len = (int)strlen(cline);
            if (len > 0 && cline[len-1] == '\n') { cline[len-1] = '\0'; len--; }

            std::string cline_str(cline);
            size_t trim_end = cline_str.find_last_not_of(" \t\r");
            if (trim_end == std::string::npos) continue;

            size_t trim_start = cline_str.find_first_not_of(" \t");
            std::string str = (trim_start == std::string::npos) ? "" : cline_str.substr(trim_start);
            if (!str.empty() && str[0] == '!') continue;

            std::string key;
            {
                size_t sp = str.find_first_of(" \t");
                key = (sp == std::string::npos) ? str : str.substr(0, sp);
            }
            toupper_str(key);
            size_t ke = key.find_last_not_of(" \t");
            if (ke != std::string::npos) key = key.substr(0, ke + 1);

            if (key == "PROFILE") {
                // READ(cline,*)key,oparm(nprofiles)%name,iblk,istr,iend,jstr,jend,kstr,kend
                nprofiles++;
                output_parm& p = s_oparm[nprofiles - 1];
                char name_buf[41];
                int iblk, istr, iend, jstr, jend, kstr, kend;
                int nread = sscanf(cline, "%*s %40s %d %d %d %d %d %d %d",
                                   name_buf, &iblk, &istr, &iend, &jstr, &jend, &kstr, &kend);
                if (nread >= 8) {
                    std::memset(p.name, ' ', 40);
                    int nlen = (int)strlen(name_buf);
                    if (nlen > 40) nlen = 40;
                    std::memcpy(p.name, name_buf, nlen);
                    p.iblk = iblk; p.istr = istr; p.iend = iend;
                    p.jstr = jstr; p.jend = jend;
                    p.kstr = kstr; p.kend = kend;
                    p.wall_dir = 0; p.wall_str = 0; p.wall_end = 0;
                }
            } else if (key == "CP_PROFILE") {
                // READ(cline,*)key,cfparm(ncp_profs)%name,iblk,...,wall_dir,wall_str,wall_end
                ncp_profs++;
                output_parm& p = s_cfparm[ncp_profs - 1];
                char name_buf[41];
                int iblk, istr, iend, jstr, jend, kstr, kend, wall_dir, wall_str, wall_end;
                int nread = sscanf(cline, "%*s %40s %d %d %d %d %d %d %d %d %d %d",
                                   name_buf, &iblk, &istr, &iend, &jstr, &jend, &kstr, &kend,
                                   &wall_dir, &wall_str, &wall_end);
                if (nread >= 11) {
                    std::memset(p.name, ' ', 40);
                    int nlen = (int)strlen(name_buf);
                    if (nlen > 40) nlen = 40;
                    std::memcpy(p.name, name_buf, nlen);
                    p.iblk = iblk; p.istr = istr; p.iend = iend;
                    p.jstr = jstr; p.jend = jend;
                    p.kstr = kstr; p.kend = kend;
                    p.wall_dir = wall_dir; p.wall_str = wall_str; p.wall_end = wall_end;
                }
            }
        }
    }

    fortran_close_unit(1235);
    return;
}



// -----------------------------------------------------------------------
// profile_plot
// -----------------------------------------------------------------------
void profile_plot(int& jdim, int& kdim, int& idim, int& nummem, int& nblk,
                  FortranArray4DRef<double> q,
                  FortranArray3DRef<double> x,
                  FortranArray3DRef<double> y,
                  FortranArray3DRef<double> z,
                  FortranArray3DRef<double> vist3d,
                  FortranArray4DRef<double> ux,
                  FortranArray4DRef<double> turre,
                  FortranArray3DRef<double> smin) {

    // IF(.NOT.ALLOCATED(oparm)) RETURN
    if (s_oparm.empty()) return;

    int n_oparm = (int)s_oparm.size();
    for (int n = 1; n <= n_oparm; n++) {
        output_parm& p = s_oparm[n - 1];
        // IF(oparm(n)%iblk/=nblk) CYCLE
        if (p.iblk != nblk) continue;

        // Build filename: TRIM(ADJUSTL(oparm(n)%name)) // ".plt"
        char name_buf[41];
        std::memcpy(name_buf, p.name, 40);
        name_buf[40] = '\0';
        int nlen = 40;
        while (nlen > 0 && name_buf[nlen-1] == ' ') nlen--;
        name_buf[nlen] = '\0';
        int nstart = 0;
        while (nstart < nlen && name_buf[nstart] == ' ') nstart++;
        std::string trimmed_name(name_buf + nstart);
        std::string plt_file = trimmed_name + ".plt";

        fortran_open_unit(1234, plt_file.c_str(), "w");

        if (nummem == 7) {
            fortran_write_unit(1234, "%s\n",
                "variables=\"x\",\"y\",\"z\",\"dist\",\"q1\", \"q2\", \"q3\", \"q4\", \"q5\", "
                "\"dudx\",\"dudy\",\"dudz\",\"dvdx\",\"dvdy\",\"dvdz\",\"dwdx\",\"dwdy\",\"dwdz\", "
                "\"mut\", \"t11\", \"t22\", \"t33\", \"t12\", \"t23\", \"t13\", \"zeta/omega/eps\"");
        } else if (nummem == 2) {
            fortran_write_unit(1234, "%s\n",
                "variables=\"x\",\"y\",\"z\",\"dist\",\"q1\", \"q2\", \"q3\", \"q4\", \"q5\", "
                "\"dudx\",\"dudy\",\"dudz\",\"dvdx\",\"dvdy\",\"dvdz\",\"dwdx\",\"dwdy\",\"dwdz\", "
                ",\"mut\", \"omega/eps/zeta\", \"tke\"");
        } else {
            // CLOSE(1234,status='delete')
            fortran_close_unit(1234);
            remove(plt_file.c_str());
            printf("profile_plot: unknown nummem %s %d\n", trimmed_name.c_str(), nummem);
            continue;
        }

        // WRITE(1234,'(A)')'TITLE="'//TRIM(ADJUSTL(oparm(n)%name))//'"'
        fortran_write_unit(1234, "TITLE=\"%s\"\n", trimmed_name.c_str());

        int istr = p.istr;
        int jstr = p.jstr;
        int kstr = p.kstr;
        int iend = p.iend;
        int jend = p.jend;
        int kend = p.kend;

        set_default(istr, iend, 1, idim - 1);
        set_default(jstr, jend, 1, jdim - 1);
        set_default(kstr, kend, 1, kdim - 1);

        // WRITE(1234,'(A)')'ZONE T="'//TRIM(ADJUSTL(oparm(n)%name))
        fortran_write_unit(1234, "ZONE T=\"%s\n", trimmed_name.c_str());

        for (int k = kstr; k <= kend; k++) {
            for (int j = jstr; j <= jend; j++) {
                for (int i = istr; i <= iend; i++) {
                    double xc = 0.125 * (x(j,k,i) + x(j+1,k,i) + x(j,k+1,i) + x(j,k,i+1) +
                                         x(j+1,k+1,i) + x(j+1,k,i+1) + x(j,k+1,i+1) + x(j+1,k+1,i+1));
                    double yc = 0.125 * (y(j,k,i) + y(j+1,k,i) + y(j,k+1,i) + y(j,k,i+1) +
                                         y(j+1,k+1,i) + y(j+1,k,i+1) + y(j,k+1,i+1) + y(j+1,k+1,i+1));
                    double zc = 0.125 * (z(j,k,i) + z(j+1,k,i) + z(j,k+1,i) + z(j,k,i+1) +
                                         z(j+1,k+1,i) + z(j+1,k,i+1) + z(j,k+1,i+1) + z(j+1,k+1,i+1));

                    // WRITE(1234,'(35(1x,es12.5))')xc,yc,zc,smin(j,k,i),q(j,k,i,1:5),
                    //      ux(j,k,i,1:9), vist3d(j,k,i), turre(j,k,i,1:nummem)
                    FILE* fp = fortran_get_unit(1234);
                    fprintf(fp, " %12.5e %12.5e %12.5e %12.5e",
                            xc, yc, zc, smin(j,k,i));
                    for (int m = 1; m <= 5; m++)
                        fprintf(fp, " %12.5e", q(j,k,i,m));
                    for (int m = 1; m <= 9; m++)
                        fprintf(fp, " %12.5e", ux(j,k,i,m));
                    fprintf(fp, " %12.5e", vist3d(j,k,i));
                    for (int m = 1; m <= nummem; m++)
                        fprintf(fp, " %12.5e", turre(j,k,i,m));
                    fprintf(fp, "\n");
                }
            }
        }
        fortran_close_unit(1234);
    }
}



// -----------------------------------------------------------------------
// cfcp_plot
// -----------------------------------------------------------------------
void cfcp_plot(int& jdim, int& kdim, int& idim, int& nummem, int& nblk,
               FortranArray4DRef<double> q,
               FortranArray3DRef<double> x,
               FortranArray3DRef<double> y,
               FortranArray3DRef<double> z,
               FortranArray3DRef<double> vist3d,
               FortranArray4DRef<double> ux,
               FortranArray4DRef<double> turre,
               FortranArray3DRef<double> smin) {

    // COMMON block aliases
    float& pr     = cmn_fluid2.pr;
    float& prt    = cmn_fluid2.prt;
    float& cbar   = cmn_fluid2.cbar;
    float& reue   = cmn_reyue.reue;
    float& tinf   = cmn_reyue.tinf;
    float& gamma_ = cmn_fluid.gamma;
    float& xmach  = cmn_info.xmach;
    float& p0     = cmn_ivals.p0;
    float& rho0   = cmn_ivals.rho0;
    (void)pr; (void)prt; // declared in COMMON but not used in this subroutine

    const float rinf = 1.0f;

    double Re = (double)reue / (double)xmach;

    // IF(.NOT.ALLOCATED(cfparm)) RETURN
    if (s_cfparm.empty()) return;

    int n_cfparm = (int)s_cfparm.size();
    for (int n = 1; n <= n_cfparm; n++) {
        output_parm& p = s_cfparm[n - 1];
        // IF(cfparm(n)%iblk/=nblk) CYCLE
        if (p.iblk != nblk) continue;

        int istr = p.istr;
        int jstr = p.jstr;
        int kstr = p.kstr;
        int iend = p.iend;
        int jend = p.jend;
        int kend = p.kend;

        set_default(istr, iend, 1, idim - 1);
        set_default(jstr, jend, 1, jdim - 1);
        set_default(kstr, kend, 1, kdim - 1);

        // select case(cfparm(n)%wall_dir)
        switch (p.wall_dir) {
        case 1:
            set_default(p.wall_str, p.wall_end, 1, idim - 1);
            break;
        case 2:
            set_default(p.wall_str, p.wall_end, 1, jdim - 1);
            break;
        case 3:
            set_default(p.wall_str, p.wall_end, 1, kdim - 1);
            break;
        }

        // Build filename: TRIM(ADJUSTL(cfparm(n)%name)) // ".plt"
        char name_buf[41];
        std::memcpy(name_buf, p.name, 40);
        name_buf[40] = '\0';
        int nlen = 40;
        while (nlen > 0 && name_buf[nlen-1] == ' ') nlen--;
        name_buf[nlen] = '\0';
        int nstart = 0;
        while (nstart < nlen && name_buf[nstart] == ' ') nstart++;
        std::string trimmed_name(name_buf + nstart);
        std::string plt_file = trimmed_name + ".plt";

        fortran_open_unit(1234, plt_file.c_str(), "w");
        fortran_write_unit(1234, "%s\n", "variables=x,y,z,cf,cp,theta,Re_theta,dstar,Re_dstar");
        // WRITE(1234,'(A)')'ZONE T="'//TRIM(ADJUSTL(cfparm(n)%name))//'"'
        fortran_write_unit(1234, "ZONE T=\"%s\"\n", trimmed_name.c_str());

        for (int i = istr; i <= iend; i++) {
            for (int k = kstr; k <= kend; k++) {
                for (int j = jstr; j <= jend; j++) {
                    double c2b  = (double)cbar / (double)tinf;
                    double c2bp = c2b + 1.0;
                    double tt   = (double)gamma_ * q(j,k,i,5) / q(j,k,i,1);
                    double fnu  = c2bp * tt * std::sqrt(tt) / (c2b + tt);

                    // vel = SIGN(SQRT(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2),q(j,k,i,2))
                    double vel_mag = std::sqrt(q(j,k,i,2)*q(j,k,i,2) +
                                               q(j,k,i,3)*q(j,k,i,3) +
                                               q(j,k,i,4)*q(j,k,i,4));
                    double vel = (q(j,k,i,2) >= 0.0) ? vel_mag : -vel_mag;

                    double tau_w = fnu * std::abs(vel) / std::abs(smin(j,k,i)) / Re;

                    double xc = 0.125 * (x(j,k,i) + x(j+1,k,i) + x(j,k+1,i) + x(j,k,i+1) +
                                         x(j+1,k+1,i) + x(j+1,k,i+1) + x(j,k+1,i+1) + x(j+1,k+1,i+1));
                    double yc = 0.125 * (y(j,k,i) + y(j+1,k,i) + y(j,k+1,i) + y(j,k,i+1) +
                                         y(j+1,k+1,i) + y(j+1,k,i+1) + y(j,k+1,i+1) + y(j+1,k+1,i+1));
                    double zc = 0.125 * (z(j,k,i) + z(j+1,k,i) + z(j,k+1,i) + z(j,k,i+1) +
                                         z(j+1,k+1,i) + z(j+1,k,i+1) + z(j,k+1,i+1) + z(j+1,k+1,i+1));

                    // idx=(/i,j,k/) — Fortran 1-based: idx(1)=i, idx(2)=j, idx(3)=k
                    int idx[3];
                    idx[0] = i;  // idx(1)
                    idx[1] = j;  // idx(2)
                    idx[2] = k;  // idx(3)

                    double dprev  = 0.0;
                    double uprev  = 0.0;
                    double ruprev = 0.0;
                    double theta  = 0.0;
                    double dstar  = 0.0;

                    // idx(cfparm(n)%wall_dir) = cfparm(n)%wall_end
                    idx[p.wall_dir - 1] = p.wall_end;

                    double uinf = 0.0;

                    // SIGN(1, cfparm(n)%wall_end - cfparm(n)%wall_str)
                    int wall_diff = p.wall_end - p.wall_str;
                    int wall_sign = (wall_diff >= 0) ? 1 : -1;

                    // First loop: pick maximum velocity as edge velocity
                    for (int k1 = p.wall_str;
                         (wall_sign > 0) ? (k1 <= p.wall_end) : (k1 >= p.wall_end);
                         k1 += wall_sign) {
                        idx[p.wall_dir - 1] = k1;
                        // q(idx(2),idx(3),idx(1),2) — idx(1)=i-dir, idx(2)=j-dir, idx(3)=k-dir
                        double uinf_cand = q(idx[1], idx[2], idx[0], 2);
                        double vel2 = std::sqrt(q(idx[1],idx[2],idx[0],2)*q(idx[1],idx[2],idx[0],2) +
                                                q(idx[1],idx[2],idx[0],3)*q(idx[1],idx[2],idx[0],3) +
                                                q(idx[1],idx[2],idx[0],4)*q(idx[1],idx[2],idx[0],4));
                        uinf = std::max(uinf, std::max(uinf_cand, vel2));
                    }

                    // Second loop: integrate theta and dstar
                    for (int k1 = p.wall_str;
                         (wall_sign > 0) ? (k1 <= p.wall_end) : (k1 >= p.wall_end);
                         k1 += wall_sign) {
                        idx[p.wall_dir - 1] = k1;

                        double vel2 = std::sqrt(q(idx[1],idx[2],idx[0],2)*q(idx[1],idx[2],idx[0],2) +
                                                q(idx[1],idx[2],idx[0],3)*q(idx[1],idx[2],idx[0],3) +
                                                q(idx[1],idx[2],idx[0],4)*q(idx[1],idx[2],idx[0],4));
                        double rho      = q(idx[1], idx[2], idx[0], 1);
                        double rhou_ave = 0.5 * (rho * vel2 + ruprev);
                        double u_ave    = 0.5 * (vel2 + uprev);

                        double ds = smin(idx[1], idx[2], idx[0]) - dprev;
                        theta += rhou_ave / ((double)rinf * uinf) *
                                 std::max(1.0 - u_ave / uinf, 0.0) * ds;
                        dstar += std::max(1.0 - rhou_ave / ((double)rinf * uinf), 0.0) * ds;

                        dprev  = smin(idx[1], idx[2], idx[0]);
                        uprev  = vel2;
                        ruprev = rho * vel2;

                        // WRITE(*,*)k1,theta,dstar,vel,xmach
                        printf("%d %12.5e %12.5e %12.5e %12.5e\n",
                               k1, theta, dstar, vel2, (double)xmach);
                    }

                    // WRITE(1234,'(14(1x,es12.5))') xc,yc,zc,tau_w/(0.5*rho0*xmach**2),
                    //      (q(j,k,i,5)-p0)/(0.5*rho0*xmach**2),theta,reue*theta,dstar,dstar*reue
                    FILE* fp = fortran_get_unit(1234);
                    fprintf(fp, " %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                            xc, yc, zc,
                            tau_w / (0.5 * (double)rho0 * (double)xmach * (double)xmach),
                            (q(j,k,i,5) - (double)p0) / (0.5 * (double)rho0 * (double)xmach * (double)xmach),
                            theta,
                            (double)reue * theta,
                            dstar,
                            dstar * (double)reue);
                }
            }
        }
        fortran_close_unit(1234);
    }
}


} // namespace module_profileout_ns
