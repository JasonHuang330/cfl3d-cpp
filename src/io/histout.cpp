// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "histout.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

namespace histout_ns {

void histout(int& ihstry, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw,
             FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw,
             FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw,
             FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw,
             FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw,
             FortranArray1DRef<double> cmzw, int& n_clcd, FortranArray3DRef<double> clcd,
             int& nblocks_clcd, FortranArray2DRef<int> blocks_clcd,
             FortranArray1DRef<double> chdw, FortranArray1DRef<double> swetw,
             FortranArray1DRef<double> fmdotw, FortranArray1DRef<double> cfttotw,
             FortranArray1DRef<double> cftmomw, FortranArray1DRef<double> cftpw,
             FortranArray1DRef<double> cftvw, FortranArray2DRef<double> rmstr,
             FortranArray2DRef<int> nneg, int& ncycmax, FortranArray4DRef<double> aehist,
             FortranArray2DRef<double> aesrfdat, int& nmds, int& maxaes,
             FortranArray1DRef<double> timekeep, int& nummem)
{
    // Local character variables
    char filename[81];
    char num[3];

    // COMMON block aliases
    float*    title    = cmn_info.title;      // title(20), 1-based → [0..19]
    float&    xmach    = cmn_info.xmach;
    float&    alpha    = cmn_info.alpha;
    float&    dt       = cmn_info.dt;
    int32_t&  ntt      = cmn_info.ntt;
    int32_t&  iclcd    = cmn_avgdata.iclcd;
    float&    radtodeg = cmn_conversion.radtodeg;
    float&    reue     = cmn_reyue.reue;
    float&    time     = cmn_unst.time;
    int32_t&  iunst    = cmn_unst.iunst;
    int32_t&  ialph    = cmn_igrdtyp.ialph;
    int32_t&  ivmx     = cmn_maxiv.ivmx;
    int32_t&  naesrf   = cmn_elastic.naesrf;
    char*     clcdp    = cmn_filenam2.clcdp;  // char[241]

    // Helper: len_trim equivalent for fixed-length char array
    auto len_trim241 = [](const char* s, int maxlen) -> int {
        int len = maxlen;
        while (len > 0 && s[len-1] == ' ') len--;
        return len;
    };

    // alphw = radtodeg * alpha
    float alphw = radtodeg * alpha;

    // Get FILE* for unit 12
    FILE* f12 = fortran_get_unit(12);

    // write(12,2) (real(title(i)),i=1,20)
    // format 2: '# ' followed by 20a4
    // title(20) is float[20]; each float's 4 bytes are printed as 4 characters (a4)
    {
        fprintf(f12, "# ");
        for (int i = 0; i < 20; i++) {
            char buf[5];
            std::memcpy(buf, &title[i], 4);
            buf[4] = '\0';
            fprintf(f12, "%4.4s", buf);
        }
        fprintf(f12, "\n");
    }


    // write(12,'(''# Mach='',e12.4,'', alpha='',e12.4,'', ReUe='',e12.4)')
    //   real(xmach), real(alphw), real(reue)
    fprintf(f12, "# Mach=%12.4e, alpha=%12.4e, ReUe=%12.4e\n",
            (float)xmach, (float)alphw, (float)reue);

    // if (real(dt) .gt. 0) then
    if ((float)dt > 0.0f) {
        // write(12,'(''# Final res='',e12.4,''   Final time (for unsteady)='',f10.4)')
        //   real(rms(ntt)), real(time)
        fprintf(f12, "# Final res=%12.4e   Final time (for unsteady)=%10.4f\n",
                (float)rms(ntt), (float)time);
    } else {
        // write(12,'(''# Final res='',e12.4)') real(rms(ntt))
        fprintf(f12, "# Final res=%12.4e\n", (float)rms(ntt));
    }

    // Branch on ihstry
    if (ihstry == 0) {
        if (ialph == 0) {
            // write(12,'(''# Final cl,cd,cy,cmy='',4e12.4)')
            //   real(clw(ntt)),real(cdw(ntt)),real(cyw(ntt)),real(cmyw(ntt))
            fprintf(f12, "# Final cl,cd,cy,cmy=%12.4e%12.4e%12.4e%12.4e\n",
                    (float)clw(ntt), (float)cdw(ntt), (float)cyw(ntt), (float)cmyw(ntt));

            if (iclcd != 1 && iclcd != 2) {
                // write(12,'(''#'',i5,'' it    log(res)        cl'',
                //   ''            cd            cy            cmy'')') ntt
                fprintf(f12, "#%5d it    log(res)        cl            cd            cy            cmy\n", ntt);
            } else {
                // write(12,'(''variables = it    log(res)        cl'',
                //   ''            cd            cz            cmz'')')
                fprintf(f12, "variables = it    log(res)        cl            cd            cz            cmz\n");
            }

            // do 3000 n=1,ntt
            for (int n = 1; n <= ntt; n++) {
                if ((float)rms(n) <= 0.0f) rms(n) = 1.0;
                // write(12,'(3x,i6,6e15.7)') n,log10(real(rms(n))),
                //   real(clw(n)),real(cdw(n)),real(cyw(n)),real(cmyw(n))
                fprintf(f12, "   %6d%15.7e%15.7e%15.7e%15.7e%15.7e\n",
                        n, (double)std::log10((float)rms(n)),
                        (float)clw(n), (float)cdw(n), (float)cyw(n), (float)cmyw(n));
            }
        } else {
            // write(12,'(''# Final cl,cd,cz,cmz='',4e12.4)')
            //   real(clw(ntt)),real(cdw(ntt)),real(czw(ntt)),real(cmzw(ntt))
            fprintf(f12, "# Final cl,cd,cz,cmz=%12.4e%12.4e%12.4e%12.4e\n",
                    (float)clw(ntt), (float)cdw(ntt), (float)czw(ntt), (float)cmzw(ntt));

            if (iclcd != 1 && iclcd != 2) {
                // write(12,'(''#'',i5,'' it    log(res)        cl'',
                //   ''            cd            cy            cmy'')') ntt
                fprintf(f12, "#%5d it    log(res)        cl            cd            cy            cmy\n", ntt);
            } else {
                // write(12,'(''variables = it    log(res)        cl'',
                //   ''            cd            cz            cmz'')')
                fprintf(f12, "variables = it    log(res)        cl            cd            cz            cmz\n");
            }

            // do 3001 n=1,ntt
            for (int n = 1; n <= ntt; n++) {
                if ((float)rms(n) <= 0.0f) rms(n) = 1.0;
                // write(12,'(3x,i6,6e15.7)') n,log10(real(rms(n))),
                //   real(clw(n)),real(cdw(n)),real(czw(n)),real(cmzw(n))
                fprintf(f12, "   %6d%15.7e%15.7e%15.7e%15.7e%15.7e\n",
                        n, (double)std::log10((float)rms(n)),
                        (float)clw(n), (float)cdw(n), (float)czw(n), (float)cmzw(n));
            }
        }

        // if (iclcd .eq. 1 .or. iclcd .eq. 2) then
        if (iclcd == 1 || iclcd == 2) {
            // do nn = 1, n_clcd
            for (int nn = 1; nn <= n_clcd; nn++) {
                // write(num,'(I2)') nn  → format integer as 2-char string
                std::snprintf(num, sizeof(num), "%2d", nn);
                // num = adjustl(num) → left-trim spaces
                int start = 0;
                while (start < 2 && num[start] == ' ') start++;
                char num_trimmed[3];
                std::strncpy(num_trimmed, num + start, 2 - start);
                num_trimmed[2 - start] = '\0';

                // filename = clcdp(1:len_trim(clcdp)) // num(1:len_trim(num)) // ".dat"
                int clcdp_len = len_trim241(clcdp, 241);
                std::snprintf(filename, sizeof(filename), "%.*s%s.dat",
                              clcdp_len, clcdp, num_trimmed);

                // open(unit=151,file=filename,form="formatted",status="unknown")
                fortran_open_unit(151, filename, "w");
                FILE* f151 = fortran_get_unit(151);

                // write(151,'(''variables=it,Cl,Cd,Clt,Cdt,Cdp,Cdv'')')
                fprintf(f151, "variables=it,Cl,Cd,Clt,Cdt,Cdp,Cdv\n");

                // do n=1,ntt
                for (int n = 1; n <= ntt; n++) {
                    // if (real(clcd(1,nn,n)) < 1e20) then
                    if ((float)clcd(1, nn, n) < 1.0e20f) {
                        // write(151,'(3x,i6,6e15.7)') n,
                        //   real(clcd(1,nn,n)), real(clcd(2,nn,n)),
                        //   real(clw(n)), real(cdw(n)),
                        //   real(cdpw(n)), real(cdvw(n))
                        fprintf(f151, "   %6d%15.7e%15.7e%15.7e%15.7e%15.7e%15.7e\n",
                                n,
                                (float)clcd(1, nn, n), (float)clcd(2, nn, n),
                                (float)clw(n), (float)cdw(n),
                                (float)cdpw(n), (float)cdvw(n));
                    }
                }

                // close(151)
                fortran_close_unit(151);
            }
        }

    } else if (ihstry == 1) {
        // write(12,'(''# Final mass_flow,cftp,cftv,cftmom='',4e12.4)')
        //   real(fmdotw(ntt)),real(cftpw(ntt)),real(cftvw(ntt)),real(cftmomw(ntt))
        fprintf(f12, "# Final mass_flow,cftp,cftv,cftmom=%12.4e%12.4e%12.4e%12.4e\n",
                (float)fmdotw(ntt), (float)cftpw(ntt), (float)cftvw(ntt), (float)cftmomw(ntt));

        if (iclcd != 1 && iclcd != 2) {
            // write(12,'(''#'',i5,'' it    log(res)      mass_flow'',
            //   ''       cftp           cftv        cftmom '')') ntt
            fprintf(f12, "#%5d it    log(res)      mass_flow       cftp           cftv        cftmom \n", ntt);
        } else {
            // write(12,'(i6,'' it    log(res)      mass_flow'',
            //   ''       cftp           cftv        cftmom '')') ntt
            fprintf(f12, "%6d it    log(res)      mass_flow       cftp           cftv        cftmom \n", ntt);
        }

        // do 3100 n=1,ntt
        for (int n = 1; n <= ntt; n++) {
            if ((float)rms(n) <= 0.0f) rms(n) = 1.0;
            // write(12,'(3x,i6,6e15.7)') n,log10(real(rms(n))),
            //   real(fmdotw(n)),real(cftpw(n)),real(cftvw(n)),real(cftmomw(n))
            fprintf(f12, "   %6d%15.7e%15.7e%15.7e%15.7e%15.7e\n",
                    n, (double)std::log10((float)rms(n)),
                    (float)fmdotw(n), (float)cftpw(n), (float)cftvw(n), (float)cftmomw(n));
        }

    } else {
        // ihstry == 2 (enhanced standard): ALL force/moment coefficients
        // write(12,'(''# Final cl,cd       ='',2e13.5)')
        //   real(clw(ntt)),real(cdw(ntt))
        fprintf(f12, "# Final cl,cd       =%13.5e%13.5e\n",
                (float)clw(ntt), (float)cdw(ntt));
        // write(12,'(''# Final cx,cy,cz    ='',3e13.5)')
        //   real(cxw(ntt)),real(cyw(ntt)),real(czw(ntt))
        fprintf(f12, "# Final cx,cy,cz    =%13.5e%13.5e%13.5e\n",
                (float)cxw(ntt), (float)cyw(ntt), (float)czw(ntt));
        // write(12,'(''# Final cmx,cmy,cmz ='',3e13.5)')
        //   real(cmxw(ntt)),real(cmyw(ntt)),real(cmzw(ntt))
        fprintf(f12, "# Final cmx,cmy,cmz =%13.5e%13.5e%13.5e\n",
                (float)cmxw(ntt), (float)cmyw(ntt), (float)cmzw(ntt));
        // write(12,*)
        fprintf(f12, "\n");
        // write(12,'(''#    it     log(res)           cl'',
        //   ''           cd           cx           cy'',
        //   ''           cz          cmx          cmy'',
        //   ''          cmz'')')
        fprintf(f12, "#    it     log(res)           cl           cd           cx           cy           cz          cmx          cmy          cmz\n");
        // write(12,*)
        fprintf(f12, "\n");
        // do n=1,ntt
        for (int n = 1; n <= ntt; n++) {
            if ((float)rms(n) <= 0.0f) rms(n) = 1.0;
            // write(12,'(i6,9e13.5)') n,log10(real(rms(n))),
            //   real(clw(n)),real(cdw(n)),real(cxw(n)),real(cyw(n)),
            //   real(czw(n)),real(cmxw(n)),real(cmyw(n)),real(cmzw(n))
            fprintf(f12, "%6d%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e\n",
                    n, (double)std::log10((float)rms(n)),
                    (float)clw(n), (float)cdw(n), (float)cxw(n), (float)cyw(n),
                    (float)czw(n), (float)cmxw(n), (float)cmyw(n), (float)cmzw(n));
        }
    }


    // output convergence history for field equation turb. model
    // if (ivmx.gt.2) then
    if (ivmx > 2) {
        FILE* f13 = fortran_get_unit(13);

        // write(13,2) (real(title(i)),i=1,20)
        {
            fprintf(f13, "# ");
            for (int i = 0; i < 20; i++) {
                char buf[5];
                std::memcpy(buf, &title[i], 4);
                buf[4] = '\0';
                fprintf(f13, "%4.4s", buf);
            }
            fprintf(f13, "\n");
        }

        // write(13,'(''# Mach='',e12.4,'', alpha='',e12.4,'', ReUe='',e12.4)')
        fprintf(f13, "# Mach=%12.4e, alpha=%12.4e, ReUe=%12.4e\n",
                (float)xmach, (float)alphw, (float)reue);

        // do l=1,nummem: write(13,'(''# Final turres'',i2,''='',e12.4)') l,real(rmstr(ntt,l))
        for (int l = 1; l <= nummem; l++) {
            fprintf(f13, "# Final turres%2d=%12.4e\n", l, (float)rmstr(ntt, l));
        }

        // if (real(dt) .gt. 0) then
        if ((float)dt > 0.0f) {
            fprintf(f13, "# Final time (for unsteady)=%10.4f\n", (float)time);
        }

        if (iclcd != 1 && iclcd != 2) {
            if (nummem == 2) {
                fprintf(f13, "#%6d it  log(turres1)  log(turres2)  nneg1  nneg2\n", ntt);
            } else if (nummem == 3) {
                fprintf(f13, "#%6d it  log(turres1)  log(turres2)  log(turres3)  nneg1  nneg2  nneg3\n", ntt);
            } else if (nummem == 4) {
                fprintf(f13, "#%6d it  log(turres1)  log(turres2)  log(turres3)  log(turres4)  nneg1  nneg2  nneg3  nneg4\n", ntt);
            } else {
                fprintf(f13, "#%6d it  log(turres1)  log(turres2)  log(turres3)  log(turres4)  log(turres5)  log(turres6)  log(turres7)  nneg1  nneg2  nneg3  nneg4  nneg5  nneg6  nneg7\n", ntt);
            }
        } else {
            if (nummem == 2) {
                fprintf(f13, "variables= it  log(turres1)  log(turres2)  nneg1  nneg2\n");
            } else if (nummem == 3) {
                fprintf(f13, "variables= it  log(turres1)  log(turres2)  log(turres3)  nneg1  nneg2  nneg3\n");
            } else if (nummem == 4) {
                fprintf(f13, "variables= it  log(turres1)  log(turres2)  log(turres3)  log(turres4)  nneg1  nneg2  nneg3  nneg4\n");
            } else {
                fprintf(f13, "variables= it  log(turres1)  log(turres2)  log(turres3)  log(turres4)  log(turres5)  log(turres6)  log(turres7)  nneg1  nneg2  nneg3  nneg4  nneg5  nneg6  nneg7\n");
            }
        }

        // do 3200 n=1,ntt
        for (int n = 1; n <= ntt; n++) {
            for (int l = 1; l <= nummem; l++) {
                if ((float)rmstr(n, l) <= 0.0f) rmstr(n, l) = 1.0;
            }
            if (nummem == 2) {
                // write(13,'(3x,i6,2e14.5,1x,i6,1x,i6)')
                fprintf(f13, "   %6d%14.5e%14.5e %6d %6d\n",
                        n,
                        (double)std::log10((float)rmstr(n, 1)),
                        (double)std::log10((float)rmstr(n, 2)),
                        nneg(n, 1), nneg(n, 2));
            } else if (nummem == 3) {
                // write(13,'(3x,i6,3e14.5,1x,i6,1x,i6,1x,i6)')
                fprintf(f13, "   %6d%14.5e%14.5e%14.5e %6d %6d %6d\n",
                        n,
                        (double)std::log10((float)rmstr(n, 1)),
                        (double)std::log10((float)rmstr(n, 2)),
                        (double)std::log10((float)rmstr(n, 3)),
                        nneg(n, 1), nneg(n, 2), nneg(n, 3));
            } else if (nummem == 4) {
                // write(13,'(3x,i6,4e14.5,1x,i6,1x,i6,1x,i6,1x,i6)')
                fprintf(f13, "   %6d%14.5e%14.5e%14.5e%14.5e %6d %6d %6d %6d\n",
                        n,
                        (double)std::log10((float)rmstr(n, 1)),
                        (double)std::log10((float)rmstr(n, 2)),
                        (double)std::log10((float)rmstr(n, 3)),
                        (double)std::log10((float)rmstr(n, 4)),
                        nneg(n, 1), nneg(n, 2), nneg(n, 3), nneg(n, 4));
            } else {
                // write(13,'(3x,i6,7e14.5,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)')
                fprintf(f13, "   %6d%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e %6d %6d %6d %6d %6d %6d %6d\n",
                        n,
                        (double)std::log10((float)rmstr(n, 1)),
                        (double)std::log10((float)rmstr(n, 2)),
                        (double)std::log10((float)rmstr(n, 3)),
                        (double)std::log10((float)rmstr(n, 4)),
                        (double)std::log10((float)rmstr(n, 5)),
                        (double)std::log10((float)rmstr(n, 6)),
                        (double)std::log10((float)rmstr(n, 7)),
                        nneg(n, 1), nneg(n, 2), nneg(n, 3), nneg(n, 4),
                        nneg(n, 5), nneg(n, 6), nneg(n, 7));
            }
        }
    } // end if (ivmx > 2)


    // output history of aeroelastic data
    // if (iunst.gt.1 .and. naesrf.gt.0) then
    if (iunst > 1 && naesrf > 0) {
        FILE* f34 = fortran_get_unit(34);

        // write(34,2) (real(title(i)),i=1,20)
        {
            fprintf(f34, "# ");
            for (int i = 0; i < 20; i++) {
                char buf[5];
                std::memcpy(buf, &title[i], 4);
                buf[4] = '\0';
                fprintf(f34, "%4.4s", buf);
            }
            fprintf(f34, "\n");
        }

        // write(34,'('' Mach='',e12.4,'', alpha='',e12.4,'', ReUe='',e12.4)')
        fprintf(f34, " Mach=%12.4e, alpha=%12.4e, ReUe=%12.4e\n",
                (float)xmach, (float)alphw, (float)reue);

        // write(34,'('' Number of aeroelastic surfaces ='',i3)') naesrf
        fprintf(f34, " Number of aeroelastic surfaces =%3d\n", naesrf);

        // do iaes=1,naesrf
        for (int iaes = 1; iaes <= naesrf; iaes++) {
            // write(34,'('' Data for aeroelastic surface '',i3)') iaes
            fprintf(f34, " Data for aeroelastic surface %3d\n", iaes);

            // nmodes = aesrfdat(5,iaes)
            int nmodes = (int)aesrfdat(5, iaes);

            // do nm=1,nmodes
            for (int nm = 1; nm <= nmodes; nm++) {
                // write(34,'(''   mode number'',i4)') nm
                fprintf(f34, "   mode number%4d\n", nm);

                // write(34,'(''       it      time        xs(2*n-1)'',
                //   ''      xs(2*n)     gforcn(2*n)'')')
                fprintf(f34, "       it      time        xs(2*n-1)      xs(2*n)     gforcn(2*n)\n");

                // do n=1,ntt
                for (int n = 1; n <= ntt; n++) {
                    // write(34,'(3x,i6,6e15.7)') n,real(timekeep(n)),
                    //   real(aehist(n,1,nm,iaes)),real(aehist(n,2,nm,iaes)),
                    //   real(aehist(n,3,nm,iaes))
                    fprintf(f34, "   %6d%15.7e%15.7e%15.7e%15.7e\n",
                            n,
                            (float)timekeep(n),
                            (float)aehist(n, 1, nm, iaes),
                            (float)aehist(n, 2, nm, iaes),
                            (float)aehist(n, 3, nm, iaes));
                }
            }
        }
    } // end if (iunst > 1 && naesrf > 0)

    return;
} // end histout

} // namespace histout_ns
