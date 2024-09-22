import pandas
import numpy
import sys

namemap = {}
namemap["wastrid_vanilla"]                = "               wastrid_vanilla"
namemap["asteroid"]                       = "                      asteroid"
namemap["aster_v1.16.3.4"]                = "               aster_v1.16.3.4"
namemap["wtreeqmc_wf_n2_refined"]         = "        wtreeqmc_wf_n2_refined"
namemap["wtreeqmc_wf_n1_refined"]         = "        wtreeqmc_wf_n1_refined"
namemap["wtreeqmc_wf_n0_refined"]         = "        wtreeqmc_wf_n0_refined"
namemap["wtreeqmc_wf_n2"]                 = "                wtreeqmc_wf_n2"
namemap["wtreeqmc_wf_n1"]                 = "                wtreeqmc_wf_n1"
namemap["wtreeqmc_wf_n0"]                 = "                wtreeqmc_wf_n0"
namemap["wtreeqmc_wf_n2_shared_refined"]  = " wtreeqmc_wf_n2_shared_refined"
namemap["wtreeqmc_wf_n1_shared_refined"]  = " wtreeqmc_wf_n1_shared_refined"
namemap["wtreeqmc_wf_n2_shared"]          = "         wtreeqmc_wf_n2_shared"
namemap["wtreeqmc_wf_n1_shared"]          = "         wtreeqmc_wf_n1_shared"

mthds = ["wastrid_vanilla",
         "asteroid",
         "aster_v1.16.3.4",
         "wtreeqmc_wf_n2_refined",
         "wtreeqmc_wf_n1_refined",
         "wtreeqmc_wf_n0_refined",
         "wtreeqmc_wf_n2",
         "wtreeqmc_wf_n1",
         "wtreeqmc_wf_n0",
         "wtreeqmc_wf_n2_shared_refined",
         "wtreeqmc_wf_n1_shared_refined",
         "wtreeqmc_wf_n2_shared",
         "wtreeqmc_wf_n1_shared"]

experiments = ["varypsiz", "varyntax", "varyngen",
               "varynbps", "varyblsc", "varymiss"]

for do in experiments:
    ntaxs = [50]
    ngens = [1000]
    nbpss = [100]
    blscs = [1.0]
    psizs = [50000000]
    misss = [0.6]  # ms is always same as mf
    repls = [repl for repl in range(0, 51)]

    if do == "varypsiz":
        psizs = [10, 50000000, 100000000, 500000000, 1000000000]
    elif do == "varyntax":
        ntaxs = [25, 75, 50, 100, 125, 150]  # dup 50
    elif do == "varyngen":
        ngens = [250, 500, 1000, 2000]  # dup 1000
    elif do == "varynbps":
        nbpss = [50, 100, 200, 500]  # dup 100
    elif do == "varyblsc":
        blscs = [0.05, 0.1, 1, 10, 100, 200]  # dup 1.0
    elif do == "varymiss":
        misss = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75]  # dup 0.6
    else:
        sys.exit()

    df = pandas.read_csv("data-" + do + "-error-and-timings.csv",
                         na_values='NA', keep_default_na=False)
    sys.stdout.write("Reporting experitment %s\n" % do)

    for ntax in ntaxs:
        for ngen in ngens:
            for nbps in nbpss:
                for blsc in blscs:
                    for psiz in psizs:
                        for miss in misss:
                            print("  %d taxa, %d genes, %d bps, %f bl scaler, %d pop size, %f missing param" \
                                  % (ntax, ngen, nbps, blsc, psiz, miss))

                            xdf = df[(df["NTAX"] == ntax) &
                                     (df["NGEN"] == ngen) &
                                     (df["NBPS"] == nbps) &
                                     (df["BLSC"] == blsc) &
                                     (df["PSIZ"] == psiz) & 
                                     (df["MISS"] == miss)]

                            for mthd in mthds:
                                name = namemap[mthd]
                                ydf = xdf[(xdf["MTHD"] == mthd)]

                                data = ydf.SEFN.values
                                if data.size != 50:
                                    sys.stdout.write("ERROR - wrong number of replicates!\n")

                                nrepl = numpy.sum(~numpy.isnan(data))

                                fnavg = numpy.mean(ydf.SEFNR.values)
                                fnavg = numpy.round(numpy.round(fnavg, 4), 3)

                                fpavg = numpy.mean(ydf.SEFPR.values)
                                fpavg = numpy.round(numpy.round(fpavg, 4), 3)

                                if (mthd == "wtreeqmc_wf_n2_refined") or (mthd == "aster_v1.16.3.4"):
                                    qsavg = numpy.mean(ydf.QSCR.values)
                                    qsavg = numpy.round(numpy.round(qsavg, 4), 3)

                                    lppavg = numpy.nanmean(ydf.AVG_LPP.values)
                                    lppavg = numpy.round(numpy.round(lppavg, 4), 3)

                                    print("    %s : %1.3f (fn) %1.3f (fp) %1.3f (qs) %1.3f (lpp)" % (name, fnavg, fpavg, qsavg, lppavg))
                                else:
                                    print("    %s : %1.3f (fn) %1.3f (fp)" % (name, fnavg, fpavg))

                            print()

"""
Reporting experitment varypsiz
  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 10 pop size, 0.600000 missing param
           wastrid_vanilla : 0.212 (fn) 0.212 (fp)
           aster_v1.16.3.4 : 0.177 (fn) 0.177 (fp) 346992.340 (qs) 0.851 (lpp)
                  asteroid : 0.144 (fn) 0.144 (fp)
    wtreeqmc_wf_n2_refined : 0.154 (fn) 0.154 (fp) 346695.400 (qs) 0.858 (lpp)
            wtreeqmc_wf_n2 : 0.158 (fn) 0.144 (fp)
     wtreeqmc_wf_n2_shared : 0.226 (fn) 0.215 (fp)
            wtreeqmc_wf_n1 : 0.161 (fn) 0.146 (fp)
     wtreeqmc_wf_n1_shared : 0.215 (fn) 0.202 (fp)
            wtreeqmc_wf_n0 : 0.189 (fn) 0.176 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.246 (fn) 0.246 (fp)
           aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp) 349514.440 (qs) 0.836 (lpp)
                  asteroid : 0.181 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_refined : 0.191 (fn) 0.191 (fp) 349295.820 (qs) 0.840 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
            wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
     wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
            wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 100000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.297 (fn) 0.297 (fp)
           aster_v1.16.3.4 : 0.262 (fn) 0.262 (fp) 335223.460 (qs) 0.810 (lpp)
                  asteroid : 0.232 (fn) 0.232 (fp)
    wtreeqmc_wf_n2_refined : 0.235 (fn) 0.235 (fp) 334825.580 (qs) 0.819 (lpp)
            wtreeqmc_wf_n2 : 0.240 (fn) 0.224 (fp)
     wtreeqmc_wf_n2_shared : 0.334 (fn) 0.321 (fp)
            wtreeqmc_wf_n1 : 0.243 (fn) 0.227 (fp)
     wtreeqmc_wf_n1_shared : 0.310 (fn) 0.297 (fp)
            wtreeqmc_wf_n0 : 0.274 (fn) 0.259 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 500000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.475 (fn) 0.475 (fp)
           aster_v1.16.3.4 : 0.462 (fn) 0.462 (fp) 271430.280 (qs) 0.738 (lpp)
                  asteroid : 0.434 (fn) 0.434 (fp)
    wtreeqmc_wf_n2_refined : 0.422 (fn) 0.422 (fp) 270790.620 (qs) 0.740 (lpp)
            wtreeqmc_wf_n2 : 0.424 (fn) 0.416 (fp)
     wtreeqmc_wf_n2_shared : 0.524 (fn) 0.518 (fp)
            wtreeqmc_wf_n1 : 0.430 (fn) 0.421 (fp)
     wtreeqmc_wf_n1_shared : 0.516 (fn) 0.509 (fp)
            wtreeqmc_wf_n0 : 0.461 (fn) 0.452 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 1000000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.585 (fn) 0.585 (fp)
           aster_v1.16.3.4 : 0.588 (fn) 0.588 (fp) 245781.540 (qs) 0.693 (lpp)
                  asteroid : 0.553 (fn) 0.553 (fp)
    wtreeqmc_wf_n2_refined : 0.560 (fn) 0.560 (fp) 244972.440 (qs) 0.691 (lpp)
            wtreeqmc_wf_n2 : 0.561 (fn) 0.553 (fp)
     wtreeqmc_wf_n2_shared : 0.669 (fn) 0.662 (fp)
            wtreeqmc_wf_n1 : 0.562 (fn) 0.554 (fp)
     wtreeqmc_wf_n1_shared : 0.647 (fn) 0.640 (fp)
            wtreeqmc_wf_n0 : 0.610 (fn) 0.603 (fp)

Reporting experitment varyntax
  25 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.269 (fn) 0.269 (fp)
           aster_v1.16.3.4 : 0.276 (fn) 0.276 (fp) 16935.940 (qs) 0.813 (lpp)
                  asteroid : 0.226 (fn) 0.226 (fp)
    wtreeqmc_wf_n2_refined : 0.234 (fn) 0.234 (fp) 16896.600 (qs) 0.822 (lpp)
            wtreeqmc_wf_n2 : 0.242 (fn) 0.223 (fp)
     wtreeqmc_wf_n2_shared : 0.288 (fn) 0.270 (fp)
            wtreeqmc_wf_n1 : 0.240 (fn) 0.222 (fp)
     wtreeqmc_wf_n1_shared : 0.284 (fn) 0.267 (fp)
            wtreeqmc_wf_n0 : 0.281 (fn) 0.265 (fp)

  75 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.244 (fn) 0.244 (fp)
           aster_v1.16.3.4 : 0.209 (fn) 0.209 (fp) 1740392.400 (qs) 0.841 (lpp)
                  asteroid : 0.183 (fn) 0.183 (fp)
    wtreeqmc_wf_n2_refined : 0.189 (fn) 0.189 (fp) 1739341.140 (qs) 0.841 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.284 (fn) 0.269 (fp)
            wtreeqmc_wf_n1 : 0.200 (fn) 0.184 (fp)
     wtreeqmc_wf_n1_shared : 0.259 (fn) 0.244 (fp)
            wtreeqmc_wf_n0 : 0.222 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.246 (fn) 0.246 (fp)
           aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp) 349514.440 (qs) 0.836 (lpp)
                  asteroid : 0.181 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_refined : 0.191 (fn) 0.191 (fp) 349295.820 (qs) 0.840 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
            wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
     wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
            wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  100 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.251 (fn) 0.251 (fp)
           aster_v1.16.3.4 : 0.217 (fn) 0.217 (fp) 5383515.060 (qs) 0.834 (lpp)
                  asteroid : 0.196 (fn) 0.196 (fp)
    wtreeqmc_wf_n2_refined : 0.198 (fn) 0.198 (fp) 5379580.480 (qs) 0.837 (lpp)
            wtreeqmc_wf_n2 : 0.203 (fn) 0.186 (fp)
     wtreeqmc_wf_n2_shared : 0.298 (fn) 0.283 (fp)
            wtreeqmc_wf_n1 : 0.206 (fn) 0.190 (fp)
     wtreeqmc_wf_n1_shared : 0.277 (fn) 0.261 (fp)
            wtreeqmc_wf_n0 : 0.236 (fn) 0.220 (fp)

  125 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.247 (fn) 0.247 (fp)
           aster_v1.16.3.4 : 0.211 (fn) 0.211 (fp) 13586273.440 (qs) 0.839 (lpp)
                  asteroid : 0.185 (fn) 0.185 (fp)
    wtreeqmc_wf_n2_refined : 0.190 (fn) 0.190 (fp) 13575613.660 (qs) 0.842 (lpp)
            wtreeqmc_wf_n2 : 0.194 (fn) 0.176 (fp)
     wtreeqmc_wf_n2_shared : 0.286 (fn) 0.270 (fp)
            wtreeqmc_wf_n1 : 0.199 (fn) 0.180 (fp)
     wtreeqmc_wf_n1_shared : 0.273 (fn) 0.257 (fp)
            wtreeqmc_wf_n0 : 0.226 (fn) 0.210 (fp)

  150 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.256 (fn) 0.256 (fp)
           aster_v1.16.3.4 : 0.217 (fn) 0.217 (fp) 27725670.780 (qs) 0.840 (lpp)
                  asteroid : 0.187 (fn) 0.187 (fp)
    wtreeqmc_wf_n2_refined : 0.192 (fn) 0.192 (fp) 27716694.380 (qs) 0.842 (lpp)
            wtreeqmc_wf_n2 : 0.198 (fn) 0.182 (fp)
     wtreeqmc_wf_n2_shared : 0.303 (fn) 0.289 (fp)
            wtreeqmc_wf_n1 : 0.202 (fn) 0.186 (fp)
     wtreeqmc_wf_n1_shared : 0.287 (fn) 0.272 (fp)
            wtreeqmc_wf_n0 : 0.231 (fn) 0.217 (fp)

Reporting experitment varyngen
  50 taxa, 250 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.471 (fn) 0.471 (fp)
           aster_v1.16.3.4 : 0.431 (fn) 0.431 (fp) 84909.340 (qs) 0.730 (lpp)
                  asteroid : 0.360 (fn) 0.360 (fp)
    wtreeqmc_wf_n2_refined : 0.389 (fn) 0.389 (fp) 84464.160 (qs) 0.734 (lpp)
            wtreeqmc_wf_n2 : 0.412 (fn) 0.364 (fp)
     wtreeqmc_wf_n2_shared : 0.501 (fn) 0.461 (fp)
            wtreeqmc_wf_n1 : 0.413 (fn) 0.365 (fp)
     wtreeqmc_wf_n1_shared : 0.478 (fn) 0.435 (fp)
            wtreeqmc_wf_n0 : 0.466 (fn) 0.424 (fp)

  50 taxa, 500 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.348 (fn) 0.348 (fp)
           aster_v1.16.3.4 : 0.313 (fn) 0.313 (fp) 169684.440 (qs) 0.786 (lpp)
                  asteroid : 0.260 (fn) 0.260 (fp)
    wtreeqmc_wf_n2_refined : 0.275 (fn) 0.275 (fp) 169349.260 (qs) 0.792 (lpp)
            wtreeqmc_wf_n2 : 0.283 (fn) 0.253 (fp)
     wtreeqmc_wf_n2_shared : 0.360 (fn) 0.336 (fp)
            wtreeqmc_wf_n1 : 0.290 (fn) 0.260 (fp)
     wtreeqmc_wf_n1_shared : 0.340 (fn) 0.312 (fp)
            wtreeqmc_wf_n0 : 0.324 (fn) 0.296 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.246 (fn) 0.246 (fp)
           aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp) 349514.440 (qs) 0.836 (lpp)
                  asteroid : 0.181 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_refined : 0.191 (fn) 0.191 (fp) 349295.820 (qs) 0.840 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
            wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
     wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
            wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 2000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.182 (fn) 0.182 (fp)
           aster_v1.16.3.4 : 0.149 (fn) 0.149 (fp) 673901.940 (qs) 0.876 (lpp)
                  asteroid : 0.131 (fn) 0.131 (fp)
    wtreeqmc_wf_n2_refined : 0.127 (fn) 0.127 (fp) 673639.540 (qs) 0.878 (lpp)
            wtreeqmc_wf_n2 : 0.131 (fn) 0.120 (fp)
     wtreeqmc_wf_n2_shared : 0.197 (fn) 0.188 (fp)
            wtreeqmc_wf_n1 : 0.136 (fn) 0.126 (fp)
     wtreeqmc_wf_n1_shared : 0.182 (fn) 0.171 (fp)
            wtreeqmc_wf_n0 : 0.158 (fn) 0.149 (fp)

Reporting experitment varynbps
  50 taxa, 1000 genes, 50 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.292 (fn) 0.292 (fp)
           aster_v1.16.3.4 : 0.295 (fn) 0.295 (fp) 302708.900 (qs) 0.802 (lpp)
                  asteroid : 0.246 (fn) 0.246 (fp)
    wtreeqmc_wf_n2_refined : 0.262 (fn) 0.262 (fp) 302194.800 (qs) 0.806 (lpp)
            wtreeqmc_wf_n2 : 0.263 (fn) 0.250 (fp)
     wtreeqmc_wf_n2_shared : 0.354 (fn) 0.343 (fp)
            wtreeqmc_wf_n1 : 0.268 (fn) 0.254 (fp)
     wtreeqmc_wf_n1_shared : 0.334 (fn) 0.324 (fp)
            wtreeqmc_wf_n0 : 0.301 (fn) 0.291 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.246 (fn) 0.246 (fp)
           aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp) 349514.440 (qs) 0.836 (lpp)
                  asteroid : 0.181 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_refined : 0.191 (fn) 0.191 (fp) 349295.820 (qs) 0.840 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
            wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
     wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
            wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 200 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.231 (fn) 0.231 (fp)
           aster_v1.16.3.4 : 0.171 (fn) 0.171 (fp) 375743.780 (qs) 0.859 (lpp)
                  asteroid : 0.164 (fn) 0.164 (fp)
    wtreeqmc_wf_n2_refined : 0.169 (fn) 0.169 (fp) 375401.740 (qs) 0.865 (lpp)
            wtreeqmc_wf_n2 : 0.175 (fn) 0.153 (fp)
     wtreeqmc_wf_n2_shared : 0.227 (fn) 0.208 (fp)
            wtreeqmc_wf_n1 : 0.176 (fn) 0.154 (fp)
     wtreeqmc_wf_n1_shared : 0.212 (fn) 0.192 (fp)
            wtreeqmc_wf_n0 : 0.193 (fn) 0.173 (fp)

  50 taxa, 1000 genes, 500 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.197 (fn) 0.197 (fp)
           aster_v1.16.3.4 : 0.146 (fn) 0.146 (fp) 364006.260 (qs) 0.879 (lpp)
                  asteroid : 0.132 (fn) 0.132 (fp)
    wtreeqmc_wf_n2_refined : 0.128 (fn) 0.128 (fp) 363838.740 (qs) 0.885 (lpp)
            wtreeqmc_wf_n2 : 0.137 (fn) 0.115 (fp)
     wtreeqmc_wf_n2_shared : 0.176 (fn) 0.156 (fp)
            wtreeqmc_wf_n1 : 0.134 (fn) 0.112 (fp)
     wtreeqmc_wf_n1_shared : 0.163 (fn) 0.144 (fp)
            wtreeqmc_wf_n0 : 0.160 (fn) 0.139 (fp)

Reporting experitment varyblsc
  50 taxa, 1000 genes, 100 bps, 0.050000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.510 (fn) 0.510 (fp)
           aster_v1.16.3.4 : 0.486 (fn) 0.486 (fp) 255664.820 (qs) 0.733 (lpp)
                  asteroid : 0.451 (fn) 0.451 (fp)
    wtreeqmc_wf_n2_refined : 0.466 (fn) 0.466 (fp) 255080.740 (qs) 0.727 (lpp)
            wtreeqmc_wf_n2 : 0.468 (fn) 0.455 (fp)
     wtreeqmc_wf_n2_shared : 0.592 (fn) 0.584 (fp)
            wtreeqmc_wf_n1 : 0.463 (fn) 0.449 (fp)
     wtreeqmc_wf_n1_shared : 0.560 (fn) 0.550 (fp)
            wtreeqmc_wf_n0 : 0.490 (fn) 0.481 (fp)

  50 taxa, 1000 genes, 100 bps, 0.100000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.378 (fn) 0.378 (fp)
           aster_v1.16.3.4 : 0.340 (fn) 0.340 (fp) 298320.600 (qs) 0.779 (lpp)
                  asteroid : 0.293 (fn) 0.293 (fp)
    wtreeqmc_wf_n2_refined : 0.318 (fn) 0.318 (fp) 297751.340 (qs) 0.780 (lpp)
            wtreeqmc_wf_n2 : 0.323 (fn) 0.309 (fp)
     wtreeqmc_wf_n2_shared : 0.434 (fn) 0.424 (fp)
            wtreeqmc_wf_n1 : 0.316 (fn) 0.302 (fp)
     wtreeqmc_wf_n1_shared : 0.414 (fn) 0.402 (fp)
            wtreeqmc_wf_n0 : 0.361 (fn) 0.350 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.246 (fn) 0.246 (fp)
           aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp) 349514.440 (qs) 0.836 (lpp)
                  asteroid : 0.181 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_refined : 0.191 (fn) 0.191 (fp) 349295.820 (qs) 0.840 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
            wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
     wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
            wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 10.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.261 (fn) 0.261 (fp)
           aster_v1.16.3.4 : 0.245 (fn) 0.245 (fp) 311449.060 (qs) 0.825 (lpp)
                  asteroid : 0.207 (fn) 0.207 (fp)
    wtreeqmc_wf_n2_refined : 0.205 (fn) 0.205 (fp) 311084.660 (qs) 0.826 (lpp)
            wtreeqmc_wf_n2 : 0.212 (fn) 0.190 (fp)
     wtreeqmc_wf_n2_shared : 0.284 (fn) 0.266 (fp)
            wtreeqmc_wf_n1 : 0.223 (fn) 0.203 (fp)
     wtreeqmc_wf_n1_shared : 0.280 (fn) 0.263 (fp)
            wtreeqmc_wf_n0 : 0.256 (fn) 0.238 (fp)

  50 taxa, 1000 genes, 100 bps, 100.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.427 (fn) 0.427 (fp)
           aster_v1.16.3.4 : 0.426 (fn) 0.426 (fp) 231697.200 (qs) 0.725 (lpp)
                  asteroid : 0.358 (fn) 0.358 (fp)
    wtreeqmc_wf_n2_refined : 0.368 (fn) 0.368 (fp) 230921.400 (qs) 0.721 (lpp)
            wtreeqmc_wf_n2 : 0.372 (fn) 0.363 (fp)
     wtreeqmc_wf_n2_shared : 0.528 (fn) 0.518 (fp)
            wtreeqmc_wf_n1 : 0.388 (fn) 0.377 (fp)
     wtreeqmc_wf_n1_shared : 0.500 (fn) 0.492 (fp)
            wtreeqmc_wf_n0 : 0.434 (fn) 0.426 (fp)

  50 taxa, 1000 genes, 100 bps, 200.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.501 (fn) 0.501 (fp)
           aster_v1.16.3.4 : 0.495 (fn) 0.495 (fp) 215429.020 (qs) 0.684 (lpp)
                  asteroid : 0.456 (fn) 0.456 (fp)
    wtreeqmc_wf_n2_refined : 0.487 (fn) 0.487 (fp) 214082.720 (qs) 0.684 (lpp)
            wtreeqmc_wf_n2 : 0.493 (fn) 0.482 (fp)
     wtreeqmc_wf_n2_shared : 0.623 (fn) 0.615 (fp)
            wtreeqmc_wf_n1 : 0.495 (fn) 0.485 (fp)
     wtreeqmc_wf_n1_shared : 0.598 (fn) 0.591 (fp)
            wtreeqmc_wf_n0 : 0.533 (fn) 0.525 (fp)

Reporting experitment varymiss
  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.500000 missing param
           wastrid_vanilla : 0.134 (fn) 0.134 (fp)
           aster_v1.16.3.4 : 0.118 (fn) 0.118 (fp) 1456960.240 (qs) 0.895 (lpp)
                  asteroid : 0.099 (fn) 0.099 (fp)
    wtreeqmc_wf_n2_refined : 0.103 (fn) 0.103 (fp) 1456341.840 (qs) 0.897 (lpp)
            wtreeqmc_wf_n2 : 0.105 (fn) 0.103 (fp)
     wtreeqmc_wf_n2_shared : 0.143 (fn) 0.141 (fp)
            wtreeqmc_wf_n1 : 0.107 (fn) 0.105 (fp)
     wtreeqmc_wf_n1_shared : 0.145 (fn) 0.143 (fp)
            wtreeqmc_wf_n0 : 0.123 (fn) 0.121 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.550000 missing param
           wastrid_vanilla : 0.173 (fn) 0.173 (fp)
           aster_v1.16.3.4 : 0.172 (fn) 0.172 (fp) 713047.660 (qs) 0.870 (lpp)
                  asteroid : 0.145 (fn) 0.145 (fp)
    wtreeqmc_wf_n2_refined : 0.143 (fn) 0.143 (fp) 712668.660 (qs) 0.869 (lpp)
            wtreeqmc_wf_n2 : 0.147 (fn) 0.140 (fp)
     wtreeqmc_wf_n2_shared : 0.213 (fn) 0.207 (fp)
            wtreeqmc_wf_n1 : 0.149 (fn) 0.143 (fp)
     wtreeqmc_wf_n1_shared : 0.206 (fn) 0.200 (fp)
            wtreeqmc_wf_n0 : 0.174 (fn) 0.170 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
           wastrid_vanilla : 0.246 (fn) 0.246 (fp)
           aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp) 349514.440 (qs) 0.836 (lpp)
                  asteroid : 0.181 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_refined : 0.191 (fn) 0.191 (fp) 349295.820 (qs) 0.840 (lpp)
            wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
     wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
            wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
     wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
            wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.650000 missing param
           wastrid_vanilla : 0.385 (fn) 0.385 (fp)
           aster_v1.16.3.4 : 0.324 (fn) 0.324 (fp) 147593.780 (qs) 0.785 (lpp)
                  asteroid : 0.280 (fn) 0.280 (fp)
    wtreeqmc_wf_n2_refined : 0.277 (fn) 0.277 (fp) 147322.720 (qs) 0.789 (lpp)
            wtreeqmc_wf_n2 : 0.290 (fn) 0.247 (fp)
     wtreeqmc_wf_n2_shared : 0.388 (fn) 0.356 (fp)
            wtreeqmc_wf_n1 : 0.300 (fn) 0.257 (fp)
     wtreeqmc_wf_n1_shared : 0.372 (fn) 0.338 (fp)
            wtreeqmc_wf_n0 : 0.340 (fn) 0.302 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.700000 missing param
           wastrid_vanilla : 0.607 (fn) 0.607 (fp)
           aster_v1.16.3.4 : 0.498 (fn) 0.498 (fp) 51733.800 (qs) 0.715 (lpp)
                  asteroid : 0.431 (fn) 0.431 (fp)
    wtreeqmc_wf_n2_refined : 0.459 (fn) 0.459 (fp) 51521.200 (qs) 0.714 (lpp)
            wtreeqmc_wf_n2 : 0.488 (fn) 0.415 (fp)
     wtreeqmc_wf_n2_shared : 0.586 (fn) 0.533 (fp)
            wtreeqmc_wf_n1 : 0.491 (fn) 0.419 (fp)
     wtreeqmc_wf_n1_shared : 0.570 (fn) 0.514 (fp)
            wtreeqmc_wf_n0 : 0.519 (fn) 0.455 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.750000 missing param
           wastrid_vanilla : 0.820 (fn) 0.820 (fp)
           aster_v1.16.3.4 : 0.734 (fn) 0.734 (fp) 15659.760 (qs) 0.638 (lpp)
                  asteroid : 0.622 (fn) 0.622 (fp)
    wtreeqmc_wf_n2_refined : 0.662 (fn) 0.662 (fp) 15613.160 (qs) 0.632 (lpp)
            wtreeqmc_wf_n2 : 0.696 (fn) 0.610 (fp)
     wtreeqmc_wf_n2_shared : 0.779 (fn) 0.724 (fp)
            wtreeqmc_wf_n1 : 0.690 (fn) 0.603 (fp)
     wtreeqmc_wf_n1_shared : 0.767 (fn) 0.708 (fp)
            wtreeqmc_wf_n0 : 0.728 (fn) 0.657 (fp)
"""
