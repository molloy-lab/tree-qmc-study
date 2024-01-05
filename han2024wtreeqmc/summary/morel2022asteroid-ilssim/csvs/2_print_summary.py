import pandas
import numpy
import sys

mthds = ["wastrid_vanilla",
         "astral_3_v5.7.8",
         "aster_v1.16.3.4",
         "asteroid",
         "wtreeqmc_wf_n2",
         "wtreeqmc_wf_n2_shared",
         "wtreeqmc_wf_n1",
         "wtreeqmc_wf_n1_shared",
         "wtreeqmc_wf_n0"]

names = ["      wastrid_vanilla",
         "      astral_3_v5.7.8",
         "      aster_v1.16.3.4",
         "             asteroid",
         "       wtreeqmc_wf_n2",
         "wtreeqmc_wf_n2_shared",
         "       wtreeqmc_wf_n1",
         "wtreeqmc_wf_n1_shared",
         "       wtreeqmc_wf_n0"]

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

                            for name, mthd in zip(names, mthds):
                                ydf = xdf[(xdf["MTHD"] == mthd)]

                                data = ydf.SEFN.values
                                if data.size != 50:
                                    sys.stdout.write("ERROR - wrong number of replicates!\n")

                                nrepl = numpy.sum(~numpy.isnan(data))

                                fnavg = numpy.nanmean(ydf.SEFNR.values)
                                fnavg = numpy.round(numpy.round(fnavg, 4), 3)

                                fpavg = numpy.nanmean(ydf.SEFPR.values)
                                fpavg = numpy.round(numpy.round(fpavg, 4), 3)

                                print("    %s : %1.3f (fn) %1.3f (fp)" % (name, fnavg, fpavg))

                            print()

"""
Reporting experitment varypsiz
  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 10 pop size, 0.600000 missing param
          wastrid_vanilla : 0.212 (fn) 0.212 (fp)
          astral_3_v5.7.8 : 0.195 (fn) 0.195 (fp)
          aster_v1.16.3.4 : 0.177 (fn) 0.177 (fp)
                 asteroid : 0.144 (fn) 0.144 (fp)
           wtreeqmc_wf_n2 : 0.158 (fn) 0.144 (fp)
    wtreeqmc_wf_n2_shared : 0.226 (fn) 0.215 (fp)
           wtreeqmc_wf_n1 : 0.161 (fn) 0.146 (fp)
    wtreeqmc_wf_n1_shared : 0.215 (fn) 0.202 (fp)
           wtreeqmc_wf_n0 : 0.189 (fn) 0.176 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.246 (fn) 0.246 (fp)
          astral_3_v5.7.8 : 0.220 (fn) 0.220 (fp)
          aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp)
                 asteroid : 0.181 (fn) 0.181 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
           wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
    wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
           wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 100000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.297 (fn) 0.297 (fp)
          astral_3_v5.7.8 : 0.267 (fn) 0.267 (fp)
          aster_v1.16.3.4 : 0.262 (fn) 0.262 (fp)
                 asteroid : 0.232 (fn) 0.232 (fp)
           wtreeqmc_wf_n2 : 0.240 (fn) 0.224 (fp)
    wtreeqmc_wf_n2_shared : 0.334 (fn) 0.321 (fp)
           wtreeqmc_wf_n1 : 0.243 (fn) 0.227 (fp)
    wtreeqmc_wf_n1_shared : 0.310 (fn) 0.297 (fp)
           wtreeqmc_wf_n0 : 0.274 (fn) 0.259 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 500000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.475 (fn) 0.475 (fp)
          astral_3_v5.7.8 : 0.469 (fn) 0.469 (fp)
          aster_v1.16.3.4 : 0.462 (fn) 0.462 (fp)
                 asteroid : 0.434 (fn) 0.434 (fp)
           wtreeqmc_wf_n2 : 0.424 (fn) 0.416 (fp)
    wtreeqmc_wf_n2_shared : 0.524 (fn) 0.518 (fp)
           wtreeqmc_wf_n1 : 0.430 (fn) 0.421 (fp)
    wtreeqmc_wf_n1_shared : 0.516 (fn) 0.509 (fp)
           wtreeqmc_wf_n0 : 0.461 (fn) 0.452 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 1000000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.585 (fn) 0.585 (fp)
          astral_3_v5.7.8 : 0.600 (fn) 0.600 (fp)
          aster_v1.16.3.4 : 0.588 (fn) 0.588 (fp)
                 asteroid : 0.553 (fn) 0.553 (fp)
           wtreeqmc_wf_n2 : 0.561 (fn) 0.553 (fp)
    wtreeqmc_wf_n2_shared : 0.669 (fn) 0.662 (fp)
           wtreeqmc_wf_n1 : 0.562 (fn) 0.554 (fp)
    wtreeqmc_wf_n1_shared : 0.647 (fn) 0.640 (fp)
           wtreeqmc_wf_n0 : 0.610 (fn) 0.603 (fp)

Reporting experitment varyntax
  25 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.269 (fn) 0.269 (fp)
          astral_3_v5.7.8 : 0.270 (fn) 0.270 (fp)
          aster_v1.16.3.4 : 0.276 (fn) 0.276 (fp)
                 asteroid : 0.226 (fn) 0.226 (fp)
           wtreeqmc_wf_n2 : 0.242 (fn) 0.223 (fp)
    wtreeqmc_wf_n2_shared : 0.288 (fn) 0.270 (fp)
           wtreeqmc_wf_n1 : 0.240 (fn) 0.222 (fp)
    wtreeqmc_wf_n1_shared : 0.284 (fn) 0.267 (fp)
           wtreeqmc_wf_n0 : 0.281 (fn) 0.265 (fp)

  75 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.244 (fn) 0.244 (fp)
          astral_3_v5.7.8 : 0.215 (fn) 0.215 (fp)
          aster_v1.16.3.4 : 0.209 (fn) 0.209 (fp)
                 asteroid : 0.183 (fn) 0.183 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.284 (fn) 0.269 (fp)
           wtreeqmc_wf_n1 : 0.200 (fn) 0.184 (fp)
    wtreeqmc_wf_n1_shared : 0.259 (fn) 0.244 (fp)
           wtreeqmc_wf_n0 : 0.222 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.246 (fn) 0.246 (fp)
          astral_3_v5.7.8 : 0.220 (fn) 0.220 (fp)
          aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp)
                 asteroid : 0.181 (fn) 0.181 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
           wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
    wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
           wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  100 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.251 (fn) 0.251 (fp)
          astral_3_v5.7.8 : 0.226 (fn) 0.226 (fp)
          aster_v1.16.3.4 : 0.217 (fn) 0.217 (fp)
                 asteroid : 0.196 (fn) 0.196 (fp)
           wtreeqmc_wf_n2 : 0.203 (fn) 0.186 (fp)
    wtreeqmc_wf_n2_shared : 0.298 (fn) 0.283 (fp)
           wtreeqmc_wf_n1 : 0.206 (fn) 0.190 (fp)
    wtreeqmc_wf_n1_shared : 0.277 (fn) 0.261 (fp)
           wtreeqmc_wf_n0 : 0.236 (fn) 0.220 (fp)

  125 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.247 (fn) 0.247 (fp)
          astral_3_v5.7.8 : 0.223 (fn) 0.223 (fp)
          aster_v1.16.3.4 : 0.211 (fn) 0.211 (fp)
                 asteroid : 0.185 (fn) 0.185 (fp)
           wtreeqmc_wf_n2 : 0.194 (fn) 0.176 (fp)
    wtreeqmc_wf_n2_shared : 0.286 (fn) 0.270 (fp)
           wtreeqmc_wf_n1 : 0.199 (fn) 0.180 (fp)
    wtreeqmc_wf_n1_shared : 0.273 (fn) 0.257 (fp)
           wtreeqmc_wf_n0 : 0.226 (fn) 0.210 (fp)

  150 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.256 (fn) 0.256 (fp)
          astral_3_v5.7.8 : 0.229 (fn) 0.229 (fp)
          aster_v1.16.3.4 : 0.217 (fn) 0.217 (fp)
                 asteroid : 0.187 (fn) 0.187 (fp)
           wtreeqmc_wf_n2 : 0.198 (fn) 0.182 (fp)
    wtreeqmc_wf_n2_shared : 0.303 (fn) 0.289 (fp)
           wtreeqmc_wf_n1 : 0.202 (fn) 0.186 (fp)
    wtreeqmc_wf_n1_shared : 0.287 (fn) 0.272 (fp)
           wtreeqmc_wf_n0 : 0.231 (fn) 0.217 (fp)

Reporting experitment varyngen
  50 taxa, 250 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.471 (fn) 0.471 (fp)
          astral_3_v5.7.8 : 0.448 (fn) 0.448 (fp)
          aster_v1.16.3.4 : 0.431 (fn) 0.431 (fp)
                 asteroid : 0.360 (fn) 0.360 (fp)
           wtreeqmc_wf_n2 : 0.412 (fn) 0.364 (fp)
    wtreeqmc_wf_n2_shared : 0.501 (fn) 0.461 (fp)
           wtreeqmc_wf_n1 : 0.413 (fn) 0.365 (fp)
    wtreeqmc_wf_n1_shared : 0.478 (fn) 0.435 (fp)
           wtreeqmc_wf_n0 : 0.466 (fn) 0.424 (fp)

  50 taxa, 500 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.348 (fn) 0.348 (fp)
          astral_3_v5.7.8 : 0.330 (fn) 0.330 (fp)
          aster_v1.16.3.4 : 0.313 (fn) 0.313 (fp)
                 asteroid : 0.260 (fn) 0.260 (fp)
           wtreeqmc_wf_n2 : 0.283 (fn) 0.253 (fp)
    wtreeqmc_wf_n2_shared : 0.360 (fn) 0.336 (fp)
           wtreeqmc_wf_n1 : 0.290 (fn) 0.260 (fp)
    wtreeqmc_wf_n1_shared : 0.340 (fn) 0.312 (fp)
           wtreeqmc_wf_n0 : 0.324 (fn) 0.296 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.246 (fn) 0.246 (fp)
          astral_3_v5.7.8 : 0.220 (fn) 0.220 (fp)
          aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp)
                 asteroid : 0.181 (fn) 0.181 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
           wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
    wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
           wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 2000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.182 (fn) 0.182 (fp)
          astral_3_v5.7.8 : 0.148 (fn) 0.148 (fp)
          aster_v1.16.3.4 : 0.149 (fn) 0.149 (fp)
                 asteroid : 0.131 (fn) 0.131 (fp)
           wtreeqmc_wf_n2 : 0.131 (fn) 0.120 (fp)
    wtreeqmc_wf_n2_shared : 0.197 (fn) 0.188 (fp)
           wtreeqmc_wf_n1 : 0.136 (fn) 0.126 (fp)
    wtreeqmc_wf_n1_shared : 0.182 (fn) 0.171 (fp)
           wtreeqmc_wf_n0 : 0.158 (fn) 0.149 (fp)

Reporting experitment varynbps
  50 taxa, 1000 genes, 50 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.292 (fn) 0.292 (fp)
          astral_3_v5.7.8 : 0.291 (fn) 0.291 (fp)
          aster_v1.16.3.4 : 0.295 (fn) 0.295 (fp)
                 asteroid : 0.246 (fn) 0.246 (fp)
           wtreeqmc_wf_n2 : 0.263 (fn) 0.250 (fp)
    wtreeqmc_wf_n2_shared : 0.354 (fn) 0.343 (fp)
           wtreeqmc_wf_n1 : 0.268 (fn) 0.254 (fp)
    wtreeqmc_wf_n1_shared : 0.334 (fn) 0.324 (fp)
           wtreeqmc_wf_n0 : 0.301 (fn) 0.291 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.246 (fn) 0.246 (fp)
          astral_3_v5.7.8 : 0.220 (fn) 0.220 (fp)
          aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp)
                 asteroid : 0.181 (fn) 0.181 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
           wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
    wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
           wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 200 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.231 (fn) 0.231 (fp)
          astral_3_v5.7.8 : 0.185 (fn) 0.185 (fp)
          aster_v1.16.3.4 : 0.171 (fn) 0.171 (fp)
                 asteroid : 0.164 (fn) 0.164 (fp)
           wtreeqmc_wf_n2 : 0.175 (fn) 0.153 (fp)
    wtreeqmc_wf_n2_shared : 0.227 (fn) 0.208 (fp)
           wtreeqmc_wf_n1 : 0.176 (fn) 0.154 (fp)
    wtreeqmc_wf_n1_shared : 0.212 (fn) 0.192 (fp)
           wtreeqmc_wf_n0 : 0.193 (fn) 0.173 (fp)

  50 taxa, 1000 genes, 500 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.197 (fn) 0.197 (fp)
          astral_3_v5.7.8 : 0.156 (fn) 0.156 (fp)
          aster_v1.16.3.4 : 0.146 (fn) 0.146 (fp)
                 asteroid : 0.132 (fn) 0.132 (fp)
           wtreeqmc_wf_n2 : 0.137 (fn) 0.115 (fp)
    wtreeqmc_wf_n2_shared : 0.176 (fn) 0.156 (fp)
           wtreeqmc_wf_n1 : 0.134 (fn) 0.112 (fp)
    wtreeqmc_wf_n1_shared : 0.163 (fn) 0.144 (fp)
           wtreeqmc_wf_n0 : 0.160 (fn) 0.139 (fp)

Reporting experitment varyblsc
  50 taxa, 1000 genes, 100 bps, 0.050000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.510 (fn) 0.510 (fp)
          astral_3_v5.7.8 : 0.512 (fn) 0.512 (fp)
          aster_v1.16.3.4 : 0.486 (fn) 0.486 (fp)
                 asteroid : 0.451 (fn) 0.451 (fp)
           wtreeqmc_wf_n2 : 0.468 (fn) 0.455 (fp)
    wtreeqmc_wf_n2_shared : 0.592 (fn) 0.584 (fp)
           wtreeqmc_wf_n1 : 0.463 (fn) 0.449 (fp)
    wtreeqmc_wf_n1_shared : 0.560 (fn) 0.550 (fp)
           wtreeqmc_wf_n0 : 0.490 (fn) 0.481 (fp)

  50 taxa, 1000 genes, 100 bps, 0.100000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.378 (fn) 0.378 (fp)
          astral_3_v5.7.8 : 0.348 (fn) 0.348 (fp)
          aster_v1.16.3.4 : 0.340 (fn) 0.340 (fp)
                 asteroid : 0.293 (fn) 0.293 (fp)
           wtreeqmc_wf_n2 : 0.323 (fn) 0.309 (fp)
    wtreeqmc_wf_n2_shared : 0.434 (fn) 0.424 (fp)
           wtreeqmc_wf_n1 : 0.316 (fn) 0.302 (fp)
    wtreeqmc_wf_n1_shared : 0.414 (fn) 0.402 (fp)
           wtreeqmc_wf_n0 : 0.361 (fn) 0.350 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.246 (fn) 0.246 (fp)
          astral_3_v5.7.8 : 0.220 (fn) 0.220 (fp)
          aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp)
                 asteroid : 0.181 (fn) 0.181 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
           wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
    wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
           wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 10.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.261 (fn) 0.261 (fp)
          astral_3_v5.7.8 : 0.269 (fn) 0.269 (fp)
          aster_v1.16.3.4 : 0.245 (fn) 0.245 (fp)
                 asteroid : 0.207 (fn) 0.207 (fp)
           wtreeqmc_wf_n2 : 0.212 (fn) 0.190 (fp)
    wtreeqmc_wf_n2_shared : 0.284 (fn) 0.266 (fp)
           wtreeqmc_wf_n1 : 0.223 (fn) 0.203 (fp)
    wtreeqmc_wf_n1_shared : 0.280 (fn) 0.263 (fp)
           wtreeqmc_wf_n0 : 0.256 (fn) 0.238 (fp)

  50 taxa, 1000 genes, 100 bps, 100.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.427 (fn) 0.427 (fp)
          astral_3_v5.7.8 : 0.447 (fn) 0.447 (fp)
          aster_v1.16.3.4 : 0.426 (fn) 0.426 (fp)
                 asteroid : 0.358 (fn) 0.358 (fp)
           wtreeqmc_wf_n2 : 0.372 (fn) 0.363 (fp)
    wtreeqmc_wf_n2_shared : 0.528 (fn) 0.518 (fp)
           wtreeqmc_wf_n1 : 0.388 (fn) 0.377 (fp)
    wtreeqmc_wf_n1_shared : 0.500 (fn) 0.492 (fp)
           wtreeqmc_wf_n0 : 0.434 (fn) 0.426 (fp)

  50 taxa, 1000 genes, 100 bps, 200.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.501 (fn) 0.501 (fp)
          astral_3_v5.7.8 : 0.518 (fn) 0.518 (fp)
          aster_v1.16.3.4 : 0.495 (fn) 0.495 (fp)
                 asteroid : 0.456 (fn) 0.456 (fp)
           wtreeqmc_wf_n2 : 0.493 (fn) 0.482 (fp)
    wtreeqmc_wf_n2_shared : 0.623 (fn) 0.615 (fp)
           wtreeqmc_wf_n1 : 0.495 (fn) 0.485 (fp)
    wtreeqmc_wf_n1_shared : 0.598 (fn) 0.591 (fp)
           wtreeqmc_wf_n0 : 0.533 (fn) 0.525 (fp)

Reporting experitment varymiss
  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.500000 missing param
          wastrid_vanilla : 0.134 (fn) 0.134 (fp)
          astral_3_v5.7.8 : 0.117 (fn) 0.117 (fp)
          aster_v1.16.3.4 : 0.118 (fn) 0.118 (fp)
                 asteroid : 0.099 (fn) 0.099 (fp)
           wtreeqmc_wf_n2 : 0.105 (fn) 0.103 (fp)
    wtreeqmc_wf_n2_shared : 0.143 (fn) 0.141 (fp)
           wtreeqmc_wf_n1 : 0.107 (fn) 0.105 (fp)
    wtreeqmc_wf_n1_shared : 0.145 (fn) 0.143 (fp)
           wtreeqmc_wf_n0 : 0.123 (fn) 0.121 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.550000 missing param
          wastrid_vanilla : 0.173 (fn) 0.173 (fp)
          astral_3_v5.7.8 : 0.169 (fn) 0.169 (fp)
          aster_v1.16.3.4 : 0.172 (fn) 0.172 (fp)
                 asteroid : 0.145 (fn) 0.145 (fp)
           wtreeqmc_wf_n2 : 0.147 (fn) 0.140 (fp)
    wtreeqmc_wf_n2_shared : 0.213 (fn) 0.207 (fp)
           wtreeqmc_wf_n1 : 0.149 (fn) 0.143 (fp)
    wtreeqmc_wf_n1_shared : 0.206 (fn) 0.200 (fp)
           wtreeqmc_wf_n0 : 0.174 (fn) 0.170 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.600000 missing param
          wastrid_vanilla : 0.246 (fn) 0.246 (fp)
          astral_3_v5.7.8 : 0.220 (fn) 0.220 (fp)
          aster_v1.16.3.4 : 0.204 (fn) 0.204 (fp)
                 asteroid : 0.181 (fn) 0.181 (fp)
           wtreeqmc_wf_n2 : 0.196 (fn) 0.181 (fp)
    wtreeqmc_wf_n2_shared : 0.258 (fn) 0.246 (fp)
           wtreeqmc_wf_n1 : 0.196 (fn) 0.182 (fp)
    wtreeqmc_wf_n1_shared : 0.252 (fn) 0.237 (fp)
           wtreeqmc_wf_n0 : 0.220 (fn) 0.207 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.650000 missing param
          wastrid_vanilla : 0.385 (fn) 0.385 (fp)
          astral_3_v5.7.8 : 0.342 (fn) 0.342 (fp)
          aster_v1.16.3.4 : 0.324 (fn) 0.324 (fp)
                 asteroid : 0.280 (fn) 0.280 (fp)
           wtreeqmc_wf_n2 : 0.290 (fn) 0.247 (fp)
    wtreeqmc_wf_n2_shared : 0.388 (fn) 0.356 (fp)
           wtreeqmc_wf_n1 : 0.300 (fn) 0.257 (fp)
    wtreeqmc_wf_n1_shared : 0.372 (fn) 0.338 (fp)
           wtreeqmc_wf_n0 : 0.340 (fn) 0.302 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.700000 missing param
          wastrid_vanilla : 0.607 (fn) 0.607 (fp)
          astral_3_v5.7.8 : 0.558 (fn) 0.558 (fp)
          aster_v1.16.3.4 : 0.498 (fn) 0.498 (fp)
                 asteroid : 0.431 (fn) 0.431 (fp)
           wtreeqmc_wf_n2 : 0.488 (fn) 0.415 (fp)
    wtreeqmc_wf_n2_shared : 0.586 (fn) 0.533 (fp)
           wtreeqmc_wf_n1 : 0.491 (fn) 0.419 (fp)
    wtreeqmc_wf_n1_shared : 0.570 (fn) 0.514 (fp)
           wtreeqmc_wf_n0 : 0.519 (fn) 0.455 (fp)

  50 taxa, 1000 genes, 100 bps, 1.000000 bl scaler, 50000000 pop size, 0.750000 missing param
          wastrid_vanilla : 0.820 (fn) 0.820 (fp)
          astral_3_v5.7.8 : 0.767 (fn) 0.767 (fp)
          aster_v1.16.3.4 : 0.734 (fn) 0.734 (fp)
                 asteroid : 0.622 (fn) 0.622 (fp)
           wtreeqmc_wf_n2 : 0.696 (fn) 0.610 (fp)
    wtreeqmc_wf_n2_shared : 0.779 (fn) 0.724 (fp)
           wtreeqmc_wf_n1 : 0.690 (fn) 0.603 (fp)
    wtreeqmc_wf_n1_shared : 0.767 (fn) 0.708 (fp)
           wtreeqmc_wf_n0 : 0.728 (fn) 0.657 (fp)
"""
