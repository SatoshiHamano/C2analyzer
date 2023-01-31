# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import math
import subprocess
import numpy as np
from pyraf import iraf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from spectools_v1_1 import openspec
from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]

if __name__ == "__main__":
    filename = sys.argv[1:]

    outputfig = filename[0]

    config = ReadConfig()
    spfile = config["Spectrum"]["target"]
    telfile = config["Spectrum"]["telluric"]
    resolution = float(config["Spectrum"]["resolution"])
    band = config["C2 parameters"]["band"]

    n_comp, vel_comp = VelocityReader(config)
    n_region, region_ranges, target_region, telluric_region, _, _, _ = RegionReader(config)

    thres = 0.7
    light_sp = scipy.constants.c * 1.e-3
    vrange = light_sp / resolution * 2.
    C2lines = C2params.bandlist[band].lines

    pp = PdfPages(outputfig)
    plt.figure(figsize=(20,6))

    for i in range(n_region):
        targetsp = spfile.rstrip("fits").rstrip(".") + "_reg%d.fits" % (i + 1)
        telsp = telfile.rstrip("fits").rstrip(".") + "_reg%d.fits" % (i + 1)
        targetsp_norm = spfile.rstrip("fits").rstrip(".") + "_reg%d_norm.fits" % (i + 1)
        iraf.scopy(spfile, targetsp, w1=region_ranges[i][0], w2=region_ranges[i][1])
        iraf.scopy(telfile, telsp, w1=region_ranges[i][0], w2=region_ranges[i][1])

        spx, spy = openspec(targetsp)
        telx, tely = openspec(telsp)
        spy_tm = transmittance_resampling(spx, spy, telx, tely)

        plt.step(spx, spy, where="mid", color="k")
        plt.step(spx, spy_tm, where="mid", color="0.7")
        plt.plot([min(spx), max(spx)], [1.15, 1.15], "g")

        prewave = 0.
        normalize_range = ""
        for j in range(1, len(spx)):
            flag = 0
            if spy_tm[j] > thres:
                for line in C2lines:
                    for k in range(n_comp):
                        if math.fabs((spx[j] - line.lamair) / line.lamair * light_sp - vel_comp[k]) < vrange:
                            # if i == 1:
                            #     print(flag, line.transition, spx[j], (spx[j] - line.lamair) / line.lamair * light_sp)
                            flag += 1
            if flag == 0:
                if prewave != spx[j - 1]:
                    if prewave != 0.:
                        normalize_range += ":%.1f" % prewave
                        endw = spx[j]
                        plt.plot([startw, endw], [1.2, 1.2], "b")
                    normalize_range += ",%.1f" % spx[j]
                    startw = spx[j]
                elif j == len(spx) - 1:
                    normalize_range += ":%.1f" % spx[j]
                    endw = spx[j]
                    plt.plot([startw, endw], [1.2, 1.2], "b")
                prewave = spx[j]

        subprocess.call("python continuum.py %s %s %s" % (targetsp, targetsp_norm, normalize_range.lstrip(",")), shell=True)

        config.set("Spectrum", "region%d target" % (i + 1), targetsp)
        config.set("Spectrum", "region%d target normfits" % (i + 1), targetsp_norm)
        config.set("Spectrum", "region%d telluric" % (i + 1), telsp)
        config.set("Spectrum", "region%d snr" % (i + 1), "")


    for line in C2lines:
        plt.scatter(line.lamair, 1.1, marker="+")

    plt.ylim(0., 1.3)
    plt.savefig(pp, format="pdf")
    OverwriteConfig(config)

    pp.close()
