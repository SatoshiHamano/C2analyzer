# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import numpy as np
from pyraf import iraf
import matplotlib.pyplot as plt

from spectools_v1_1 import openspec
from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]


if __name__ == "__main__":
    filename = sys.argv[1:]

    outputfig = filename[0]

    if len(filename) == 2:
        thres = float(filename[1])
    else:
        thres = 0.7

    config = ReadConfig()
    spfile = config["Spectrum"]["target"]
    telfile = config["Spectrum"]["telluric"]
    resolution = float(config["Spectrum"]["resolution"])
    band = config["C2 parameters"]["band"]
    detection = config["C2 parameters"]["detection"]

    linew = C2params.bandlist[band].minimumlambda() / resolution

    vel_comp = []

    if detection == "nondetection":
        print("C2 band is not detected for this object.")
        sys.exit()
    elif detection == "detection":
        n_comp, vel_comp = VelocityReader(config)
        if config.has_option("Spectrum", "region number"):
            n_region = int(config["Spectrum"]["region number"])
            for i in range(n_region):
                if config.has_option("Spectrum", "region%d range" % (i + 1)):
                    config.remove_option("Spectrum", "region%d range" % (i + 1))
    else:
        detection = alternativequestion("C2 band detection or nondetection?: ", "detection", "nondetection")
        config.set("C2 parameters", "detection", detection)
        if detection == "detection":
            n_comp = int(input("Number of velocity components: "))
            config.set("Cloud parameters", "n_comp", str(n_comp))
            for i in range(n_comp):
                vel = float(input("Velocity of comp. %d: " % (i + 1)))
                config.set("Cloud parameters", "vel%d" % (i + 1), str(vel))
                vel_comp.append(vel)
        else:
            print("C2 band is not detected for this object.")
            sys.exit()

    spx, spy = openspec(spfile)
    telx, tely = openspec(telfile)
    minlam, maxlam = min(spx), max(spx)
    spy_tm = transmittance_resampling(spx, spy, telx, tely)

    C2lines = C2params.bandlist[band].lines
    c2min, c2max = C2params.bandlist[band].minimumlambda(), C2params.bandlist[band].maximumlambda()
    baffer = 20.

    spx_cut = spx[(c2min - baffer < spx) & (spx < baffer + c2max)]
    spy_cut = spy[(c2min - baffer < spx) & (spx < baffer + c2max)]
    spy_tmcut = spy_tm[(c2min - baffer < spx) & (spx < baffer + c2max)]

    regions_start = []
    regions_end = []
    flagtm = False
    for i in range(len(spx_cut)):
        if spy_tmcut[i] > thres:
            if not flagtm:
                regions_start.append(spx_cut[i])
            flagtm = True
        else:
            if flagtm:
                regions_end.append(spx_cut[i])
            flagtm = False
        if i == len(spx_cut) - 1 and flagtm:
            regions_end.append(spx_cut[i])


    united_regions_start = []
    united_regions_end = []
    separation_min = 0.3

    regionid = 0
    while regionid < len(regions_start):
        united_regions_start.append(regions_start[regionid])
        for i in range(regionid, len(regions_start) - 1):
            if regions_start[i + 1] - regions_end[i] < separation_min:
                regionid += 1
            else:
                break
        united_regions_end.append(regions_end[regionid])
        regionid += 1


    fig = plt.figure(figsize=(30, 10))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    spectrum_plot(ax, [spx_cut, spx_cut], [spy_cut, spy_tmcut], [min(spx_cut), max(spx_cut)], [0., 1.4],
                  colors=["cyan", "0.7"], xaxis_label="Wavelength ($\AA$)")
    ax.plot([min(spx_cut), max(spx_cut)], [thres, thres], "r")
    ax.text(min(spx_cut), thres, "  Threshold=%.2f  " % thres, color="r", ha="right", va="center")

    ax.plot([min(spx_cut) + 5., min(spx_cut) + 7.], [1.27, 1.27], "k")
    ax.text(min(spx_cut) + 9., 1.27, "2$\AA$", ha="left", va="center")
    ax.plot([min(spx_cut) + 5., min(spx_cut) + 10.], [1.3, 1.3], "k")
    ax.text(min(spx_cut) + 12., 1.3, "5$\AA$", ha="left", va="center")
    ax.plot([min(spx_cut) + 5., min(spx_cut) + 15.], [1.33, 1.33], "k")
    ax.text(min(spx_cut) + 17., 1.33, "10$\AA$", ha="left", va="center")
    ax.plot([min(spx_cut) + 5., min(spx_cut) + 35.], [1.36, 1.36], "k")
    ax.text(min(spx_cut) + 37., 1.36, "30$\AA$", ha="left", va="center")

    regions_c2include = []
    flagc2 = False

    c2included = [False for i in range(len(C2lines))]

    for i in range(len(united_regions_start)):
        for j in range(len(C2lines)):
            for k in range(n_comp):
                if united_regions_start[i] < C2lines[j].lamair * (1. + vel_comp[k] / (scipy.constants.c * 1.e-3)) < \
                        united_regions_end[i]:
                    flagc2 = True
                    c2included[j] = True
                    break
        if flagc2:
            regions_c2include.append([united_regions_start[i], united_regions_end[i]])
        flagc2 = False

    config.set("Spectrum", "region number", str(len(regions_c2include)))
    for i in range(len(regions_c2include)):
        config.set("Spectrum", "region%d range" % (i + 1),
                   "%.4f,%.4f" % (regions_c2include[i][0], regions_c2include[i][1]))
        ax.plot([regions_c2include[i][0], regions_c2include[i][1]], [1.2, 1.2], "y")
        ax.text(np.average(regions_c2include[i]), 1.22, "Reg. %d" % (i + 1), va="bottom", ha="center")

    for i in range(len(C2lines)):
        if c2included[i]:
            ax.plot([C2lines[i].lamair * (1. + vel_comp[0] / (scipy.constants.c * 1.e-3)),
                     C2lines[i].lamair * (1. + vel_comp[0] / (scipy.constants.c * 1.e-3))
                     ], [1.1, 1.15], color="b")
            config.set("C2 parameters", C2lines[i].transition, "on")
        else:
            ax.plot([C2lines[i].lamair * (1. + vel_comp[0] / (scipy.constants.c * 1.e-3)),
                     C2lines[i].lamair * (1. + vel_comp[0] / (scipy.constants.c * 1.e-3))
                     ], [1.1, 1.15], color="r")
            config.set("C2 parameters", C2lines[i].transition, "off")

    plt.savefig(outputfig)

    OverwriteConfig(config)
