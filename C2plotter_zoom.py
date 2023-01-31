# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from spectools_v1_1 import openspec
from C2_parameters import *
from C2param_parser import *
from spectra_plotter import *
import scipy.constants

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]

if __name__ == "__main__":

    filename = sys.argv[1:]

    outputfig = filename[0]

    config = ReadConfig()
    spfile = config["Spectrum"]["target"]
    band = config["C2 parameters"]["band"]

    spx, spy = openspec(spfile)
    minlam, maxlam = min(spx), max(spx)

    figrow = 5

    figw = 0.25
    figh = 1. / (figrow + 1)
    gap_w = 0.06
    gap_h = 0.06 * 2 / (figrow - 1)
    orig_w = [0.08 + (gap_w + figw) * i for i in range(3)]
    orig_h = [0.05 + (gap_h + figh) * (figrow - 1 - i) for i in range(figrow)]

    velmin = -60
    velmax = 60
    ymin = 0.75
    ymax = 1.1

    C2band = [C2params.bandlist[band].branchlist("R"),
              C2params.bandlist[band].branchlist("Q"),
              C2params.bandlist[band].branchlist("P")]

    linenum = [len(C2band[i]) for i in range(3)]

    pp = PdfPages(outputfig)
    fig = plt.figure(figsize=(12, 16))

    currentrow = 0

    while currentrow < max(linenum):
        drawcount = 0

        for i in range(3):
            for j in C2band[i]:
                if (currentrow * 2 <= j.rotJ < (currentrow + figrow) * 2) and (minlam < j.lamair < maxlam):
                    ax = plt.axes([orig_w[i], orig_h[int(j.rotJ / 2) % figrow], figw, figh])
                    velx = (spx - j.lamair) / j.lamair * scipy.constants.c * 1.e-3
                    spectrum_plot(ax, [velx], [spy], [velmin, velmax], [ymin, ymax],
                                  xticks_val=range(velmin, velmax + 1, 20),
                                  yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1],
                                  xaxis_label="Velocity (km s$^{-1}$)", grid_flag=True)
                    ax.text(velmax, ymin + 0.02, "%s   " % j.transition, ha="right", va="bottom", fontsize=15)
                    drawcount += 1

        currentrow += figrow

        if drawcount > 0:
            plt.savefig(pp, format="pdf")
            plt.clf()

    pp.close()
