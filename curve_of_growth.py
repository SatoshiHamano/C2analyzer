# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
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

    config = ReadConfig()
    J, logN, logNerr = ColumnDensityReader(config)
    lines, EW, EWerr = EWReader(config)
    n_comp, vel_comp = VelocityReader(config)

    plt.figure()

    pp = PdfPages(filename[0])

    colors = ["k", "b", "r", "g"]

    for i in range(len(lines)):
        for j in range(len(J)):
            if lines[i].rotJ == J[j]:
                for k in range(n_comp):
                    if logN[j][k] != 0:
                        Nerr = (10**(logN[j][k] + logNerr[j][k]) - 10**(logN[j][k] + logNerr[j][k])) / 2.
                        plt.errorbar(EW[i][k], 10**(logN[j][k]-12) * lines[i].fvalue, xerr=EW[i][k],
                                     yerr=logNerr[j][k] * lines[i].fvalue,
                                     color=colors[k], ecolor=colors[k])
                break

    plt.xlabel("EW ($\AA$)")
    plt.ylabel("$f N$ (cm$^{-2}$)")

    plt.xscale("log")
    plt.yscale("log")

    plt.grid()

    plt.savefig(pp, format="pdf")
    pp.close()
