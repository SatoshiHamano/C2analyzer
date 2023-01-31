# coding: UTF-8

__author__ = "Satoshi Hamano"

import sys
import os
import numpy as np
import math
from matplotlib.backends.backend_pdf import PdfPages

from spectools_v1_1 import *
from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

if __name__ == "__main__":
    filename = sys.argv[1:]

    config = ReadConfig()
    n_comp, vel_comp = VelocityReader(config)

    rf = open(filename[0], "r")
    rl = rf.readlines()
    rf.close()

    for i in range(len(rl)):
        if rl[i].find("Best fit parameters") != -1:
            paramline = i

    J = []
    v = []
    v_err = []
    b = []
    b_err = []
    logN = []
    logN_err = []
    counter = -1
    for i in range(paramline, len(rl)):
        if rl[i].find("C2J") != -1:
            counter += 1
            J.append(int(rl[i].split()[0].split("J")[1]))
            v.append([])
            v_err.append([])
            b.append([])
            b_err.append([])
            logN.append([])
            logN_err.append([])
            for k in range(n_comp):
                rl1 = rl[i+k+1].split()
                v[counter].append(rl1[2])
                v_err[counter].append(rl1[4])
                b[counter].append(rl1[5])
                b_err[counter].append(rl1[7])
                logN[counter].append(rl1[8])
                logN_err[counter].append(rl1[10])

    config.set("C2 parameters", "min J", str(min(J)))
    config.set("C2 parameters", "max J", str(max(J)))


    for i in range(n_comp):
        config.set("Cloud parameters", "vel%d" % (i+1), v[0][i])
        config.set("Cloud parameters", "vel%d_err" % (i+1), v_err[0][i])
        config.set("Cloud parameters", "b%d" % (i+1), b[0][i])
        config.set("Cloud parameters", "b%d_err" % (i+1), b_err[0][i])

    for i in range(len(J)):
        for j in range(n_comp):
            config.set("Cloud parameters", "logN_J%d_comp%d" % (J[i], j+1), logN[i][j])
            config.set("Cloud parameters", "logN_J%d_comp%d_err" % (J[i], j+1), logN_err[i][j])

    OverwriteConfig(config)

