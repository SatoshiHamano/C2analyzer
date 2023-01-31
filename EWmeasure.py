# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import numpy as np
from pyraf import iraf
import datetime

from spectools_v1_1 import openspec
from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]


def gaussianfunc(x, peak, l, w):
    y = peak * numpy.exp(-(x - l) ** 2 / (2 * w ** 2))
    return y


def multigaussian(vel, res, x, comp, p, c, wn):
    y = [1. for i in range(len(x))]

    for i in range(len(x)):
        tmpy = 1.
        numpeak = 0
        for ii in range(comp):
            for j in range(len(wn)):
                l = wn[j] * (1.0 + (vel[ii]) / c)
                w = l / (c / res[ii])
                tmpy *= 1. - p[numpeak] * math.exp(-(x[i] - l) ** 2 / (2 * w ** 2))
                numpeak += 1

        y[i] = tmpy

    return numpy.array(y)


def residue1(p, y, x, vel, res, comp, c, wn, fluxerror):
    res = (y - multigaussian(vel, res, x, comp, p, c, wn)) / fluxerror
    return (res)


if __name__ == "__main__":
    filename = sys.argv[1:]

    config = ReadConfig()
    band = config["C2 parameters"]["band"]
    resolution = float(config["Spectrum"]["resolution"])
    light_sp = scipy.constants.c * 1.e-3
    rp = math.sqrt(math.pi)
    vres = light_sp / resolution

    n_comp, vel_comp, _, b_comp, _ = CloundParameterReader(config)
    n_region, region_ranges, target_region, telluric_region, target_region_normfits, target_region_normtxt, region_snr = RegionReader(
        config)

    C2lines = C2params.bandlist[band].lines

    res = [math.sqrt((b_comp[i] / 2 ** 0.5) ** 2. + (vres / (2 * (2 * math.log(2)) ** 0.5)) ** 2.) for i in
           range(n_comp)]

    rmsfactor = [-0.5, 0., 0.5]

    output = open(filename[0], "a")

    output.write("### Time: %s ###\n\n" % (datetime.datetime.today()))

    output.write("# initial input parameters for C2 fitting\n\n")
    output.write("\tC2 band (1-0 or 0-0):\t%s\n" % band)
    output.write("\tNumber of velocity components:\t%d\n" % n_comp)
    for i in range(n_comp):
        output.write("\tVelocity of component %d:\t%.2f\n" % (i + 1, vel_comp[i]))
        output.write("\tDopper width of component %d:\t%.2f\n" % (i + 1, b_comp[i]))

    output.write("\n\n")

    for i in range(n_region):
        spx, spy, fluxerror, maskregion = numpy.loadtxt(target_region_normtxt[i],
                                                        dtype=[("col1", "f8"), ("col2", "f8"), ("col3", "f8"),
                                                               ("col4", "int16")], unpack=True)
        minlam, maxlam = np.amin(spx), np.amax(spx)

        lineinclude = []
        linelambda = []
        for line in C2lines:
            if minlam < line.lamair < maxlam:
                lineinclude.append(line)
                linelambda.append(line.lamair)
        n_line = len(linelambda)

        output.write("Input spectrum file: %s\n" % target_region_normtxt[i])
        output.write("Wavelength range in vacuum: %.1f - %.1f\n\n" % (minlam, maxlam))
        output.write("Fitted C2 bands (transition and wavelength in vacuum)\n\n")
        for line in lineinclude:
            output.write("\t%s\t%.3f\n" % (line.transition, line.lamair))

        output.write("\n\n")

        ewrms = []
        ewerrrms = []

        for j in range(len(rmsfactor)):
            spy_shift = spy + 1. / region_snr[i] * rmsfactor[j]

            param_output1 = scipy.optimize.leastsq(residue1, p0, args=(
                spy_shift, spx, vel_comp, res, n_comp, light_sp, linelambda, fluxerror),
                                                   full_output=True)
            fit_peaks = [param_output1[0][ii] for ii in range(n_comp * n_line)]
            fit_errors = [param_output1[1][ii][ii] ** 0.5 for ii in range(n_comp * n_line)]
            fit_ews = [
                [fit_peaks[jj + ii * n_line] * linelambda[jj] * res[ii] / light_sp * rp * math.sqrt(2) * 1e+3
                 for ii in range(n_comp)] for jj in range(n_line)]
            fit_ews_errors = [
                [fit_errors[jj + ii * n_line] * linelambda[jj] * res[ii] / light_sp * rp * math.sqrt(2) * 1e+3
                 for ii in range(n_comp)] for jj in range(n_line)]

            ewrms.append(fit_ews)
            ewerrrms.append(fit_errors)

        for j in range(n_comp):
            for k in range(n_line):
                error_stat = ewerrrms[1][k][j]
                error_cont = math.fabs(ewrms[0][k][j] - ewrms[2][k][j]) / 2.
                error_total = (error_cont ** 2. + error_stat ** 2.) ** 0.5
                output.write(
                    "%s Comp. %d EW = %.2f pm %.2f\n" % (lineinclude[k].transition, j + 1, ewrms[1][k][j], error_total))
                config.set("Cloud parameters", "ew_%s_comp%d" % (lineinclude[k].transition, j+1), "%.2f" % ewrms[1][k][j])
                config.set("Cloud parameters", "ew_%s_comp%d_err" % (lineinclude[k].rotJ, j+1), "%.2f" % error_total)

    output.close()
    OverwriteConfig(config)