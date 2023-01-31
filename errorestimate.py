# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import numpy as np
from pyraf import iraf

from spectools_v1_1 import openspecfits
from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]


def savespecascii(objfits, telfits, output, snr, band, vel_comp, resolution):
    spxtel, spytel, _, _, _ = openspecfits(telfits)
    spx, spy, _, _, _ = openspecfits(objfits)
    C2lines = C2params.bandlist[band].lines
    spy_tm = transmittance_resampling(spx, spy, spxtel, spytel)

    thres = 0.6
    thres_snr = 0.9

    lightvel = scipy.constants.c * 1.e-3
    vel_low = min(vel_comp) - lightvel / resolution * 2.
    vel_high = max(vel_comp) + lightvel / resolution * 2.

    maskregion = np.ones(spx.size)
    maskregion[spy_tm < thres] = 0

    maskregion_snr = np.array([True for i in range(spx.size)])
    maskregion_snr[spy_tm < thres_snr] = False
    maskregion_snr[spy == 0.0] = False
    maskregion_snr[spy > 1.01] = False
    maskregion_snr[spy < 0.99] = False

    for line in C2lines:
        w_low = line.lamair * (1.0 + vel_low / lightvel)
        w_high = line.lamair * (1.0 * vel_high / lightvel)
        maskregion_snr[(w_low < spxtel) & (w_high > spxtel)] = False

    if snr == "INDEF":
        snr = 1. / (np.std(spy[maskregion_snr]))
    errarray = np.ones(spx.size) / snr / (spy_tm ** 0.5)

    deleteid = []
    for i in range(spx.size):
        if spy[i] == 0.0:
            deleteid.append(i)
    for i in range(len(deleteid)):
        spx = np.delete(spx, deleteid[len(deleteid) - 1 - i])
        spy = np.delete(spy, deleteid[len(deleteid) - 1 - i])
        errarray = np.delete(errarray, deleteid[len(deleteid) - 1 - i])
        maskregion = np.delete(maskregion, deleteid[len(deleteid) - 1 - i])

    np.savetxt(output, np.vstack((spx, spy, errarray, maskregion)).T, fmt="%.10e",
                  header="Object: %s\nTelluric: %s\nThreshold of mask: %.2f\nSignal-to-noise: %d" % (
                      objfits, telfits, thres, snr))

    return snr


if __name__ == "__main__":
    config = ReadConfig()
    spfile = config["Spectrum"]["target"]
    telfile = config["Spectrum"]["telluric"]
    resolution = float(config["Spectrum"]["resolution"])
    band = config["C2 parameters"]["band"]


    n_region, region_ranges, target_region, telluric_region, target_region_normfits, _, region_snr = RegionReader(
        config)
    n_comp, vel_comp = VelocityReader(config)

    target_region_normtxt = [target_region_normfits[i].replace(".fits", ".txt") for i in range(n_region)]

    if region_snr == []:
        region_snr = ["INDEF" for i in range(n_region)]

    for i in range(n_region):
        snr = savespecascii(target_region_normfits[i], telluric_region[i], target_region_normtxt[i], region_snr[i], band,
                      vel_comp, resolution)
        if region_snr[i] == "INDEF":
            config.set("Spectrum", "region%d snr" % (i + 1), str(snr))

    OverwriteConfig(config)
