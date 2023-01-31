# coding: UTF-8

__author__ = "Satoshi Hamano"

import sys
import VoigtFit
import numpy as np
import pickle
from scipy.signal import fftconvolve, gaussian
import VoigtFit.voigt as voigt
from VoigtFit.dataset import Line
import matplotlib.pyplot as plt
import math
from matplotlib.backends.backend_pdf import PdfPages

from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]

if __name__ == "__main__":
    filename = sys.argv[1:]

    config = ReadConfig()
    spfile = config["Spectrum"]["target"]
    telfile = config["Spectrum"]["telluric"]
    resolution = float(config["Spectrum"]["resolution"])
    band = config["C2 parameters"]["band"]

    n_comp, vel_comp = VelocityReader(config)
    n_region, region_ranges, target_region, telluric_region, target_region_normfits, target_region_normtxt, region_snr = RegionReader(
        config)

    C2lines = C2params.bandlist[band].lines
    C2bool = [config.getboolean("C2 parameters", c2line.transition) for c2line in C2lines]
    C2available = []
    C2fitlevel = []
    for i in range(len(C2lines)):
        if C2bool[i]:
            C2available.append(C2lines[i])
            C2fitlevel.append(C2lines[i].level)

    light_sp = scipy.constants.c * 1.e-3
    vres = light_sp / resolution
    sqtwopi = math.sqrt(2. * math.pi)
    gfwhm_factor = 2 * math.sqrt(2. * math.log(2))
    velspan = 150

    VFdataset = VoigtFit.DataSet(0.0)
    for i in range(n_region):
        spx, spy, fluxerror, maskregion = \
            np.loadtxt(target_region_normtxt[i],
                       dtype=[("col1", "f8"), ("col2", "f8"), ("col3", "f8"), ("col4", "int16")], unpack=True)
        spx = airtovac(spx)
        VFdataset.add_data(spx, spy, vres, err=fluxerror, normalized=True, mask=maskregion)

    VFdataset.cheb_order = -1
    for c2line in C2available:
        VFdataset.add_line(c2line.linetag, velspan=velspan)

    VFdataset.prepare_dataset(norm=False, active_only=True)
    VoigtFit.SaveDataSet("test.dataset", VFdataset)
    VFdataset.reset_components()

    doppler_b = 1.
    logN = 13.
    coeff = 1.13e+20

    for v in vel_comp:
        for c2line in C2lines:
            VFdataset.add_component_velocity(c2line.level, v, doppler_b, logN)
            baselevel = c2line.level
            break

    copiedlevel = [baselevel]
    for c2line in C2lines:
        if c2line.level not in copiedlevel and c2line.level in C2fitlevel:
            print(c2line.level)
            VFdataset.copy_components(c2line.level, baselevel)
            copiedlevel.append(c2line.level)

    VFdataset.prepare_dataset(norm=False, active_only=True)
    popt, chi2 = VFdataset.fit(verbose=True)
    VFdataset.plot_fit(filename=filename[0], loc="right")
    VFdataset.print_results(velocity=True)

    wf = open(filename[1], "w")

    # fit結果を変数に
    params = VFdataset.best_fit

    # データのインデックス。一つしかない場合は0を設定
    index = 0

    for c2line in C2available:
        # DataSetインスタンスからfitに使った各ラインについてスペクトルデータ(Region)を抽出
        regions_of_line = VFdataset.find_line(c2line.linetag)
        region = regions_of_line[index]
        x, y, err, mask = region.unpack()
        res = region.res

        N_pix = len(x) * 2
        dx = np.diff(x)[0]

        # DataSetインスタンスからLine情報を抽出
        line = VFdataset.lines[c2line.linetag]
        l0, f, gam = line.get_properties()
        ion = line.ion

        # DataSetインスタンスから成分の数を抽出
        n_comp = len(VFdataset.components[ion])

        # コンポーネントについてループ処理
        for n in range(n_comp):

            # プロファイル計算のための波長情報を作成
            wl_line = np.logspace(np.log10(x.min()), np.log10(x.max()), N_pix)
            pxs = np.diff(wl_line)[0] / wl_line[0] * 299792.458
            front_pad = np.arange(x.min() - 50 * dx, x.min(), dx)
            end_pad = np.arange(x.max(), x.max() + 50 * dx, dx)
            wl_line = np.concatenate([front_pad, wl_line, end_pad])

            # tauを計算するための変数を準備
            tau = np.zeros_like(wl_line)

            # ベストフィットの結果を抽出
            z = params['z%i_%s' % (n, ion)].value
            b = params['b%i_%s' % (n, ion)].value
            logN = params['logN%i_%s' % (n, ion)].value
            logNerr = params['logN%i_%s' % (n, ion)].stderr
            Nfit = 10 ** logN
            if logNerr < 0.5:
                Nfit_err = (10 ** (logN + logNerr) - 10 ** (logN - logNerr)) / 2.
            else:
                Nfit_err = 999e+12

            # Voigt profileを計算
            tau += voigt.Voigt(wl_line, l0, f,
                               10 ** logN, 1.e5 * b, gam, z=z)

            # プロファイルを作成
            profile_int = np.exp(-tau)
            fwhm_instrumental = res
            sigma_instrumental = fwhm_instrumental / 2.35482 / pxs
            LSF = gaussian(len(wl_line) / 2, sigma_instrumental)
            LSF = LSF / LSF.sum()
            profile_broad = fftconvolve(profile_int, LSF, 'same')
            profile = profile_broad[50:-50]
            wl_line = wl_line[50:-50]
            vel_profile = (wl_line - l0) / l0 * 299792.458

            # EWを計算
            ew = 0.
            for k in range(1, len(wl_line) - 1):
                if wl_line[k] < l0 + 1. and wl_line[k] > l0 - 1.:
                    ew += ((wl_line[k + 1] + wl_line[k]) / 2. - (wl_line[k] + wl_line[k - 1]) / 2.) * (1. - profile[k])

            peak = np.amax(1. - profile)
            g_sigma = ew / sqtwopi / peak / l0 * light_sp

            # optically thin limitにおける柱密度を計算
            Noptthin = math.log10(coeff * ew / (l0 ** 2) / f)

            if logNerr < 0.5:
                ew_upperr = 10 ** (Noptthin + logNerr) / coeff * (l0 ** 2) * f - ew
                ew_lowerr = 10 ** (Noptthin - logNerr) / coeff * (l0 ** 2) * f - ew
            else:
                ew_upperr, ew_lowerr = 0.999, -0.999

            # print linedict[line_tag],"Comp.",n+1,"EW=%.2f mA, " % (ew * 1e+3), "N=%.2e cm^{-2}" % Noptthin
            wf.write("%s Comp. %d\t%.2f\t%.2f\t%.2f\t%.3f\t%.3f\t(%.3f)\t%.3e\t%.3e\t%.3e\t%.3e\n" %
                     (c2line.transition, n + 1, ew * 1e+3, ew_upperr * 1e+3, ew_lowerr * 1e+3, logN, logNerr, Noptthin,
                      Nfit, Nfit_err, peak, g_sigma))



    wf.close()