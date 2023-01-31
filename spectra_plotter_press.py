# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import math
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from spectools_v1_1 import openspecfits
import scipy.constants
from scipy import signal, interpolate
import matplotlib as mpl


def spectrum_plot(ax, spx, spy, xrange, yrange, spx_func="NA", spy_func="NA", colors=["k", "b", "r", "g"], linew=[1.],
                  lines=["-"], labels=[""], colors_func=["k", "b", "r", "g"], linew_func=[1.],
                  lines_func=["--"], labels_func=[""], yshift=[0], yfactor=[1.], xlshift=[0], xvshift=[0],
                  xticks_val="NA", yticks_val="NA", xticks_label="NA", yticks_label="NA", xaxis_label="NA",
                  yaxis_label="NA", grid_flag=False, legend_flag=False, legend_loc=0, fst=10, fsl=12):
    # spx, spy: list of numpy.array
    # xrange, yrange: [xmin, xmax], [ymin, ymax]

    # colors, linew, lines, labels: list
    # yshift, yfactor, xlinear, xvshift: list of floats
    # xticks_val, yticks_val, xticks_label, yticks_label: list
    # grid_flag, legend_flag: boolean
    # legend_loc: string or int

    # plt.rcParams['font.family'] = 'YuGothic'

    if type(spx) is not list:
        print("Parameter \"spx\" is not list.")
        sys.exit()
    if type(spy) is not list:
        print("Parameter \"spy\" is not list.")
        sys.exit()

    Nspec = len(spx)

    colors = colors * math.ceil(Nspec / len(colors)) if len(colors) < Nspec else colors
    linew = linew * math.ceil(Nspec / len(linew)) if len(linew) < Nspec else linew
    lines = lines * math.ceil(Nspec / len(lines)) if len(lines) < Nspec else lines
    labels = labels * math.ceil(Nspec / len(labels)) if len(labels) < Nspec else labels
    yshift = yshift + [0] * (Nspec - len(yshift)) if len(yshift) < Nspec else yshift
    yfactor = yfactor + [1] * (Nspec - len(yfactor)) if len(yfactor) < Nspec else yfactor
    xlshift = xlshift + [0] * (Nspec - len(xlshift)) if len(xlshift) < Nspec else xlshift
    xvshift = xvshift + [0] * (Nspec - len(xvshift)) if len(xvshift) < Nspec else xvshift

    for i in range(Nspec):
        spx_plot = (spx[i] + xlshift[i]) * (1. + xvshift[i] / (scipy.constants.c * 1.e-3))
        spy_plot = spy[i] * yfactor[i] + yshift[i]
        ax.step(spx_plot, spy_plot, where="mid", c=colors[i], label=labels[i], lw=linew[i])

    if spx_func != "NA" and spy_func != "NA":
        if type(spx_func) is not list:
            print("Parameter \"spx_func\" is not list.")
            sys.exit()
        if type(spy_func) is not list:
            print("Parameter \"spy_func\" is not list.")
            sys.exit()

        Nspec_func = len(spx_func)

        colors_func = colors_func * math.ceil(Nspec_func / len(colors_func)) if len(
            colors_func) < Nspec_func else colors_func
        linew_func = linew_func * math.ceil(Nspec_func / len(linew_func)) if len(
            linew_func) < Nspec_func else linew_func
        lines_func = lines_func * math.ceil(Nspec_func / len(lines_func)) if len(
            lines_func) < Nspec_func else lines_func
        labels_func = labels_func * math.ceil(Nspec_func / len(labels_func)) if len(
            labels_func) < Nspec_func else labels_func

        for i in range(Nspec_func):
            ax.step(spx_func, spy_func, where="mid", c=colors_func[i], label=labels_func[i], lw=linew_func[i])

    if xticks_val != "NA":
        ax.set_xticks(xticks_val)
        if xticks_label != "NA":
            ax.set_xticklabels(xticks_label, fontsize=fst)
    if yticks_val != "NA":
        ax.set_yticks(yticks_val)
        if yticks_label != "NA":
            ax.set_yticklabels(yticks_label, fontsize=fst)

    if xaxis_label != "NA":
        ax.set_xlabel(xaxis_label, fontsize=fsl)
    if yaxis_label != "NA":
        ax.set_ylabel(yaxis_label, fontsize=fsl)

    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])

    if grid_flag:
        ax.grid()
    if legend_flag:
        ax.legend(loc=legend_loc)

def transmittance_resampling(spx, spy, telx, tely):

    telx_subp = np.linspace(min(telx), max(telx), 10*telx.size)

    it = interpolate.interp1d(telx, tely, kind="cubic")
    tely_subp = it(telx_subp)

    spy_tm = np.copy(spy)
    for i in range(len(spx)):
        spy_tm[i] = tely_subp[np.argmin(np.absolute(telx_subp - spx[i]))]

    return spy_tm


if __name__ == "__main__":
    fig = plt.figure(figsize=(10, 6))

    # [origin_x, origin_y, width, height]
    ax1 = plt.axes([0.08, 0.1, 0.9, 0.15])
    ax2 = plt.axes([0.08, 0.3, 0.9, 0.6])
    ax3 = plt.axes([0.38, 0.35, 0.45, 0.4])

    # open spectrum file
    spx, spy, _, _, _ = openspecfits("VI_Cyg_12_comb_00band_lamdif_m144add_norm.fits")

    # plot spectra
    spectrum_plot(ax1, [spx], [spy], [12100, 12300], [0.5, 1.2], xaxis_label="X jiku", yaxis_label="Y jiku ($\mu$)")

    spectrum_plot(ax2, [spx, spx], [spy, spy], [12100, 12120], [0.5, 1.2], colors=["y", "r"], lines=["--", "-"],
                  linew=[3., 1.], legend_flag=True, labels=["base", "(x,y) (+3,+0.1)"], xlshift=[0., 3.], yshift=[0., 0.1],
                  grid_flag=True)

    spectrum_plot(ax3, [spx, spx, spx], [spy, spy, spy], [12050, 12090], [0.9, 1.1], xvshift=[0., -40., 50.],
                  yshift=[-0.01, 0., 0.01], labels=["base", "-40 km/s shift", "50 km/s shift"], legend_flag=True,
                  xticks_val=[12050, 12080], yticks_val=[0.8, 1.0], legend_loc="upper left")


    plt.savefig("test.png")
