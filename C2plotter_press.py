# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy.constants
from spectools_v1_1 import openspec
from C2_parameters import *
from C2param_parser import *
from spectra_plotter_press import *
import matplotlib as mpl

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]




def C2marks(ax, band, vshift, textxshift, textyshift=0.003, ycoord=[1.11, 1.08, 1.05, 1.025], marklen=0.01, linew=1.,
            mcolor="k"):
    # ycoord: P, Q, Rlong, Rshort

    C2band = [C2params.bandlist[band].branchlist("P"),
              C2params.bandlist[band].branchlist("Q"),
              C2params.bandlist[band].branchlist("R")]
    branchs = ["$P$", "$Q$", "$R$"]
    lineco = ["k", "k", "k"]

    vf = (1. + vshift / (scipy.constants.c * 1.e-3))
    theta_array = np.arange(0, 1, 0.01) * math.pi + math.pi / 2.

    [lowlim, upplim] = ax.get_xlim()

    for i in range(3):
        lamtmp = []
        if i != 2:
            for line in C2band[i]:
                if lowlim < line.lamair < upplim:
                    plt.plot([line.lamair * vf, line.lamair * vf], [ycoord[i], ycoord[i] + marklen],
                             lw=linew, color=lineco[i])
                    # plt.text(line.lamair * vf + textxshift, ycoord[i] + marklen + textyshift, line.rotJ,
                    #          fontname="Arial", fontsize=10, va="bottom", ha="center")
                lamtmp.append(line.lamair * vf)
            plt.plot([min(lamtmp), max(lamtmp)], [ycoord[i] + marklen, ycoord[i] + marklen], lw=linew, color=lineco[i])
            plt.text(min(lamtmp) - 3., ycoord[i] - 0.0, branchs[i], color=lineco[i])
        else:
            for line in C2band[i]:
                if lowlim < line.lamair < upplim:
                    if line.rotJ <= 6:
                        plt.plot([line.lamair * vf, line.lamair * vf],
                                 [ycoord[i + 1], ycoord[i + 1] + marklen], lw=linew, color=lineco[i])
                        # plt.text(line.lamair * vf + textxshift, ycoord[i + 1] + marklen + textyshift,
                        #          line.rotJ,
                        #          fontname="Arial", fontsize=10, va="bottom", ha="center")
                    else:
                        plt.plot([line.lamair * vf, line.lamair * vf], [ycoord[i], ycoord[i] + marklen],
                                 lw=linew, color=lineco[i])
                        # plt.text(line.lamair * vf + textxshift, ycoord[i] + marklen + textyshift, line.rotJ,
                        #          fontname="Arial", fontsize=10, va="bottom", ha="center", color=lineco[i])
                lamtmp.append(line.lamair * vf)

            plt.plot([min(lamtmp), max(lamtmp)], [ycoord[i] + marklen, ycoord[i] + marklen], lw=linew,
                     color=lineco[i])
            plt.plot([min(lamtmp), lamtmp[0]], [ycoord[i + 1] + marklen, ycoord[i + 1] + marklen], lw=linew,
                     color=lineco[i])
            plt.plot(3. * np.cos(theta_array) + min(lamtmp),
                     0.5 * (ycoord[2] - ycoord[3]) * np.sin(theta_array) + (ycoord[2] + ycoord[3]) / 2. + marklen,
                     color=lineco[i], lw=linew)
            plt.text(min(lamtmp) - 6., (ycoord[2] + ycoord[3]) / 2.  - 0.0, branchs[i], color=lineco[i])


if __name__ == "__main__":

    filename = sys.argv[1:]

    [spfile, telfile, outputfig] = filename[0:3]
    band = "AX00"
    if len(filename) > 4:
        options = filename[4:]
    else:
        options = []

    options.append("one")
    textshift = 0.

    oneflag = False
    if "one" in options:
        oneflag = True

    if band not in bandlist:
        print("Registered bands keywords are: ")
        for i in bandlist:
            print(i)
        sys.exit()

    xlim = {bandlist[0]: [12060, 12180, 12300], bandlist[1]: [10120, 10260, 10400]}
    ylim = [0.75, 1.15]
    thres = 0.6

    spx, spy = openspec(spfile)
    telx, tely = openspec(telfile)
    spy_tm = transmittance_resampling(spx, spy, telx, tely)
    spy_mask = np.ma.masked_where(spy_tm < thres, spy)

    if "setting" in options:
        config = ReadConfig(band=band)
        config.set("Spectrum", "target", spfile)
        config.set("Spectrum", "telluric", telfile)
        config.set("C2 parameters", "band", band)
        OverwriteConfig(config)

    # figw = 0.87
    # figh_tel = 0.0
    # figh_obj = 0.375
    # orig_w = 0.08
    # orig_h = 0.08
    # gap_telobj = 0.0
    # gap_objtel = 0.1
    figw = 0.87
    figh_tel = 0.07
    figh_obj = 0.28
    orig_w = 0.09
    orig_h = 0.08
    gap_telobj = 0.02
    gap_objtel = 0.1
    inset_w = 0.15
    inset_h = 0.17
    inset_orig_w = 0.60
    inset_orig_h = 0.40
    inset_gap = 0.05

    lightsp = scipy.constants.c * 1.e-3

    r0 = C2params.bandlist[band].branchlist("R")[0].lamair
    isolam = vactoair(12094.3509)

    if outputfig.split(".")[-1] == "pdf":
        pp = PdfPages(outputfig)

    if oneflag:
        fig = plt.figure(figsize=(6, 3.5), facecolor="k")
        ax1 = plt.axes([orig_w, orig_h * 2., figw, figh_tel * 2.])
        spectrum_plot(ax1, [telx], [tely], [xlim[band][0], xlim[band][1]], [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                      yaxis_label="Transmittance  ", yticks_val=[0., 0.5, 1.0], linew=[1.], colors=["0.3"], fst=10, fsl=12)

        fig.patch.set_facecolor("k")


        ax2 = plt.axes([orig_w, (orig_h + figh_tel + gap_telobj) * 2., figw, figh_obj * 2.])
        spectrum_plot(ax2, [spx, spx], [spy, spy_mask], [xlim[band][0], xlim[band][1]], ylim, colors=["0.8", "k"],
                      xticks_val=range(xlim[band][0], xlim[band][1], 20),
                      xticks_label=["" for n in range(xlim[band][0], xlim[band][1], 20)],
                      yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1], yaxis_label="Normalized flux", fst=10, fsl=12)
        C2marks(ax2, band, -10., textshift, linew=1.)
        ax2.text(12062, 1.16, "C$_2$ $A$-$X$ (0,0) band", fontsize=13, color="k", ha="left", va="bottom")
        ax2.text(12062, 0.81, "Cyg OB2 No.12\nWINERED HIRES-J mode ($R=68,000$)", fontsize=10, color="k", ha="left", va="top")

        ax3 = plt.axes([inset_orig_w, inset_orig_h, inset_w, inset_h])
        # spectrum_plot(ax3, [spx], [spy], [r0-1.3, r0+0.8], ylim, colors=["k"],
        #               yticks_val=[0.8, 1.0], yticks_label=[0.8, 1.0], xticks_val=[12085, 12087], xticks_label=[12085, 12087])
        spectrum_plot(ax3, [spx], [spy], [r0-1.3, r0+0.8], ylim, colors=["k"],
                      yticks_val=[0.8, 1.0], yticks_label=[0.8, 1.0], xticks_val=[12085, 12086, 12087], xticks_label=["12085", "86", "87"], fst=7)
        ax3.text(r0-1.3, ylim[0]+0.01, " $^{12}$C$_2$", ha="left", va="bottom")
        ax3.plot([r0 * (1. - 15./lightsp), r0 * (1. - 15./lightsp)], [1.02, 1.06], "0.5")
        ax3.plot([r0 * (1. - 10./lightsp), r0 * (1. - 10./lightsp)], [1.02, 1.06], "0.5")
        ax3.plot([r0 * (1. - 4./lightsp), r0 * (1. - 4./lightsp)], [1.02, 1.06], "0.5")

        ax4  = plt.axes([inset_orig_w + inset_w + inset_gap, inset_orig_h, inset_w, inset_h])
        # spectrum_plot(ax4, [spx], [spy], [isolam-1.3, isolam+0.8], [0.98, 1.015], colors=["k"],
        #               yticks_val=[0.99, 1.0, 1.01], yticks_label=[0.99, 1.0, 1.01], xticks_val=[12090, 12091.4], xticks_label=[12090, 12091.4])
        spectrum_plot(ax4, [spx], [spy], [isolam-1.3, isolam+0.8], [0.98, 1.015], colors=["k"],
                      yticks_val=[0.99, 1.0, 1.01], yticks_label=[".99", "1.0", ".01"], xticks_val=[12090, 12091], xticks_label=["12090", "91"], fst=7)
        ax4.text(isolam-1.3, 0.9806, " $^{12}$C$^{13}$C", ha="left", va="bottom")
        ax4.plot([isolam * (1. - 10./lightsp), isolam * (1. - 10./lightsp)], [1.003, 1.006], "0.5")
        ax4.plot([isolam * (1. - 4./lightsp), isolam * (1. - 4./lightsp)], [1.003, 1.006], "0.5")

        # ax2.patch.set_facecolor("k")
        # ax2.tick_params(axis="x", colors="white")
        # ax2.tick_params(axis="y", colors="white")
        # ax3.patch.set_facecolor("k")
        # ax3.tick_params(axis="x", colors="white")
        # ax3.tick_params(axis="y", colors="white")
        # ax4.patch.set_facecolor("k")
        # ax4.tick_params(axis="x", colors="white")
        # ax4.tick_params(axis="y", colors="white")
        #
    else:
        fig = plt.figure(figsize=(10, 10))
        ax1 = plt.axes([orig_w, orig_h, figw, figh_tel])
        spectrum_plot(ax1, [telx], [tely], [xlim[band][1], xlim[band][2]], [0.0, 1.1],
                      yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])

        ax2 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj, figw, figh_obj])
        spectrum_plot(ax2, [spx, spx], [spy, spy_mask], [xlim[band][1], xlim[band][2]], ylim, colors=["0.8", "orange"],
                      xticks_val=range(xlim[band][1], xlim[band][2], 20),
                      xticks_label=["" for i in range(xlim[band][1], xlim[band][2], 20)],
                      yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1])
        C2marks(ax2, band, 0., textshift)

        ax3 = plt.axes([orig_w, orig_h + figh_tel + gap_telobj + gap_objtel + figh_obj, figw, figh_tel])
        spectrum_plot(ax3, [telx], [tely], [xlim[band][0], xlim[band][1]], [0.0, 1.1], xaxis_label="Wavelength ($\AA$)",
                      yaxis_label="Transmittance", yticks_val=[0., 0.5, 1.0], linew=[1.5], colors=["0.3"])
        ax3.plot()

        ax4 = plt.axes(
            [orig_w, orig_h + figh_tel + gap_telobj + gap_objtel + figh_tel + figh_obj + gap_telobj, figw, figh_obj])
        spectrum_plot(ax4, [spx, spx], [spy, spy_mask], [xlim[band][0], xlim[band][1]], ylim, colors=["0.8", "k"],
                      xticks_val=range(xlim[band][0], xlim[band][1], 20),
                      xticks_label=["" for i in range(xlim[band][0], xlim[band][1], 20)],
                      yticks_val=[0.8, 0.9, 1.0, 1.1], yticks_label=[0.8, 0.9, 1.0, 1.1])
        C2marks(ax4, band, 0., textshift)

    if outputfig.split(".")[-1] == "pdf":
        plt.savefig(pp, format="pdf")
        pp.close()
    else:
        plt.savefig(outputfig)
