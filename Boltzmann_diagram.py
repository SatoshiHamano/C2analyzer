# coding: UTF-8

__author__ = "Satoshi Hamano"

import os
import sys
import np as np
from pyraf import iraf
import datetime

from spectools_v1_1 import openspec
from C2param_parser import *
from C2_parameters import *
from spectra_plotter import *

C2params = moleclinelist("C2")
bandlist = ["AX00", "AX10", "AX20"]


def convert_to_BD_n2(rot_qn, NJ_obs):
    for i in range(len(rot_qn)):
        if rot_qn[i] == 2:
            N2 = NJ_obs[i]

    return -np.log(5. * NJ_obs / (2 * rot_qn + 1) / N2), N2


def convert_to_BDerr(rot_qn, NJ_obs, NJerr_obs):
    for i in range(len(rot_qn)):
        if rot_qn[i] == 2:
            N2 = NJ_obs[i]

    return np.fabs(-np.log(5. * (NJ_obs + NJerr_obs) / (2 * rot_qn + 1) / N2) + np.log(
        5. * (NJ_obs - NJerr_obs) / (2 * rot_qn + 1) / N2))


if __name__ == "__main__":
    filename = sys.argv[1:]

    config = ReadConfig()
    band = config["C2 parameters"]["band"]
    resolution = float(config["Spectrum"]["resolution"])
    light_sp = scipy.constants.c * 1.e-3
    rp = math.sqrt(math.pi)
    vres = light_sp / resolution

    n_comp, vel_comp, _, b_comp, _ = CloudParameterReader(config)
    J, logN, logN_err = ColumnDensityReader(config)

    comp = int(filename[0])

    transition_obs = []
    cd_obs = []
    cderr_obs = []

    for i in range(len(J)):
        if logN[i][comp - 1] != 0:
            cd_obs.append(logN[i][comp - 1])
            cderr_obs.append(logN_err[i][comp - 1])
            transition_obs.append(J[i])

    modelf = open("vanDishoeck_model_mesh_sigma4WJ_I1.dat", "r")
    modell = modelf.readlines()
    modelf.close()

    transition_model = np.array([int(modell[2].split()[j]) for j in range(4, len(modell[2].split()))])
    rote_model = np.array([float(modell[3].split()[j]) for j in range(4, len(modell[2].split()))])

    temp_model = np.array([float(modell[j].split()[1]) for j in range(5, len(modell))])
    density_model = np.array([float(modell[j].split()[0]) for j in range(5, len(modell))])

    solutions = []

    for i in range(5, len(modell)):
        tmp_sol = []
        for j in range(3, len(modell[i].split())):
            tmp_sol.append(float(modell[i].split()[j]))
        tmp_sol_array = np.array(tmp_sol)
        tmp_sol_array[tmp_sol_array < 0.] = 1.0e-10
        solutions.append(convert_to_BD(transition_model, tmp_sol_array))

    solutions_array = np.array(solutions)

    B_rot = 1.8114  # cm^-1
    D_cen = 0.00000692  # cm^-1

    obsdata_rotqn = np.array(transition_obs)
    obsdata_rot_e = np.array(
        [B_rot * obsdata_rotqn[i] * (obsdata_rotqn[i] + 1) - D_cen * obsdata_rotqn[i] ** 2 * (obsdata_rotqn[i] + 1) ** 2
         for i in range(len(obsdata_rotqn))])
    obsdata_cd, obsdata_cd_n2 = convert_to_BD_n2(obsdata_rotqn, np.array(cd_obs))
    obsdata_cderr = convert_to_BDerr(obsdata_rotqn, np.array(cd_obs), np.array(cderr_obs))

    rot_qn_obs = []
    num_data = 0
    for i in range(len(transition_model)):
        rot_qn_obs.append(transition_model[i] in obsdata_rotqn)
        if transition_model[i] in obsdata_rotqn:
            num_data += 1

    rot_qn_obs = np.array(rot_qn_obs)

    wf = open(filename[1], "w")
    reduced_chi = []

    for i in range(len(solutions_array)):
        chi_sq = np.sum(((obsdata_cd - solutions_array[i][rot_qn_obs]) / obsdata_cderr) ** 2)
        wf.write("%.3e\t%.3e\t%.3e\t%.3e\n" % (temp_model[i], density_model[i], chi_sq, chi_sq / (num_data - 3)))
        reduced_chi.append(chi_sq / (num_data - 3))

    smi = reduced_chi.index(np.amin(reduced_chi))
    print("T=%.2e and n=%.3e" % (temp_model[smi], density_model[smi]))
    smchi = np.amin(reduced_chi)

    temp_chi1 = []
    density_chi1 = []
    solution_range = []

    for i in range(len(reduced_chi)):
        if reduced_chi[i] - 1.15 < smchi:
            temp_chi1.append(temp_model[i])
            density_chi1.append(density_model[i])
            solution_range.append(solutions_array[i])

    temp_min = np.min(temp_chi1)
    temp_max = np.max(temp_chi1)
    density_min = np.min(density_chi1)
    density_max = np.max(density_chi1)

    solution_range_array = np.array(solution_range).T
    solution_range_min = []
    solution_range_max = []
    for i in range(len(transition_model)):
        solution_range_min.append(np.min(solution_range_array[i]))
        solution_range_max.append(np.max(solution_range_array[i]))

    rotqn_plot = np.arange(0., 40.1, 2.)
    rote_plot = B_rot * rotqn_plot * (rotqn_plot + 1) - D_cen * rotqn_plot ** 2 * (rotqn_plot + 1) ** 2
    rotqn_plot_odd = np.arange(0., 40., 1.)
    rote_plot_odd = B_rot * rotqn_plot_odd * (rotqn_plot_odd + 1) - D_cen * rotqn_plot_odd ** 2 * (
            rotqn_plot_odd + 1) ** 2

    cd_model = (np.exp(-solutions_array[smi])) * (2 * rotqn_plot + 1) * obsdata_cd_n2 / 5. / 1.e+12
    cd_model_max = (np.exp(-np.array(solution_range_min))) * (2 * rotqn_plot + 1) * obsdata_cd_n2 / 5. / 1.e+12
    cd_model_min = (np.exp(-np.array(solution_range_max))) * (2 * rotqn_plot + 1) * obsdata_cd_n2 / 5. / 1.e+12

    print("Total column density of 12C2 (1e+12 cm-2): %.2f (%.2f - %.2f)" % (
        np.sum(cd_model), np.sum(cd_model_min), np.sum(cd_model_max)))
    print("Summed column density of 12C2 (1e+12 cm-2): %.2f " % (np.sum(cd_obs) / 1.e+12))

    isratio = 50.
    totalcd_12c13c = numpy.sum(cd_model) / isratio * 2.
    print("Total column density of 12C13C (1e+12 cm-2): ", numpy.sum(totalcd_12c13c))

    print("")
    print("Pop. distribution for 12C13C:")

    cd_model_odd_norm = []
    for i in range(len(rotqn_plot_odd) - 1):
        if i % 2 == 0:
            tmpodd = solutions_array[smi][i / 2]
        if i % 2 == 1:
            tmpodd = (solutions_array[smi][(i - 1) / 2] * (rote_plot_odd[i + 1] - rote_plot_odd[i]) +
                      solutions_array[smi][(i + 1) / 2] * (rote_plot_odd[i] - rote_plot_odd[i - 1])) / (
                             rote_plot_odd[i + 1] - rote_plot_odd[i - 1])
        print(i, tmpodd, (numpy.exp(-tmpodd)) * (2 * rotqn_plot_odd[i] + 1) / 5.)
        cd_model_odd_norm.append((numpy.exp(-tmpodd)) * (2 * rotqn_plot_odd[i] + 1) / 5.)

    print("")
    print("Column densities for 12C13C:")

    cd_model_odd = [cd_model_odd_norm[i] / numpy.sum(cd_model_odd_norm) * totalcd_12c13c for i in
                    range(len(cd_model_odd_norm))]
    for i in range(len(cd_model_odd_norm)):
        print("J=%d\t%.2e" % (i, cd_model_odd[i] * 1.e+12))

    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.65])

    ymin = -1
    ymax = 7
    xmin = -100
    xmax = 1550

    Tans = temp_model[smi]
    Terr_upp = numpy.amax(temp_chi1) - Tans
    Terr_low = Tans - numpy.amin(temp_chi1)
    nans = density_model[smi]
    nerr_upp = numpy.amax(density_chi1) - nans
    nerr_low = nans - numpy.amin(density_chi1)

    ax1.plot(rote_model, solutions_array[smi], "k", lw=0.5)
    ax1.plot(rote_model, solution_range_min, "k--", lw=0.5)
    ax1.plot(rote_model, solution_range_max, "k--", lw=0.5)
    ax1.errorbar(obsdata_rot_e, obsdata_cd, yerr=obsdata_cderr, fmt="o", capsize=3., ms=5., mfc="k", mec="k",
                 ecolor="k")
    ax1.text((xmin + xmax) / 2., ymin - (ymax - ymin) / 8. * 1.2, "$E(J'')$ [K]", ha="center", va="center", fontsize=10)
    ax1.text(xmin - 150., (ymin + ymax) / 2., "$-\log (5N(J'') / N(2) / (2J''+1))$", rotation=90, ha="center",
             va="center", fontsize=10)
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)

    ax2 = ax1.twiny()
    ax2.set_xticks(rote_plot)  # , rotqn_plot.astype(numpy.int64))
    ax2.set_xticklabels(["0", "", "", "", "", "10", "", "", "", "", "20", "", "", "", "", "30", "", "", "", ""])
    ax2.set_xlim(xmin, xmax)
    ax2.text((xmin + xmax) / 2., ymax + (ymax - ymin) / 8. * 1., "$J''$", ha="center", va="center", fontsize=10)
    plt.text(750, 1.2, "Comp. %d" % comp, fontsize=12)
    if Terr_low == Terr_upp:
        plt.text(800, 0.5, "$T=%d \pm %d$ K" % (Tans, Terr_upp), fontsize=11)
    else:
        plt.text(800, 0.5, "$T=%d^{+%d}_{-%d}$ K" % (Tans, Terr_upp, Terr_low), fontsize=11)
    if nerr_low == nerr_upp:
        plt.text(800, -0.2, "$n=%d \pm %d$ cm$^{-3}$" % (nans, nerr_upp), fontsize=11)
    else:
        plt.text(800, -0.2, "$n=%d^{+%d}_{-%d}$ cm$^{-3}$" % (nans, nerr_upp, nerr_low), fontsize=11)

    pp = PdfPages(filename[2])
    plt.savefig(pp, format="pdf")
    plt.clf()
    pp.close()

    plt.scatter(temp_chi1, density_chi1)
    plt.savefig("%s_T_n.png" % filename[2].rstrip(".pdf"))

    wf.close()
