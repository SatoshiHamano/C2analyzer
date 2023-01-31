from spectra_plotter import *
from spectools_v1_1 import openspecfits
import matplotlib.pyplot as plt


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
