import numpy as np
# from pyraf import iraf
import scipy.optimize
import astropy.io.fits as fits

#Description:
#   This script contains some basic and useful functions to manipulate the spectrum files.
#
#Updates:
#   ver1.0 made by Hamano in 2017.01.25
#
#   ver1.1 updated by Hamano in 2018.05.24
#       - pyfits --> astropy.io.fits

# iraf.noao()
# iraf.onedspec()

def openspecfits(fitsfile):
    if fitsfile.find("fits") == -1:
        fitsfile += ".fits"
    
    spfits = fits.open(fitsfile)
    splength = spfits[0].header["NAXIS1"]
    spdata = spfits[0].data

    rcrval1 = float(spfits[0].header["CRVAL1"])
    rcdelt1 = float(spfits[0].header["CDELT1"])
    rcrpix1 = float(spfits[0].header["CRPIX1"])
    
    lamx = np.array([rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)])
    spfits.close()
    
    return lamx, spdata, rcrval1, rcdelt1, rcrpix1

def openspectxt(txtfile):
    spx, spy, fluxerror, maskregion = np.loadtxt(txtfile, dtype=[("col1", "f8"), ("col2", "f8"), ("col3", "f8"),
                                                                   ("col4", "int16")], unpack=True)

    return spx, spy, fluxerror, maskregion

def openspec(spfile):
    if spfile.split(".")[-1] == "fits":
        spx, spy, crval1, cdelt1, crpix1 = openspecfits(spfile)
        return spx, spy
    else:
        spx, spy, fluxerror, maskregion = openspectxt(spfile)
        return spx, spy


def binning_spec(lambdax, fluxy, binning_size):

    lambdax_bin = []
    fluxy_bin = []

    for i in range(len(lambdax)//binning_size):
        tmp_x = 0.
        tmp_y = 0.

        for j in range(binning_size):
            tmp_x += lambdax[binning_size*i+j]
            tmp_y += fluxy[binning_size*i+j]

        lambdax_bin.append(tmp_x/binning_size)
        fluxy_bin.append(tmp_y/binning_size)

    return np.array(lambdax_bin), np.array(fluxy_bin)

#
# def PyScombine(inputlist,output):
#
#     flist = ""
#
#     for i in range(len(inputlist)):
#         flist += inputlist[i]+","
#
#     flist = flist.rstrip(",")
#
#     iraf.scombine(flist, output, combine="average", group="all")
