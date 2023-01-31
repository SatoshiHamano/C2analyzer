# -*- coding:utf-8 -*-

from pyraf import iraf
import math
import sys
import numpy

iraf.imred()
iraf.echelle()

if __name__ == "__main__":
    filename = sys.argv[1:]

    iraf.continuum(filename[0], filename[1], sample=filename[2], override="yes", interactive="yes", func="legendre",
                   high_rej="3", order="3")
