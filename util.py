import os
import glob
import re

from .epic_extractor import *

def list_obsids():
    #
    # Return a list of the OBSIDs directories in the current location
    # (returns list of directories whose names are just 10 digits)
    #
    return sorted([d for d in next(os.walk('.'))[1] if re.match('^[0-9]+$', d)])


def abbrev_obsid(files):
    file_list = sorted(glob.glob(files))
    obsids = list(set([ re.search('^[0-9]+', f).group(0) for f in file_list ]))

    for num, obsid in enumerate(obsids):
        these_files = [ f for f in file_list if obsid in f ]
        for f in these_files:
            newname = f.replace(obsid, str(num))
            os.rename(f, newname)


def get_spectra(instrument='pn'):
    for obsid in list_obsids():
        extractor = EPICExtractor(obsid, instrument=instrument)
        extractor.get_spectrum()


def get_lightcurves(tbin=100, instrument='pn', nogaps=False):
    for obsid in list_obsids():
        extractor = EPICExtractor(obsid, instrument=instrument)

        if nogaps:
            extractor.evls = extractor.filt_evls

        extractor.get_lightcurve(tbin)
