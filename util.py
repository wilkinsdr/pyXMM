import os
import glob
import re

def list_obsids():
    #
    # Return a list of the OBSIDs directories in the current location
    # (returns list of directories whose names are just 10 digits)
    #
    return sorted([d for d in next(os.walk('.'))[1] if re.match('^[0-9]{10}$', d)])


def abbrev_obsid(files):
    file_list = sorted(glob.glob(files))
    obsids = list(set([ re.search('^[0-9]+', f).group(0) for f in file_list ]))

    for num, obsid in enumerate(obsids):
        these_files = [ f for f in file_list if obsid in f ]
        for f in these_files:
            newname = f.replace(obsid, str(num))
            os.rename(f, newname)
