from .epic_extractor import *

lagen_bins = ([300,400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000],
             [400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000,10000])

lagen_bins_coarse = ([300,400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000],
                     [400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000,10000])

lagen_bins_finefek = ([300,400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,6000,7000,8000],
                      [400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,6000,7000,8000,10000])

lagfreq_bins = ([300,300,1000,1200,4000],
                [1000,800,4000,4000,7000])


def get_lc_for_lagen(tbin=10, enbins=lagen_bins, nogaps=False, lcdir='for_lagen'):
    for o in list_obsids():
        print("Processing OBSID " + str(obsid) + "...")

        extractor = EPICExtractor(obsid)

        if nogaps:
            extractor.evls = extractor.filt_evls

        for emin, emax in zip(enbins[0], enbins[1]):
            print(obsid + ": Getting light curves for energy range {0}-{1}".format(emin, emax))

            extract_dir = extractor.lcdir + '/' + lcdir
            extractor.get_energy_lightcurve(emin, emax, tbin, extract_dir)


def get_lc_for_lagfreq(tbin=10, enbins=lagfreq_bins, nogaps=False, lcdir='for_lagfreq'):
    for o in list_obsids():
        print("Processing OBSID " + str(obsid) + "...")

        extractor = EPICExtractor(obsid)

        if nogaps:
            extractor.evls = extractor.filt_evls

        for emin, emax in zip(enbins[0], enbins[1]):
            print(obsid + ": Getting light curves for energy range {0}-{1}".format(emin, emax))

            extract_dir = extractor.lcdir + '/' + lcdir
            extractor.get_energy_lightcurve(emin, emax, tbin, extract_dir)
