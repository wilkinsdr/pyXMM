import os
import glob
import re
import subprocess
import tarfile
try:
    import astropy.io.fits as pyfits
except:
    import pyfits
import numpy as np

from .spec_util import *
from .gti import *

class NustarExtractor(object):
    #
    # class to extract data products from XMM-Newton EPIC pn observations
    #

    def __init__(self, obsdir, region_file='regions.sh', run_reduction=False, saacalc=3, saamode=None, tentacle=None, suffix=None, region_suffix=None):
        self.obsdir = obsdir

        self.stem = 'nu%s' % obsdir

        self.evlsdir = self.obsdir + '/reproc' + ('_%s' % suffix if suffix is not None else '')
        self.extractdir = self.obsdir + '/products' + ('_%s' % suffix if suffix is not None else '') + ('_%s' % region_suffix if region_suffix is not None else '')
        self.regiondir = self.obsdir + '/regions'
        self.gtidir = self.obsdir + '/gti'

        self.saacalc = saacalc
        self.saamode = saamode
        self.tentacle = tentacle

        #
        # check that the observation is ready to go and if not, run the reduction
        #
        # then odfingest needs to have been run
        if(not self._check_reproc()):
            if run_reduction:
                self.reprocess(saacalc, saamode, tentacle)
            else:
                pass
                #raise AssertionError("nupipeline has not been run on OBSID " + obsdir)


        # the region specification for each pointing/event list
        self.regions = self._get_regions(region_suffix)

        # and lastly, the start and stop times of each pointing/event list
        self.start_time = []
        self.stop_time = []
        #self._get_start_stop_times()


    #-- Observation status checks --------------------------------------------

    def _check_reproc(self):
        return len(glob.glob(self.evlsdir + '*01_cl.evt')) > 0

    def _ready_to_extract(self):
        #
        # check if this observation is ready to extract products
        #
        if(not self._check_odfingest()):
            print(
                "Not ready to extract products: odfingest has not been run on these ODFs")
            return False
        if(not self._check_cif()):
            print("Not ready to extract products: CIF has not been built")
            return False
        if(len(self.evls) < 1):
            print("Not ready to extract products: Could not find filtered event list")
            return False
        if(len(self.regions) < 1):
            print("Not ready to extract products: No regions defined")
            return False

        return True

    #-- Initial data reduction -----------------------------------------------

    def reprocess(self, saacalc=3, saamode=None, tentacle=None):
        #
        # argument array for Popen
        #
        args = ['nupipeline',
                'indir=' + self.obsdir,
                'steminputs=' + self.stem,
                'outdir=' + self.evlsdir
                ]

        if saamode is not None:
            args += ['saamode=%s' % saamode,
                     'saacalc=%d' % saacalc]
        if tentacle is not None:
            args += ['tentacle=%s' % tentacle]
        #
        # and execute it
        #
        proc = subprocess.Popen(args).wait()


    #-- Region lookup  -------------------------------------------------------
    def _get_regions(self, region_suffix=None):
        regions = {'FPMA' : {'src': self.regiondir + '/src_FPMA%s.reg' % ('_%s' % region_suffix if region_suffix is not None else ''),
                            'bkg': self.regiondir + '/bkg_FPMA%s.reg' % ('_%s' % region_suffix if region_suffix is not None else '')},
                   'FPMB' : {'src': self.regiondir + '/src_FPMB%s.reg' % ('_%s' % region_suffix if region_suffix is not None else ''),
                            'bkg': self.regiondir + '/bkg_FPMB%s.reg' % ('_%s' % region_suffix if region_suffix is not None else '')}}
        return regions


    #-- Observation header data lookups --------------------------------------

    def _get_start_stop_times(self):
        #
        # return the start time of the observation (in satellite time)
        # specified in the event list FITS header
        #
        for evl in self.evls:
            if(not os.path.exists(evl)):
                print("StartTime ERROR: Could not open event list")
                return

            f = pyfits.open(evl)
            self.start_time.append(f[1].header['TSTART'])
            self.stop_time.append(f[1].header['TSTOP'])
            f.close()

    #-- Product extraction routines -----------------------------------------
    def nuproducts(self, extractdir=None, instruments=None, **kwargs):
        if extractdir is None:
            extractdir = self.extractdir

        if instruments is None:
            instruments = ['FPMA', 'FPMB']
        elif isinstance(instruments, str):
            instruments = [instruments]

        for inst in instruments:
            #
            # argument array for Popen
            #
            args = ['nuproducts',
                    'srcregionfile=' + self.regions[inst]['src'],
                    'bkgregionfile=' + self.regions[inst]['bkg'],
                    'indir=' + self.evlsdir,
                    'outdir=' + extractdir,
                    'instrument=' + inst,
                    'steminputs=' + self.stem,
                    ]
            if self.regions[inst]['bkg'] is not None:
                args += ['bkgextract=yes']

            args += ['%s=%s' % (arg, str(kwargs[arg])) for arg in kwargs]

            #
            # and execute it
            #
            proc = subprocess.Popen(args).wait()

    def extract_products(self, extractdir=None, grptype='opt', grpscale=None, usebkg=False, **kwargs):
        if extractdir is None:
            extractdir = self.extractdir

        self.nuproducts(extractdir, **kwargs)
        if grptype is not None:
            self.group_spectra(extractdir, grptype, grpscale, usebkg)

    def group_spectra(self, extractdir=None, grptype='opt', grpscale=None, usebkg=False):
        if extractdir is None:
            extractdir = self.extractdir

        for inst in ['A', 'B']:
            srcspec = extractdir + self.stem + inst + '01_sr.pha'
            bkgspec = (extractdir + self.stem + inst + '01_sr.pha') if usebkg else None
            rmf = extractdir + self.stem + inst + '01_sr.rmf'
            grpspec = extractdir + self.stem + inst + '01_sr.grp'

            group_spec(grpspec, srcspec, bkgfile=bkgspec, rmffile=rmf, grptype=grptype, grpscale=grpscale, usebkg=usebkg)

    def extract_products_time(self, tstart, tstop, label=None, extractdir=None, **kwargs):
        if label is None:
            label = 'time%g-%g' % (tstart, tstop)

        if extractdir is None:
            extractdir = self.extractdir + '_' + label

        gtifile = self.gtidir + label + '.gti'

        g = GTI()
        g.add_row(tstart, tstop)
        g.write(gtifile)

        self.extract_products(extractdir=extractdir, usrgtifile=gtifile, **kwargs)

    def extract_products_gti(self, gtifile, label=None, extractdir=None, **kwargs):
        if label is None:
            label = os.path.splitext(os.path.basename(gtifile))[0]

        if extractdir is None:
            extractdir = self.extractdir + '_' + label

        self.extract_products(extractdir=extractdir, usrgtifile=gtifile, **kwargs)

    def get_lightcurve(self, tbin=10, corr_energy=5, lcfiles=None, bkglcfiles=None, extractdir=None, **kwargs):
        if lcfiles is None:
            lcfiles = ["%s_src_%s_tbin%d.lc" % (self.stem, inst, tbin) for inst in ['FPMA', 'FPMB']]
        if lcfiles is None:
            bkglcfiles = ["%s_bkg_%s_tbin%d.lc" % (self.stem, inst, tbin) for inst in ['FPMA', 'FPMB']]

        for inst, lcfile, bkglcfile in zip(['FPMA', 'FPMB'], lcfiles, bkglcfiles:
            self.nuproducts(extractdir, instruments=[inst], lcfile=lcfile, bkglcfile=bkglcfile, imagefile='NONE', phafile='NONE', runmkarf='no', runmkrmf='no', binsize=tbin, lcenergy=corr_energy, **kwargs)


    #-- Batch light curve generation ------------------------------------
    def batch_energy_lightcurve(self, tbin=10, en0=300, enmax=10000, Nen=10, enbins=None, subdir='energy', **kwargs):
        if enbins is None:
            enbins = np.round(np.logspace(np.log10(en0), np.log10(enmax), Nen)).astype(int)
        elif enbins == 'lagen':
            enbins = [300,400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000,10000]
        elif enbins == 'lagen_coarse':
            enbins = [300,500,700,1000,1400,2000,3000,5000,7000,10000]
        elif enbins == 'lagfreq':
            enbins = ([300,300,1000,1200,4000], [800,1000,4000,4000,7000])

        if isinstance(enbins, tuple):
            enstart = enbins[0]
            enend = enbins[1]
        else:
            enstart = enbins[:-1]
            enend = enbins[1:]

        extract_dir = self.lcdir + '/' + subdir

        for s, e in zip(enstart, enend):
            self.get_energy_lightcurve(s, e, tbin, extract_dir, **kwargs)

