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

class NicerExtractor(object):
    #
    # class to extract data products from NuSTAR observations
    #

    def __init__(self, obsdir, region_file='regions.sh', run_reduction=False, suffix=None):
        self.obsdir = obsdir

        self.stem = 'nu%s' % obsdir

        self.evlsdir = self.obsdir + '/xti/event_cl' + ('_%s' % suffix if suffix is not None else '')
        self.specdir = self.obsdir + '/xti/spectras' + ('_%s' % suffix if suffix is not None else '') + ('_%s' % region_suffix if region_suffix is not None else '')
        self.lcdir = self.obsdir + '/xti/lightcurves' + ('_%s' % suffix if suffix is not None else '') + ('_%s' % region_suffix if region_suffix is not None else '')
        self.gtidir = self.obsdir + '/gti'

        #
        # check that the observation is ready to go and if not, run nupipeline
        #
        if(not self._check_reproc()):
            if run_reduction:
                self.reprocess()
            else:
                raise AssertionError("nupipeline has not been run on OBSID " + obsdir)

        # and lastly, the start and stop times of each pointing/event list
        self.start_time = []
        self.stop_time = []
        #self._get_start_stop_times()

        self._ready_to_extract()


    #-- Observation status checks --------------------------------------------

    def _check_reproc(self):
        return len(glob.glob(self.evlsdir + '/*01_cl.evt')) > 0

    def _ready_to_extract(self):
        #
        # check if this observation is ready to extract products
        #
        if(not self._check_reproc()):
            print("Not ready to extract products: Could not find reprocessed event lists.")
            return False
        if(not self._check_regions()):
            print("Not ready to extract products: Regions have not been defined.")
            return False

        return True

    #-- Initial data reduction -----------------------------------------------

    def reprocess(self):
        #
        # argument array for Popen
        #
        args = ['nicerl2',
                'indir=' + self.obsdir,
                'cldir=' + self.evlsdir,
                'clobber=YES'
                ]

        args += ['%s=%s' % (arg, str(kwargs[arg])) for arg in kwargs]

        #
        # and execute it
        #
        proc = subprocess.Popen(args).wait()

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
    def get_spectrum(self, extractdir=None, bkg='scorpeon', suffix=None, **kwargs):
        if not self._ready_to_extract():
            return

        #
        # argument array for Popen
        #
        args = ['nicerl3-spect',
                self.obsdir,
                'cldir=' + self.evlsdir,
                'bkgmodeltype=' + bkg,
                'clobber=YES'
                ]
        if suffix is not None:
            args.append('suffix=_' + suffix)

        args += ['%s=%s' % (arg, str(kwargs[arg])) for arg in kwargs]

        #
        # and execute it
        #
        proc = subprocess.Popen(args).wait()


    def get_lightcurve(self, tbin=10., pirange=(30,800), bkg=None, suffix=None, **kwargs):
        if not self._ready_to_extract():
            return

        if extractdir is None:
            extractdir = self.lcdir

        if not os.path.exists(extractdir):
            os.mkdir(extractdir)

            #
            # argument array for Popen
            #
            args = ['nicerl3-lc',
                    self.obsdir,
                    'cldir=' + self.evlsdir,
                    'pirange=%d-%d' % pirange,
                    'timebin=%g' % tbin,
                    'clobber=YES'
                    ]
            if bkg is not None:
                args.append('bkgmodeltype=%s' % bkg)
            if suffix is not None:
                args.append('suffix=_' + suffix)

            args += ['%s=%s' % (arg, str(kwargs[arg])) for arg in kwargs]

            #
            # and execute it
            #
            proc = subprocess.Popen(args).wait()

    def get_energy_lightcurve(self, enmin, enmax, tbin=10, corr_energy=5, lcfiles=None, bkglcfiles=None, sumlcfile=None,
                              sumbkglcfile=None, bkgsublcfiles=None, sumbkgsublcfile=None, extractdir=None, **kwargs):
        if extractdir is None:
            extractdir = self.lcdir

        if lcfiles is None:
            lcfiles = ["%s_src_%s_tbin%d_en%d-%d.lc" % (self.stem, inst, tbin, enmin, enmax) for inst in ['FPMA', 'FPMB']]
        if bkglcfiles is None:
            bkglcfiles = ["%s_bkg_%s_tbin%d_en%d-%d.lc" % (self.stem, inst, tbin, enmin, enmax) for inst in ['FPMA', 'FPMB']]
        if bkgsublcfiles is None:
            bkgsublcfiles = ["%s_src_bkgsub_%s_tbin%d_en%d-%d.lc" % (self.stem, inst, tbin, enmin, enmax) for inst in ['FPMA', 'FPMB']]
        if sumlcfile is None:
            sumlcfile = "%s_src_%s_tbin%d_en%d-%d.lc" % (self.stem, 'fpmsum', tbin, enmin, enmax)
        if sumbkglcfile is None:
            sumbkglcfile = "%s_bkg_%s_tbin%d_en%d-%d.lc" % (self.stem, 'fpmsum', tbin, enmin, enmax)
        if sumbkgsublcfile is None:
            sumbkgsublcfile = "%s_src_bkgsub_%s_tbin%d_en%d-%d.lc" % (self.stem, 'fpmsum', tbin, enmin, enmax)


        corr_energy = 0.5 * (enmin + enmax) / 1000. # apply PSF vignetting correction at the midpoint of the energy band

        self.get_lightcurve(tbin, corr_energy=corr_energy, lcfiles=lcfiles, bkglcfiles=bkglcfiles,
                            sumlcfile=sumlcfile, sumbkglcfile=sumbkglcfile, bkgsublcfiles=bkgsublcfiles,
                            sumbkgsublcfile=sumbkgsublcfile, extractdir=extractdir,
                            pilow=self.pichan(enmin), pihigh=self.pichan(enmax), **kwargs)

    #-- Utility functions -----------------------------------------------
    @staticmethod
    def pichan(energy):
        """
        Convert energy in eV to NuSTAR PHA channel
        """
        return int(energy / 10)

    #-- Batch light curve generation ------------------------------------
    def batch_energy_lightcurve(self, tbin=10, en0=300, enmax=10000, Nen=10, enbins=None, subdir='energy', extractdir=None, **kwargs):
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

        if extractdir is None:
            extractdir = self.lcdir + '/' + subdir

        for s, e in zip(enstart, enend):
            self.get_energy_lightcurve(s, e, tbin, extractdir=extractdir, **kwargs)

