"""
pyXMM OM Utilities
"""
import os
import glob
import re
import subprocess
import tarfile
import numpy as np
try:
    import astropy.io.fits as pyfits
except:
    import pyfits

class OMExtractor(object):
    #
    # class to extract data products from XMM-Newton EPIC pn observations
    #

    def __init__(self, obsdir, run_reduction=False, bkgfilt=True, subdir=None, mode='image'):
        self.obsdir = obsdir

        self.odfdir = self.obsdir + '/odf'
        self.procdir = self.obsdir + '/proc'

        self.omdir = self.procdir + '/om'
        if subdir is not None:
            self.omdir += '/' + subdir

        self.cif = self.procdir + '/ccf.cif'

        #
        # read the current environment and add the required variables for SAS
        #
        self.envvars = dict(os.environ)
        self.envvars['SAS_ODF'] = os.path.abspath(self.odfdir)
        self.envvars['SAS_CCF'] = os.path.abspath(self.cif)

        #
        # check that the observation is ready to go and if not, run the reduction
        #
        # first we need the CIF
        if run_reduction and not self._check_odfs_extracted():
            self.extract_odfs()
        if(not self._check_cif()):
            if run_reduction:
                if not os.path.exists(self.procdir):
                    os.mkdir(self.procdir)
                if not os.path.exists(self.instdir):
                    os.mkdir(self.instdir)
                self.cifbuild()
            else:
                raise AssertionError("CIF has not been built for OBSID " + obsdir)
        # then odfingest needs to have been run
        if(not self._check_odfingest()):
            if run_reduction:
                self.odfingest()
            else:
                raise AssertionError("odfingest has not been run on OBSID " + obsdir)

        if not os.path.exists(self.omdir):
            os.mkdir(self.omdir)

        if run_reduction:
            self.proc_om(mode)

        self.find_files()

    #-- Observation status checks --------------------------------------------

    def _check_odfs_extracted(self):
        return len(glob.glob(self.odfdir + '/*.FIT')) > 0

    def _check_odfingest(self):
        #
        # check odfingest has been run to initialise the observation for reduction
        #
        return (len(glob.glob(self.odfdir + '/*SUM.SAS')) > 0)

    def _check_cif(self):
        #
        # check the calibration information file has been built
        #
        return os.path.exists(self.cif)

    #-- Initial data reduction -----------------------------------------------

    def cifbuild(self):
        #
        # argument array for Popen
        #
        args = ['cifbuild',
                'fullpath=yes'
                ]
        #
        # and execute it
        #
        print("Building calibration index...")

        proc = subprocess.Popen(args, cwd=self.procdir, env=self.envvars).wait()

    def extract_odfs(self):
        #
        # extract ODFs from tarball if needed
        #
        print("Extracting ODFs...")
        if len(glob.glob(self.odfdir + '/*.FIT')) == 0:
            for odftar in glob.glob(self.odfdir + '/*.tar*'):
                tar = tarfile.open(odftar)
                tar.extractall(self.odfdir)
                tar.close()
            # then need to extract the next level of tarballs
            for odftar2 in glob.glob(self.odfdir + '/*.TAR'):
                tar = tarfile.open(odftar2)
                tar.extractall(self.odfdir)
                tar.close()

    def odfingest(self):
        #
        # argument array for Popen
        #
        print("Ingesting ODFs...")

        args = ['odfingest',
                'odfdir=' + self.odfdir,
                'outdir=' + self.odfdir
                ]
        #
        # and execute it
        #
        proc = subprocess.Popen(args, env=self.envvars).wait()

    #-- EPIC pipeline reduction ----------------------------------------------

    def proc_om(self, mode='image'):
        if mode == 'image':
            args = ['omichain']
        elif mode == 'fast':
            args = ['omfchain']

        print("Running %s to process RGS data..." % args[0])

        proc = subprocess.Popen(args, cwd=self.omdir, env=self.envvars).wait()

        self.find_files()

    def find_files(self):
        self.lcfiles = sorted(glob.glob(self.omdir + '/*TIMESR*.FIT'))

    def concatenate_lc(self, om_filter=None, rebin=None, outfile=None):
        from pylag.lightcurve import LightCurve

        lclist = []

        for lcfile in self.lcfiles:
            with pyfits.open(lcfile) as f:
                if om_filter is not None and f[0].header['FILTER'] != om_filter:
                    continue
                t = np.array(f['RATE'].data['TIME'])
                r = np.array(f['RATE'].data['RATE'])
                e = np.array(f['RATE'].data['ERROR'])

            this_lc = LightCurve(t=t, r=r, e=e)
            if rebin is not None:
                this_lc = this_lc.rebin(rebin)
            lclist.append(this_lc)

        lc = LightCurve().concatenate(lclist)

        if outfile is None:
            namearr = [self.obsdir, 'om', om_filter]
            name = '_'.join(filter(None, namearr))
            outfile = self.omdir + '/' + name + '.lc'
        lc.write_fits(outfile)

