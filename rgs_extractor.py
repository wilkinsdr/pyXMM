"""
pyXMM RGS utilities

Temproary instructions
1) initialise RGSExtractor
2) proc_rgs()
3) filter()
4) find_spectra()
5) combine_spectra()
"""
import os
import glob
import re
import subprocess
import tarfile
try:
    import astropy.io.fits as pyfits
except:
    import pyfits

class RGSExtractor(object):
    #
    # class to extract data products from XMM-Newton EPIC pn observations
    #

    def __init__(self, obsdir, run_reduction=False, bkgfilt=True, bkg_rate=0.2, subdir=None):
        self.obsdir = obsdir

        self.odfdir = self.obsdir + '/odf'
        self.procdir = self.obsdir + '/proc'

        self.rgsdir = self.procdir + '/rgs'
        if subdir is not None:
            self.rgsdir += '/' + subdir

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

        if not os.path.exists(self.rgsdir):
            os.mkdir(self.rgsdir)

        self.find_files()
        if len(self.evls) == 0:
            print("RGS event lists not found. Have you run proc_rgs()?")
        if len(self.src_lists) == 0:
            print("RGS source lists not found. Have you run proc_rgs()?")

        self.find_spectra()

        if len(self.src_spectra_o1) == 0 and run_reduction:
            self.proc_rgs()
            if bkg_filt:
                self.filter(rate=bkg_rate)
            self.find_spectra()
            self.combine_spectra()

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

    def proc_rgs(self):
        args = ['rgsproc']

        print("Running %s to process RGS data..." % args[0])

        proc = subprocess.Popen(args, cwd=self.rgsdir, env=self.envvars).wait()

        self.find_files()

    def find_files(self):
        self.evls = sorted(glob.glob(self.rgsdir + '/*EVENLI0000.FIT'))
        self.src_lists = sorted(glob.glob(self.rgsdir + '/*SRCLI_0000.FIT'))


    #-- Event filtering ------------------------------------------------------

    def filter(self, rate=0.2, time=None):
        #
        # run the RGS filters based on a GTI created from the background rate
        #
        print("Filtering RGS data...")

        # just filter based on the first event list (should be RGS1)
        evl = self.evls[0]
        srclist = self.src_lists[0]

        filt = ['(CCDNR==9)', '(REGION(%s:RGS1_BACKGROUND,M_LAMBDA,XDSP_CORR))' % srclist]
        expr = '&&'.join(filt)

        lcfile = evl.replace('EVENLI0000.FIT', '_bkg.lc')
        gtifile = evl.replace('EVENLI0000.FIT', '_bkg.gti')

        args = ['evselect',
                'table='+evl,
                'timebinsize=100',
                'makeratecolumn=yes',
                'maketimecolumn=yes',
                'expression='+expr,
                'rateset='+lcfile]
        proc = subprocess.Popen(args, env=self.envvars).wait()

        gtiexpr = '(RATE<%g)' % rate

        if isinstance(time, tuple):
            gtiexpr += '&&(TIME>=%g)&&(TIME<%g)' % time

        args = ['tabgtigen',
                'table='+lcfile,
                'gtiset='+gtifile,
                'expression='+gtiexpr]
        proc = subprocess.Popen(args, env=self.envvars).wait()

        args = ['rgsproc',
                'entrystage=3:filter',
                'auxgtitables='+os.path.basename(gtifile)]
        proc = subprocess.Popen(args, cwd=self.rgsdir, env=self.envvars).wait()

    # -- Spectra -------------------------------------------------------------

    def find_spectra(self):
        self.src_spectra_o1 = sorted(glob.glob(self.rgsdir + '/*SRSPEC1001.FIT'))
        self.src_spectra_o2 = sorted(glob.glob(self.rgsdir + '/*SRSPEC2001.FIT'))
        self.bkg_spectra_o1 = sorted(glob.glob(self.rgsdir + '/*BGSPEC1001.FIT'))
        self.bkg_spectra_o2 = sorted(glob.glob(self.rgsdir + '/*BGSPEC2001.FIT'))
        self.rsp_o1 = sorted(glob.glob(self.rgsdir + '/*RSPMAT1001.FIT'))
        self.rsp_o2 = sorted(glob.glob(self.rgsdir + '/*RSPMAT2001.FIT'))
    def combine_spectra(self):
        pha = " ".join([os.path.basename(s) for s in self.src_spectra_o1])
        bkg = " ".join([os.path.basename(s) for s in self.bkg_spectra_o1])
        rsp = " ".join([os.path.basename(r) for r in self.rsp_o1])

        comb_pha = "%s_src_rgs_comb_o1.pha" % self.obsdir
        comb_bkg = "%s_bkg_rgs_comb_o1.pha" % self.obsdir
        comb_rsp = "%s_src_rgs_comb_o1.rsp" % self.obsdir

        args = ['rgscombine',
                'pha='+pha,
                'rmf='+rsp,
                'bkg='+bkg,
                'filepha='+comb_pha,
                'filermf='+comb_rsp,
                'filebkg='+comb_bkg]
        proc = subprocess.Popen(args, cwd=self.rgsdir, env=self.envvars).wait()

        pha = " ".join([os.path.basename(s) for s in self.src_spectra_o2])
        bkg = " ".join([os.path.basename(s) for s in self.bkg_spectra_o2])
        rsp = " ".join([os.path.basename(r) for r in self.rsp_o2])

        comb_pha = "%s_src_rgs_comb_o2.pha" % self.obsdir
        comb_bkg = "%s_bkg_rgs_comb_o2.pha" % self.obsdir
        comb_rsp = "%s_src_rgs_comb_o2.rsp" % self.obsdir

        args = ['rgscombine',
                'pha='+pha,
                'rmf='+rsp,
                'bkg='+bkg,
                'filepha='+comb_pha,
                'filermf='+comb_rsp,
                'filebkg='+comb_bkg]
        proc = subprocess.Popen(args, cwd=self.rgsdir, env=self.envvars).wait()


def rgs_collect_spectra(dest='combined/rgs/spectra'):
    import shutil
    from .util import list_obsids

    if not os.path.exists(os.path.dirname(dest)):
        os.mkdir(os.path.dirname(dest))
    if not os.path.exists(dest):
        os.mkdir(dest)

    for o in list_obsids():
        x = RGSExtractor(o)
        for f in x.src_spectra_o1:
            shutil.copy(f, dest)
        for f in x.bkg_spectra_o1:
            shutil.copy(f, dest)
        for f in x.rsp_o1:
            shutil.copy(f, dest)
        for f in x.src_spectra_o2:
            shutil.copy(f, dest)
        for f in x.bkg_spectra_o2:
            shutil.copy(f, dest)
        for f in x.rsp_o2:
            shutil.copy(f, dest)

def rgs_combine_spectra(dir='combined/rgs/spectra', prefix=''):
    pha = " ".join([os.path.basename(s) for s in sorted(glob.glob(dir + '/P%s*SRSPEC1001.FIT' % prefix))])
    bkg = " ".join([os.path.basename(s) for s in sorted(glob.glob(dir + '/P%s*BGSPEC1001.FIT' % prefix))])
    rsp = " ".join([os.path.basename(s) for s in sorted(glob.glob(dir + '/P%s*RSPMAT1001.FIT' % prefix))])

    comb_pha = "%ssrc_rgs_comb_o1.pha" % (prefix + '_' if prefix != '' else '')
    comb_bkg = "%sbkg_rgs_comb_o1.pha" % (prefix + '_' if prefix != '' else '')
    comb_rsp = "%ssrc_rgs_comb_o1.rsp" % (prefix + '_' if prefix != '' else '')

    args = ['rgscombine',
            'pha=' + pha,
            'rmf=' + rsp,
            'bkg=' + bkg,
            'filepha=' + comb_pha,
            'filermf=' + comb_rsp,
            'filebkg=' + comb_bkg]
    proc = subprocess.Popen(args, cwd=dir).wait()

    pha = " ".join([os.path.basename(s) for s in sorted(glob.glob(dir + '/P%s*SRSPEC2001.FIT' % prefix))])
    bkg = " ".join([os.path.basename(s) for s in sorted(glob.glob(dir + '/P%s*BGSPEC2001.FIT' % prefix))])
    rsp = " ".join([os.path.basename(s) for s in sorted(glob.glob(dir + '/P%s*RSPMAT2001.FIT' % prefix))])

    comb_pha = "%ssrc_rgs_comb_o2.pha" % (prefix + '_' if prefix != '' else '')
    comb_bkg = "%sbkg_rgs_comb_o2.pha" % (prefix + '_' if prefix != '' else '')
    comb_rsp = "%ssrc_rgs_comb_o2.rsp" % (prefix + '_' if prefix != '' else '')

    args = ['rgscombine',
            'pha=' + pha,
            'rmf=' + rsp,
            'bkg=' + bkg,
            'filepha=' + comb_pha,
            'filermf=' + comb_rsp,
            'filebkg=' + comb_bkg]
    proc = subprocess.Popen(args, cwd=dir).wait()
