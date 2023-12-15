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

class EPICExtractor(object):
    #
    # class to extract data products from XMM-Newton EPIC pn observations
    #

    def __init__(self, obsdir, instrument='pn', region_file='regions.sh', run_reduction=False, bkgfilt=True, pileup_corr=False, ccdnr=4):
        self.obsdir = obsdir

        self.odfdir = self.obsdir + '/odf'
        self.procdir = self.obsdir + '/proc'

        self.instrument = instrument

        if(instrument == 'pn'):
            self.instdir = self.procdir + '/pn'
        elif(instrument == 'mos'):
            self.instdir = self.procdir + '/mos'

        self.evlsdir = self.instdir + '/evl'
        self.regiondir = self.instdir + '/regions'
        self.specdir = self.instdir + '/spectra'
        self.lcdir = self.instdir + '/lightcurves'
        self.imagedir = self.instdir + '/images'
        self.epatdir = self.instdir + '/epat'
        self.gtidir = self.instdir + '/gti'
        
        if pileup_corr:
            self.evlsdir += '_pileupcorr'
            self.specdir += '_pileupcorr'

        self.pileup_corr = pileup_corr
        self.ccdnr = ccdnr

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

        self.get_evls(bkgfilt=bkgfilt, run_reduction=run_reduction, pileupcorr=pileup_corr)

        # the region specification for each pointing/event list
        self._get_regions(region_file)

        # and lastly, the start and stop times of each pointing/event list
        self.start_time = []
        self.stop_time = []
        self._get_start_stop_times()


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

    #-- Find event lists and reduce/filter them if necessary -----------------
    def get_evls(self, bkgfilt=True, run_reduction=False, pileupcorr=False):
        # first check if event lists exusr
        self.unfilt_evls = sorted(glob.glob(self.evlsdir + '/*ImagingEvts.ds'))
        self.filt_evls = sorted(glob.glob(self.evlsdir + '/*filtered.evl'))
        self.bkgfilt_evls = sorted(glob.glob(self.evlsdir + '/*filtered?gti.evl'))
        self.raw_evls = sorted(glob.glob(self.evlsdir + '/*PileupEvts.ds'))

        # read off the section of the filename that identifies the instrument and pointing
        # we do this from the filtered list first in case any have been disabled by renaming the files
        if (len(self.bkgfilt_evls) > 0):
            self.pointings = [re.search('(E[A-Z]+[0-9]_)*[A-Z][0-9]{3}', f).group(0) for f in self.bkgfilt_evls]
        else:
            self.pointings = [re.search('(E[A-Z]+[0-9]_)*[A-Z][0-9]{3}', f).group(0) for f in self.unfilt_evls]

        if (len(self.bkgfilt_evls) < 1):
            # if we dont' have filtered event lists, see what we do have and build
            # them if necessary
            if (len(self.unfilt_evls) < 1):
                if run_reduction:
                    if not os.path.exists(self.instdir):
                        os.mkdir(self.instdir)
                    if not os.path.exists(self.evlsdir):
                        os.mkdir(self.evlsdir)
                    self.proc_evl(pileupcorr=pileupcorr)
                    if not os.path.exists(self.regiondir):
                        os.mkdir(self.regiondir)
                    # read off the section of the new filenames that identifies the instrument and pointing
                    self.unfilt_evls = sorted(glob.glob(self.evlsdir + '/*ImagingEvts.ds'))
                    self.raw_evls = sorted(glob.glob(self.evlsdir + '/*PileupEvts.ds'))
                    self.pointings = [re.search('(E[A-Z]+[0-9]_)*[A-Z][0-9]{3}', f).group(0) for f in
                                      sorted(glob.glob(self.evlsdir + '/*ImagingEvts.ds'))]
                else:
                    raise AssertionError("Event lists have not been reduced for OBSID " + obsdir)
            else:
                self.evls = self.unfilt_evls
            # and the first stange filtered event lists, before background flare removal
            if (len(self.filt_evls) < 1):
                if run_reduction:
                    self.filter_evl()
                    self.filt_evls = sorted(glob.glob(self.evlsdir + '/*filtered.evl'))
                else:
                    print("WARNING: Event lists have not been filtered for OBSID " + obsdir)

            if bkgfilt:
                if (len(self.bkgfilt_evls) < 1):
                    if run_reduction:
                        self.remove_bkg_flares()
                        self.bkgfilt_evls = sorted(glob.glob(self.evlsdir + '/*filtered?gti.evl'))
                        self.evls = self.bkgfilt_evls
                    else:
                        print("WARNING: Background flaring intervals have not been removed for OBSID " + obsdir)

        self.evls = self.bkgfilt_evls if bkgfilt else self.filt_evls


    #-- EPIC pipeline reduction ----------------------------------------------

    def proc_evl(self, pileupcorr=False):
        if self.instrument == 'mos':
            args = ['emproc']
        elif self.instrument == 'pn':
            args = ['epproc']
            if pileupcorr:
                args += ['pileuptempfile=yes', 'runepxrlcorr=yes']

        print("Running %s to produce event lists..." % args[0])

        proc = subprocess.Popen(args, cwd=self.evlsdir, env=self.envvars).wait()

    #-- Event filtering ------------------------------------------------------

    def filter_evl(self, filter_terms=[]):
        #
        # apply the standard event filtering criteria to the event lists
        # and create new, filtered event lists
        #
        print("Filtering event lists...")

        if self.instrument == 'mos':
            filt = ['(PATTERN<=12)','(PI in [200:12000])','(FLAG==0)','#XMMEA_EM'] + filter_terms
        elif self.instrument == self.instrument:
            filt = ['(PATTERN<=4)','(PI in [200:15000])','(FLAG==0)','#XMMEA_EP'] + filter_terms
        expr = '&&'.join(filt)

        for pointing, evl in zip(self.pointings, self.unfilt_evls):
            namearr = [self.obsdir, pointing, 'filtered']
            filt_evl_name = '_'.join(namearr) + '.evl'
            filt_evl = self.evlsdir + '/' + filt_evl_name
            #
            # argument array for Popen
            #
            args = ['evselect',
                    'table='+evl,
                    'expression='+expr,
                    'withfilteredset=yes',
                    'filteredset='+filt_evl]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def remove_bkg_flares(self):
        #
        # remove background flaring intervals from the event list (create a new
        # filtered event list) by creating a GTI filter based on the rate in the
        # highest energy channels
        #
        print("Removing background flaring intervals from event lists...")

        enmin = 10000
        enmax = 12000
        tbin = 10
        if self.instrument == 'mos':
            ratelim = 0.35
        elif self.instrument == self.instrument:
            ratelim = 0.4

        for pointing, evl in zip(self.pointings, self.filt_evls):
            lcnamearr = [self.obsdir, pointing, 'tbin%d' % tbin, 'en%d-%d' % (enmin, enmax)]
            lcname = '_'.join(lcnamearr) + '.lc'
            lc = self.evlsdir + '/' + lcname

            gtinamearr = [self.obsdir, pointing, 'bkgflare_filt']
            gtiname = '_'.join(gtinamearr) + '.gti'
            gti = self.evlsdir + '/' + gtiname

            filtnamearr = [self.obsdir, pointing, 'filtered', 'gti']
            filtname = '_'.join(filtnamearr) + '.evl'
            filt_evl = self.evlsdir + '/' + filtname

            # first, we need a light curve in the high energy channels
            filt = ['(PI in [%d:%d])' % (enmin, enmax),'(PATTERN==0)']
            expr = '&&'.join(filt)
            #
            # argument array for Popen, first to get a light curve
            #
            args = ['evselect',
                    'table='+evl,
                    'energycolumn=PI',
                    'timecolumn=TIME',
                    'expression='+expr,
                    'withrateset=yes',
                    'rateset='+lc,
                    'timebinsize=%d'%tbin,
                    'maketimecolumn=yes',
                    'makeratecolumn=yes']
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

            # now we need to make a GTI
            expr = '(RATE<%g)' % (ratelim)
            args = ['tabgtigen',
                    'table='+lc+':RATE',
                    'gtiset='+gti,
                    'timecolumn=TIME',
                    'expression='+expr]
            proc = subprocess.Popen(args, env=self.envvars).wait()

            # last step, filter the event list with the GTI
            filt = ['gti(%s,TIME)' % gti]
            expr = '&&'.join(filt)
            args = ['evselect',
                    'table='+evl,
                    'expression='+expr,
                    'withfilteredset=yes',
                    'filteredset='+filt_evl]
            proc = subprocess.Popen(args, env=self.envvars).wait()

    #-- Region lookup  -------------------------------------------------------
    def _get_regions(self, region_file='regions.sh'):
        self.regionfiles = []
        for p in self.pointings:
            this_region = glob.glob(self.regiondir + '/' + p + '_' + region_file)
            if(len(this_region) == 1):
                self.regionfiles.append(this_region[0])
            else:
                this_region = glob.glob(self.regiondir + '/' + region_file)
                if(len(this_region) == 1):
                    self.regionfiles.append(this_region[0])
                else:
                    print("WARNING: could not locate region file for " + p)
                    self.regionfiles.append('')
        self.regions = []
        for r in self.regionfiles:
            self.regions.append(self._getregions_sh(r))

    def _getregions_sh(self, regfile):
        #
        # read the source and background regions from variables in a shell script
        #
        if(not os.path.exists(regfile)):
            print("WARNING: Could not read region file. Has it been created yet?")
            return

        f = open(regfile, 'r')
        regionstr = f.read()
        f.close()

        regions = {}

        m = re.search("SRCREGION='(.*?)'", regionstr)
        regions['src'] = m.groups(1)[0]
        m = re.search("BKGREGION='(.*?)'", regionstr)
        regions['bkg'] = m.groups(1)[0]

        return regions

    def find_regions(self, region_size=None, bkg_region=None):
        #
        # try to find the optimal point source extraction region
        # if region_size is not specified, will use the optimal size found by eregionanalyse
        #
        # the regions are saved into the region files and will be automatically loaded next time
        #
        region_re = re.compile(r"SASCIRCLE: (\(X,Y\) in CIRCLE\(.*?\))")
        region_coord_re = re.compile(r"CIRCLE\((.*?),(.*?),(.*?)\)")

        if self.instrument == 'pn':
            src_start = '(DETX, DETY) in CIRCLE(639, -769, 1200)'
        else:
            src_start = '(DETX, DETY) in CIRCLE(0, 0, 1500)'

        if bkg_region is None and self.instrument=='pn':
            bkg_region = '((DETX,DETY) IN circle(-1031.23,2151.2,700))'
        elif bkg_region is None and self.instrument == 'mos':
            bkg_region = '((DETX,DETY) IN circle(-3798,-4540.5,700))'

        self.extract_image(extract_dir=self.regiondir)


        for pointing in self.pointings:
            namearr = [self.obsdir, pointing, self.instrument]
            name = '_'.join(filter(None, namearr))
            imagefile = self.regiondir + '/' + name + '.fits'

            reg_output = subprocess.check_output(['eregionanalyse', 'imageset=%s' % imagefile,
                                     'srcexp=%s' % src_start,
                                     'backexp=%s' % bkg_region], env=self.envvars)

            src_region = region_re.search(str(reg_output)).group(1)

            print("found region: ", src_region)

            if region_size is not None:
                region_coord = region_coord_re.search(src_region)
                x = float(region_coord.group(1))
                y = float(region_coord.group(2))
                r = float(region_coord.group(3))
                src_region = "(X,Y) in CIRCLE(%g,%g,%g)" % (x, y, region_size)

            with open(self.regiondir + "/%s_regions.sh" % pointing, 'w') as f:
                f.write("SRCREGION='%s'\n" % src_region)
                f.write("BKGREGION='%s'\n" % bkg_region)

        self._get_regions()

    def excise_psf(self, radius=25):
        #
        # excise the core of the PSF from the source extraction region
        # specify radius in arcsec
        #
        # region change is temporary and not saved in the region files
        #
        pix_radius = radius / 0.05
        region_re = re.compile(r"circle\(([0-9\.]+),([0-9\.]+),([0-9\.]+)\)")

        for pointing_reg in self.regions:
            pointing_reg['src'] = region_re.sub(r"annulus(\1,\2,%g,\3)" % pix_radius, pointing_reg['src'])

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

    #-- Spectral extraction routines -----------------------------------------

    def extract_spectrum(self, regionkey, extract_dir='', specid='', filter_terms=[], pointings=[]):
        #
        # extract a spectrum
        # uses the region specified by regionkey (usually 'src' or 'bkg'
        # creates the spectrum in extract_dir with identifier spec_id
        # add terms (joined by &&) to the evselect filter expression with filter_terms
        #
        if(not self._ready_to_extract()):
            print("ERROR: Not ready to extract spectrum")
            return

        if(extract_dir == ''):
            extract_dir = self.specdir

        if(not os.path.exists(extract_dir)):
            os.makedirs(extract_dir)

        if(self.instrument == self.instrument):
            specchannelmax = 20479
        elif(self.instrument == 'mos'):
            specchannelmax = 11999

        for pointing, evl, regions in zip(self.pointings, self.evls, self.regions):

            if(len(pointings) > 0 and pointing not in pointings):
                continue

            namearr = [self.obsdir, pointing, regionkey, self.instrument, specid]
            name = '_'.join(filter(None, namearr))

            specfile = extract_dir + '/' + name + '.pha'

            if(os.path.exists(specfile)):
                os.remove(specfile)

            #
            # The SAS filter expression
            #
            expr = '&&'.join([regions[regionkey], '(FLAG==0)'] + filter_terms)
            #
            # Argument array for Popen
            #
            args = ['evselect',
                    'table=' + evl,
                    'expression=' + expr,
                    'withspectrumset=yes',
                    'spectrumset=' + specfile,
                    'energycolumn=PI',
                    'spectralbinsize=5',
                    'withspecranges=yes',
                    'specchannelmin=0',
                    'specchannelmax=' + str(specchannelmax)]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

            #
            # and run backscale
            #
            args = ['backscale',
                    'spectrumset=' + specfile,
                    'badpixlocation=' + evl,
                    ]
            #
            # execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

            return specfile

    def create_rmf(self, regionkey, extract_dir='', specid='', pointings=[]):
        #
        # create the RMF for the spectrum with a specified specid in extract_dir
        #
        if(extract_dir == ''):
            extract_dir = self.specdir

        for pointing in self.pointings:
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            namearr = [self.obsdir, pointing, regionkey, self.instrument, specid]
            name = '_'.join(filter(None, namearr))

            specfile = extract_dir + '/' + name + '.pha'
            rmffile = extract_dir + '/' + name + '.rmf'

            if(not os.path.exists(specfile)):
                print("create_rmf ERROR: Spectrum file not found")
                return

            if(os.path.exists(rmffile)):
                os.remove(rmffile)

            #
            # Argument array for Popen
            #
            args = ['rmfgen',
                    'rmfset=' + rmffile,
                    'spectrumset=' + specfile]

            if self.pileup_corr:
                raw_evl = [evl for evl in self.raw_evls if '%s_%02d' % (pointing,self.ccdnr) in evl][0]
                args += ['correctforpileup=yes', 'raweventfile=%s' % raw_evl]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def create_arf(self, regionkey, extract_dir='', specid='', pointings=[]):
        #
        # create the ARF for the spectrum with a specified specid in extract_dir
        #
        if(extract_dir == ''):
            extract_dir = self.specdir

        for pointing, evl in zip(self.pointings, self.evls):
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            namearr = [self.obsdir, pointing, regionkey, self.instrument, specid]
            name = '_'.join(filter(None, namearr))

            specfile = extract_dir + '/' + name + '.pha'
            rmffile = extract_dir + '/' + name + '.rmf'
            arffile = extract_dir + '/' + name + '.arf'

            if(not os.path.exists(specfile)):
                print("create_arf ERROR: Spectrum file not found")
                return
            if(not os.path.exists(rmffile)):
                print("create_arf ERROR: RMF not found")
                return

            if(os.path.exists(arffile)):
                os.remove(arffile)

            #
            # Argument array for Popen
            #
            args = ['arfgen',
                    'arfset=' + arffile,
                    'spectrumset=' + specfile,
                    'withrmfset=yes',
                    'rmfset=' + rmffile,
                    'badpixlocation=' + evl]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def group_spectrum(self, extract_dir='', specid='', grpmin=20, srckey='src', bkgkey='bkg', pointings=[]):
        #
        # group a spectrum
        # takes the source and background spectra with specified id in extract_dir
        #
        if(extract_dir == ''):
            extract_dir = self.specdir

        for pointing in self.pointings:
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            srcnamearr = [self.obsdir, pointing, srckey, self.instrument, specid]
            srcname = '_'.join(filter(None, srcnamearr))
            bkgnamearr = [self.obsdir, pointing, bkgkey, self.instrument, specid]
            bkgname = '_'.join(filter(None, bkgnamearr))

            specfile = extract_dir + '/' + srcname + '.pha'
            bkgfile = bkgname + '.pha'
            rmffile = srcname + '.rmf'
            arffile = srcname + '.arf'

            grpfile = extract_dir + '/' + srcname + '.grp'

            if(not os.path.exists(specfile)):
                print("group_spectrum ERROR: Source spectrum file not found")
                return
            if(not os.path.exists(extract_dir + '/' + bkgfile)):
                print("group_spectrum WARNING: Background spectrum file not found")
            if(not os.path.exists(extract_dir + '/' + rmffile)):
                print("group_spectrum WARNING: RMF not found")
            if(not os.path.exists(extract_dir + '/' + arffile)):
                print("group_spectrum WARNING: ARF not found")

            if(os.path.exists(grpfile)):
                os.remove(grpfile)

            #
            # command string for grppha
            #
            grpcomm = "reset grouping & group min %d & chkey backfile %s & chkey respfile %s & chkey ancrfile %s & exit" % (
                grpmin, bkgfile, rmffile, arffile)

            #
            # Argument array for Popen
            #
            args = ['grppha',
                    'infile=' + specfile,
                    'outfile=' + grpfile,
                    'comm=' + grpcomm]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def get_spectrum(self, extract_dir='', specid='', grpmin=20, filter_terms=[], pointings=[]):
        #
        # automated extraction of the spectrum
        # extract the source and background spectra, creates the RMF and ARF, then groups the spectrum
        #
        if(extract_dir == ''):
            extract_dir = self.specdir

        self.extract_spectrum('src', extract_dir, specid,
                              filter_terms, pointings=pointings)
        self.extract_spectrum('bkg', extract_dir, specid,
                              filter_terms, pointings=pointings)
        self.create_rmf('src', extract_dir, specid, pointings=pointings)
        self.create_arf('src', extract_dir, specid, pointings=pointings)
        self.group_spectrum(extract_dir, specid, grpmin, pointings=pointings)

    def get_time_spectrum(self, tstart, tend, extract_dir='', specid='', grpmin=20, filter_terms=[], pointings=[], from_start=False)	:
        #
        # extract the spectrum between specfied start and stop time
        #
        if(extract_dir == ''):
            extract_dir = self.specdir
        if(specid == ''):
            specid = "time%g-%g" % (tstart, tend)

        if(from_start):
            tstart = tstart + self.start_time[0]
            if(tend > 0):
                tend = tend + self.start_time[0]

        for pointing, point_start, point_stop in zip(self.pointings, self.start_time, self.stop_time):
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            # only extract if this pointing is during the selected time period
            if((tstart < point_start and tend < point_start) or (tstart > point_stop and tend > point_stop)):
                continue

            if(tstart < point_start):
                tstart = point_start
            if(tend <= 0 or tend > point_stop):
                tend = point_stop

            filt = ["(TIME in [%d:%d])" % (tstart, tend)] + filter_terms
            self.get_spectrum(extract_dir, specid, grpmin, filt, [pointing])

    def get_gti_spectrum(self, gti, extract_dir='', specid='', grpmin=20, filter_terms=[], pointings=[]):
        #
        # extract a spectrum accumulated from time intervals specified in a GTI file
        #
        if(extract_dir == ''):
            extract_dir = self.specdir
        if(specid == ''):
            specid = os.path.splitext(os.path.basename(gti))[0]

        filt = ["gti(%s,TIME)" % gti] + filter_terms
        self.get_spectrum(extract_dir, specid, grpmin, filt, pointings)

    def get_ratecut_spectrum(self, ratemin=0, ratemax=0, timebin=10, extract_dir='', specid='', grpmin=20, gti_dir='', lccorr=True, makelc=True, filter_terms=[], pointings=[]):
        #
        # extract a spectrum accumulated from time intervals specified in a GTI file
        #
        if(extract_dir == ''):
            extract_dir = self.specdir
        if(specid == ''):
            specid = 'rate%g-%g' % (ratemin, ratemax)
        if (gti_dir == ''):
            gti_dir = self.gtidir

        for pointing in self.pointings:
            if (len(pointings) > 0 and pointing not in pointings):
                continue

            self.make_gti(ratemin, ratemax, timebin, lcdir='', lcid='', lckey='src', lccorr=lccorr, makelc=makelc, pointings=[pointing])

            ratestr = '%g-%g' % (ratemin, ratemax)
            gtinamearr = [self.obsdir, pointing, 'src', self.instrument,
                          'tbin' + str(timebin), 'rate' + ratestr, '']
            gtiname = '_'.join(filter(None, gtinamearr))
            gti = gti_dir + '/' + gtiname + '.gti'

            if not os.path.exists(gti):
                continue

            filt = ["gti(%s,TIME)" % gti] + filter_terms
            self.get_spectrum(extract_dir, specid, grpmin, filt, pointings=[pointing])

    #-- Lightcurve extraction routines ---------------------------------------

    def extract_lightcurve(self, timebin, regionkey, extract_dir='', lcid='', filter_terms=[], pointings=[], skip_if_exists=True):
        #
        # extract a spectrum
        # uses the region specified by regionkey (usually 'src' or 'bkg'
        # creates the spectrum in extract_dir with identifier spec_id
        # add terms (joined by &&) to the evselect filter expression with filter_terms
        #
        if(not self._ready_to_extract()):
            print("ERROR: Not ready to extract lightcurve")
            return

        if(extract_dir == ''):
            extract_dir = self.lcdir

        if(not os.path.exists(extract_dir)):
            os.makedirs(extract_dir)

        for pointing, evl, regions in zip(self.pointings, self.evls, self.regions):

            if(len(pointings) > 0 and pointing not in pointing):
                continue

            namearr = [self.obsdir, pointing, regionkey,
                       self.instrument, 'tbin' + str(timebin), lcid]
            name = '_'.join(filter(None, namearr))

            lcfile = extract_dir + '/' + name + '.lc'

            if os.path.exists(lcfile):
                if skip_if_exists:
                    continue
                else:
                    os.remove(lcfile)

            #
            # The SAS filter expression
            #
            expr = '&&'.join([regions[regionkey], '(FLAG==0)'] + filter_terms)
            #
            # Argument array for Popen
            #
            args = ['evselect',
                    'table=' + evl,
                    'expression=' + expr,
                    'withrateset=yes',
                    'rateset=' + lcfile,
                    'timecolumn=TIME',
                    'timebinsize=' + str(timebin),
                    'maketimecolumn=yes',
                    'makeratecolumn=yes']
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def correct_lightcurve(self, timebin, extract_dir='', lcid='', srckey='src', bkgkey='bkg', pointings=[], skip_if_exists=True):

        if(extract_dir == ''):
            extract_dir = self.lcdir

        for pointing, evl in zip(self.pointings, self.evls):

            if(len(pointings) > 0 and pointing not in pointing):
                continue

            srcnamearr = [self.obsdir, pointing, srckey,
                          self.instrument, 'tbin' + str(timebin), lcid]
            srcname = '_'.join(filter(None, srcnamearr))
            bkgnamearr = [self.obsdir, pointing, bkgkey,
                          self.instrument, 'tbin' + str(timebin), lcid]
            bkgname = '_'.join(filter(None, bkgnamearr))

            srclc = extract_dir + '/' + srcname + '.lc'
            bkglc = extract_dir + '/' + bkgname + '.lc'
            corrlc = extract_dir + '/' + srcname + '.lc_corr'

            if(not os.path.exists(srclc)):
                print(
                    "correct_lightcurve ERROR: Could not find source lightcurve to correct")
                return
            if(bkgkey != None and bkgkey != '' and not os.path.exists(bkglc)):
                print(
                    "correct_lightcurve WARNING: Background lightcurve could not be found")

            if os.path.exists(corrlc):
                if skip_if_exists:
                    continue
                else:
                    os.remove(corrlc)

            #
            # Argument array for Popen
            #
            args = ['epiclccorr',
                    'srctslist=' + srclc,
                    'eventlist=' + evl,
                    'outset=' + corrlc,
                    'applyabsolutecorrections=yes']

            if(bkgkey != None and bkgkey != '' and os.path.exists(bkglc)):
                args.append('withbkgset=yes')
                args.append('bkgtslist=' + bkglc)
            else:
                print("correct_lightcurve WARNING: No background subtraction")
                args.append('withbkgset=no')

            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def get_lightcurve(self, timebin, extract_dir='', lcid='', filter_terms=[], pointings=[], skip_if_exists=True):
        #
        # automated extraction of a lightcurve
        # extract the source and background lightcurves then apply correction
        #
        if(extract_dir == ''):
            extract_dir = self.lcdir

        self.extract_lightcurve(timebin, 'src', extract_dir,
                                lcid, filter_terms, pointings=pointings, skip_if_exists=skip_if_exists)
        self.extract_lightcurve(timebin, 'bkg', extract_dir,
                                lcid, filter_terms, pointings=pointings, skip_if_exists=skip_if_exists)
        self.correct_lightcurve(timebin, extract_dir, lcid, pointings=pointings, skip_if_exists=skip_if_exists)

    def get_energy_lightcurve(self, enmin, enmax, timebin, extract_dir='', lcid='', filter_terms=[], pointings=[], skip_if_exists=True):
        #
        # extract the lightcurve in a specified energt range (in eV)
        #
        if(extract_dir == ''):
            extract_dir = self.lcdir
        if(lcid == ''):
            lcid = "en%d-%d" % (enmin, enmax)

        filt = ["(PI in [%d:%d])" % (enmin, enmax)] + filter_terms
        self.get_lightcurve(timebin, extract_dir, lcid,
                            filt, pointings=pointings, skip_if_exists=skip_if_exists)

    def get_time_lightcurve(self, tstart, tend, timebin, extract_dir='', lcid='', filter_terms=[], pointings=[], from_start=False)	:
        #
        # extract the lightcurve between specfied start and stop time
        #
        if(extract_dir == ''):
            extract_dir = self.lcdir
        if(lcid == ''):
            lcid = "time%g-%g" % (tstart, tend)

        if(from_start):
            tstart = tstart + self.start_time[0]
            if(tend > 0):
                tend = tend + self.start_time[0]

        for pointing, point_start, point_stop in zip(self.pointings, self.start_time, self.stop_time):
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            # only extract if this pointing is during the selected time period
            if((tstart < point_start and tend < point_start) or (tstart > point_stop and tend > point_stop)):
                continue

            if(tstart < point_start):
                tstart = point_start
            if(tend <= 0 or tend > point_stop):
                tend = point_stop

            filt = ["(TIME in [%d:%d])" % (tstart, tend)] + filter_terms
            self.get_lightcurve(timebin, extract_dir, lcid, filt, [pointing])

    def get_time_energy_lightcurve(self, tstart, tend, enmin, enmax, timebin, extract_dir='', lcid='', filter_terms=[], pointings=[], from_start=False)	:
        #
        # extract the lightcurve between specfied start and stop time in a specific energy range (in eV)
        #
        if(extract_dir == ''):
            extract_dir = self.lcdir
        if(lcid == ''):
            lcid = "time%g-%g_en%d-%d" % (tstart, tend, enmin, enmax)

        if(from_start):
            tstart = tstart + self.start_time[0]
            if(tend > 0):
                tend = tend + self.start_time[0]

        for pointing, point_start, point_stop in zip(self.pointings, self.start_time, self.stop_time):
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            # only extract if this pointing is during the selected time period
            if((tstart < point_start and tend < point_start) or (tstart > point_stop and tend > point_stop)):
                continue

            if(tstart < point_start):
                tstart = point_start
            if(tend <= 0 or tend > point_stop):
                tend = point_stop

            filt = ["(TIME in [%d:%d])" % (tstart, tend),
                    "(PI in [%d:%d])" % (enmin, enmax)] + filter_terms
            self.get_lightcurve(timebin, extract_dir, lcid, filt, [pointing])

    def make_gti(self, ratemin=0, ratemax=0, timebin=10, gti_dir='', extract_dir='', lcdir='', lcid='', lckey='src', lccorr=False, makelc=True, pointings=[]):

        if(extract_dir == ''):
            extract_dir = self.lcdir
        if(lcdir == ''):
            lcdir = self.lcdir
        if(gti_dir == ''):
            gti_dir = self.gtidir

        if not os.path.exists(gti_dir):
            os.mkdir(gti_dir)

        for pointing, evl in zip(self.pointings, self.evls):

            if(len(pointings) > 0 and pointing not in pointing):
                continue

            if lccorr:
                lcextn = '.lc_corr'
            else:
                lcextn = '.lc'
            lcnamearr = [self.obsdir, pointing, lckey,
                       self.instrument, 'tbin' + str(timebin), lcid]
            lcname = '_'.join(filter(None, lcnamearr))
            lc = lcdir + '/' + lcname + lcextn

            if os.path.exists(lc):
                print("Using existing lightcurve %s" % lc)
            elif makelc:
                print("make_gti: Light curve %s does not exist - creating..." % lc)
                if lccorr:
                    lc = self.get_lightcurve(timebin, lcdir, lcid, pointings=[pointing])
                else:
                    lc = self.extract_lightcurve(timebin, lckey, lcdir, lcid, pointings=[pointing])
            else:
                continue

            ratestr = '%g-%g' % (ratemin, ratemax)
            gtinamearr = [self.obsdir, pointing, lckey, self.instrument,
                          'tbin' + str(timebin), 'rate' + ratestr, lcid]
            gtiname = '_'.join(filter(None, gtinamearr))
            gti = gti_dir + '/' + gtiname + '.gti'

            if(os.path.exists(gti)):
                os.remove(gti)

            expr_terms = []
            if ratemin > 0:
                expr_terms.append('(RATE>=%g)' % ratemin)
            if ratemax > 0:
                expr_terms.append('(RATE<%g)' % ratemax)
            expr = '&&'.join(expr_terms)

            print("DEBUG::lc=%s" % lc )
            print('table=%s:RATE' % lc)

            #
            # Argument array for Popen
            #
            args = ['tabgtigen',
                    'table=%s:RATE' % lc,
                    'gtiset=%s' % gti,
                    'timecolumn=TIME',
                    'expression=' + expr]

            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    #-- Batch light curve generation ------------------------------------
    def batch_energy_lightcurve(self, tbin=10, en0=300, enmax=10000, Nen=10, enbins=None, subdir='energy'):
        if enbins is None:
            enbins = np.round(np.logspace(np.log10(en0), np.log10(enmax), Nen)).astype(int)
        elif enbins == 'lagen':
            enbins = [300,400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000,10000]
        elif enbins == 'lagen_coarse':
            enbins = [300,400,500,600,800,1000,1300,1600,2000,2500,3000,4000,5000,7000]
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
            self.get_energy_lightcurve(s, e, tbin, extract_dir)

    #-- Image extraction routines ---------------------------------------

    def extract_image(self, regionkey=None, extract_dir='', imageid='', filter_terms=[], pointings=[], raimagecenter=None, decimagecenter=None):
        #
        # extract a spectrum
        # uses the region specified by regionkey (usually 'src' or 'bkg'
        # creates the spectrum in extract_dir with identifier spec_id
        # add terms (joined by &&) to the evselect filter expression with filter_terms
        #
        if(not self._ready_to_extract()):
            print("ERROR: Not ready to extract image")
            return

        if(extract_dir == ''):
            extract_dir = self.imagedir

        if(not os.path.exists(extract_dir)):
            os.makedirs(extract_dir)

        for pointing, evl, regions in zip(self.pointings, self.evls, self.regions):

            if(len(pointings) > 0 and pointing not in pointing):
                continue

            namearr = [self.obsdir, pointing, regionkey,
                       self.instrument, imageid]
            name = '_'.join(filter(None, namearr))

            imagefile = extract_dir + '/' + name + '.fits'

            if(os.path.exists(imagefile)):
                os.remove(imagefile)

            #
            # The SAS filter expression
            #
            filter_list = filter_terms
            if regionkey is not None:
                filter_list = filter_list + [regions[regionkey]]

            expr = '&&'.join(filter_list)
            #
            # Argument array for Popen
            #
            args = ['evselect',
                    'table=' + evl,
                    'expression=' + expr,
                    'withimageset=yes',
                    'imageset=' + imagefile,
                    'xcolumn=X',
                    'ycolumn=Y',
                    ]

            if raimagecenter is not None and decimagecenter is not None:
                args += ['raimagecenter='+str(raimagecenter),
                         'decimagecenter=' + str(decimagecenter)]

            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    def extract_time_image(self, tstart, tend, imageid='', filter_terms=[], pointings=[], from_start=False, **kwargs)	:
        #
        # extract the image between specfied start and stop time
        #
        if(imageid == ''):
            imageid = "time%g-%g" % (tstart, tend)

        if(from_start):
            tstart = tstart + self.start_time[0]
            if(tend > 0):
                tend = tend + self.start_time[0]

        for pointing, point_start, point_stop in zip(self.pointings, self.start_time, self.stop_time):
            if(len(pointings) > 0 and pointing not in pointings):
                continue

            # only extract if this pointing is during the selected time period
            if((tstart < point_start and tend < point_start) or (tstart > point_stop and tend > point_stop)):
                continue

            if(tstart < point_start):
                tstart = point_start
            if(tend <= 0 or tend > point_stop):
                tend = point_stop

            filt = ["(TIME in [%d:%d])" % (tstart, tend)] + filter_terms
            self.extract_image(imageid=imageid, filter_terms=filt, pointings=[pointing], **kwargs)


    # -- Event list extraction routines ---------------------------------------

    def extract_event_list(self, regionkey=None, extract_dir='', evl_id='', filter_terms=[], pointings=[], update_exposure=False, orig_evls=None):
        #
        # extract an event list based on some filter criteria
        #
        if (extract_dir == ''):
            extract_dir = self.evlsdir

        if (not os.path.exists(extract_dir)):
            os.makedirs(extract_dir)

        if orig_evls is None:
            orig_evls = self.evls

        for pointing, evl, regions in zip(self.pointings, orig_evls, self.regions):

            if(len(pointings) > 0 and pointing not in pointings):
                continue

            namearr = [self.obsdir, pointing, regionkey, self.instrument, evl_id]
            name = '_'.join(filter(None, namearr))

            evlfile = extract_dir + '/' + name + '.evl'

            if(os.path.exists(evlfile)):
                os.remove(evlfile)

            #
            # The SAS filter expression
            #
            filt = list(filter_terms)
            if regionkey is not None:
                filt.append(regions[regionkey])
            expr = '&&'.join(filt)
            #
            # argument array for Popen
            #
            args = ['evselect',
                    'table=' + evl,
                    'expression=' + expr,
                    'withfilteredset=yes',
                    'filteredset=' + evlfile]
            if update_exposure:
                args.append('updateexposure=yes')
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

    # -- Pile-up test ---------------------------------------------

    def epat_plot(self, epat_dir='', pointings=[], use_gti=True):
        if (epat_dir == ''):
            epat_dir = self.epatdir

        if (not os.path.exists(epat_dir)):
            os.makedirs(epat_dir)

        for pointing, evl, regions in zip(self.pointings, self.evls, self.regions):

            if(len(pointings) > 0 and pointing not in pointings):
                continue

            namearr = [self.obsdir, pointing, self.instrument, 'src', 'epat']
            name = '_'.join(filter(None, namearr))
            plotfile = epat_dir + '/' + name + '.ps'

            namearr_bkg = [self.obsdir, pointing, self.instrument, 'src', 'epat', 'withbkg']
            name_bkg = '_'.join(filter(None, namearr_bkg))
            plotfile_bkg = epat_dir + '/' + name_bkg + '.ps'

            if os.path.exists(plotfile):
                os.remove(plotfile)
            if os.path.exists(plotfile_bkg):
                os.remove(plotfile_bkg)

            filt = []
            if use_gti:
                gtinamearr = [self.obsdir, pointing, 'bkgflare_filt']
                gtiname = '_'.join(gtinamearr) + '.gti'
                gti = self.evlsdir + '/' + gtiname
                filt.append('gti(%s,TIME)' % gti)

                if not os.path.exists(gti):
                    raise AssertionError('Could not find flare removal GTI file')

            self.extract_event_list(regionkey='src', extract_dir=epat_dir, filter_terms=filt,
                                              pointings=[pointing], update_exposure=True, orig_evls=self.unfilt_evls)
            self.extract_event_list(regionkey='bkg', extract_dir=epat_dir, filter_terms=filt,
                                              pointings=[pointing], update_exposure=True, orig_evls=self.unfilt_evls)

            srcnamearr = [self.obsdir, pointing, 'src', self.instrument, '']
            srcname = '_'.join(filter(None, srcnamearr))
            src_evl = epat_dir + '/' + srcname + '.evl'

            bkgnamearr = [self.obsdir, pointing, 'bkg', self.instrument, '']
            bkgname = '_'.join(filter(None, bkgnamearr))
            bkg_evl = epat_dir + '/' + bkgname + '.evl'

            #
            # argument array for Popen
            #
            args = ['epatplot',
                    'set=' + src_evl,
                    'useplotfile=yes',
                    'plotfile=' + plotfile]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

            #
            # and run the test with the background
            #
            args = ['epatplot',
                    'set=' + src_evl,
                    'useplotfile=yes',
                    'plotfile=' + plotfile_bkg,
                    'withbackgroundset=yes',
                    'backgroundset=' + bkg_evl]
            #
            # and execute it
            #
            proc = subprocess.Popen(args, env=self.envvars).wait()

