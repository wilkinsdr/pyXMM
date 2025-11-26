"""
XRISM Resolve extractor
"""
import glob
import subprocess
from .xselect import Xselect
import os
import re
from astropy.io import fits

rsl_grade_def = {'Hp': 0, 'Mp': 1, 'Ms': 2, 'Lp': 3, 'Ls': 4}

rsl_pixel_region =['box(4,3,1,1,0)', # pixel 0
                'box(6,3,1,1,0)', # pixel 1
                'box(5,3,1,1,0)', # pixel 2
                'box(6,2,1,1,0)', # pixel 3
                'box(5,2,1,1,0)', # pixel 4
                'box(6,1,1,1,0)', # pixel 5
                'box(5,1,1,1,0)', # pixel 6
                'box(4,2,1,1,0)', # pixel 7
                'box(4,1,1,1,0)', # pixel 8
                'box(1,3,1,1,0)', # pixel 9
                'box(2,3,1,1,0)', # pixel 10
                'box(1,2,1,1,0)', # pixel 11
                '', # no pixel 12
                'box(2,2,1,1,0)', # pixel 13
                'box(2,1,1,1,0)', # pixel 14
                'box(3,2,1,1,0)', # pixel 15
                'box(3,1,1,1,0)', # pixel 16
                'box(3,3,1,1,0)', # pixel 17
                'box(3,4,1,1,0)', # pixel 18
                'box(1,4,1,1,0)', # pixel 19
                'box(2,4,1,1,0)', # pixel 20
                'box(1,5,1,1,0)', # pixel 21
                'box(2,5,1,1,0)', # pixel 22
                'box(1,6,1,1,0)', # pixel 23
                'box(2,6,1,1,0)', # pixel 24
                'box(3,5,1,1,0)', # pixel 25
                'box(3,6,1,1,0)', # pixel 26
                'box(6,4,1,1,0)', # pixel 27
                'box(5,4,1,1,0)', # pixel 28
                'box(6,5,1,1,0)', # pixel 29
                'box(6,6,1,1,0)', # pixel 30
                'box(5,5,1,1,0)', # pixel 31
                'box(5,6,1,1,0)', # pixel 32
                'box(4,5,1,1,0)', # pixel 33
                'box(4,6,1,1,0)', # pixel 34
                'box(4,4,1,1,0)' # pixel 35
                ]

class ResolveExtractor(object):
    #
    # class to extract data products from XRISM Resolve observations
    #

    def __init__(self, obsdir, filt_level=2, run_reduction=False, suffix=None):
        self.obsdir = obsdir
        self.rsldir = obsdir + '/resolve'

        self.stem = 'xa%s' % obsdir

        self.suffix = suffix

        self.evlsdir = self.rsldir + '/event_cl' + ('_%s' % suffix if suffix is not None else '')
        self.specdir = self.rsldir + '/spectra' + ('_%s' % suffix if suffix is not None else '')
        self.lcdir = self.rsldir + '/lightcurves' + ('_%s' % suffix if suffix is not None else '')
        self.regiondir = self.rsldir + '/regions'
        self.gtidir = self.rsldir + '/gti'
        self.ufdir = self.rsldir + '/event_uf'
        self.scratchdir = self.rsldir + '/scratch'

        if not os.path.exists(self.scratchdir):
            os.mkdir(self.scratchdir)

        self.ehk = glob.glob(self.obsdir + '/auxil/%s.ehk*' % self.stem)[0]

        self.evls = self.find_evls(filt_level=filt_level)
        if len(self.evls) == 0:
            print('Filter level %d event list not available, falling back to level 1' % filt_level)
            self.evls = self.find_evls(filt_level='')
            self.filt_level = 1

        with fits.open(self.evls[0]) as f:
            self.ra_nom = f[0].header['RA_NOM']
            self.dec_nom = f[0].header['DEC_NOM']

    def find_evls(self, filt_level=''):
        return [evl for evl in sorted(glob.glob(self.evlsdir + '/%srsl_p0px????_cl%s.evt*' % (self.stem, filt_level))) if 'px5000' not in evl and 'px0000' not in evl]

    def filter(self, rise_time=True, anomalous_Ls=True, frame_events=True, adj_channels=False, filt_level=2):
        criteria = ['(PI>=600)']
        if rise_time:
            if anomalous_Ls:
                # filter on rise time and exclude anomalous Ls events
                criteria.append('(((((RISE_TIME+0.00075*DERIV_MAX)>46)&&((RISE_TIME+0.00075*DERIV_MAX)<58))&&ITYPE<4))')
            else:
                # filter on rise time but keep anomalous Ls events
                criteria.append('(((((RISE_TIME+0.00075*DERIV_MAX)>46)&&((RISE_TIME+0.00075*DERIV_MAX)<58))&&ITYPE<4)||(ITYPE==4))')
        if frame_events:
            # STATUS[4] indicates particle events absorbed into silicon frame
            criteria.append('STATUS[4]==b0')
        if adj_channels:
            # STATUS[13] indicates events in electrically adjacent channels - use this for bright sources
            criteria.append('STATUS[13]==b0')

        filt_str = '&&'.join(criteria)

        for evl in self.evls:
            outfile = evl.replace('_cl.evt', '_cl%d.evt' % filt_level).replace('.gz','')
            args = ['ftcopy',
                    'infile=%s[EVENTS][%s]' % (evl, filt_str),
                    'outfile=%s' % outfile,
                    'copyall=yes',
                    'clobber=yes',
                    'history=yes'
                    ]

            proc = subprocess.Popen(args).wait()

        self.evls = self.find_evls(filt_level=filt_level)

    def extract_spectrum(self, evl=None, spec_file=None, grade=['Hp'], pixels='0:11,13:26,28:35', extract_evl=False, suffix=None):
        if isinstance(grade, str):
            grade = [grade]
        grade_sel = [rsl_grade_def[g] for g in grade]
        grade_filt = '%d:%d' % (min(grade_sel), max(grade_sel))

        if evl is None:
            evl = self.evls[0]
        if spec_file is None:
            name_arr = ['%srsl' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            name_arr.append(''.join(grade))
            if suffix is not None:
                name_arr.append(suffix)

            spec_filename = '_'.join(name_arr) + '.pha'
            spec_file = self.specdir + '/' + spec_filename

        with Xselect() as xsl:
            xsl.read_event(evl)
            if extract_evl:
                # if we're filtering events by time, we need to extract an event list with the time filter applied
                # for rslmkrmf
                xsl.command('extract events')
                xsl.command('save events %s' % extract_evl)
            if pixels is not None:
                xsl.command('filter column "PIXEL=%s"' % pixels)
            xsl.command('filter GRADE %s' % grade_filt)
            xsl.command('extract spectrum')
            xsl.command('save spectrum %s' % spec_file)

    def make_rmf(self, outroot=None, evl=None, whichrmf='X', grade=['Hp'], pixels='0-11,13-26,28-35', split=True):
        resolist = ','.join(['%d' % rsl_grade_def[g] for g in grade])

        if evl is None:
            evl = self.evls[0]

        name_arr = ['%srsl' % self.stem]
        if len(self.evls) > 1:
            evl_name = os.path.basename(evl).split('_')[1]
            name_arr.append(evl_name)
        if outroot is None:
            outroot = self.specdir + '/' + '_'.join(name_arr) + '_%s' % ''.join(grade) + '_%s' % whichrmf

        args = ['rslmkrmf',
                'infile=%s' % evl,
                'outfileroot=%s' % outroot,
                'regmode=DET',
                'whichrmf=%s' % whichrmf,
                'resolist=%s' % resolist,
                'regionfile=None',
                'pixlist=%s' % pixels]

        if split:
            args.append('splitrmf=yes')
            args.append('splitcomb=yes')

        proc = subprocess.Popen(['punlearn', 'rslmkrmf']).wait()
        proc = subprocess.Popen(args).wait()

    def exposure_map(self, outfile=None, evl=None):
        if evl is None:
            evl = self.evls[0]

        name_arr = ['%srsl' % self.stem]
        if len(self.evls) > 1:
            evl_name = os.path.basename(evl).split('_')[1]
            name_arr.append(evl_name)
        if outfile is None:
            outfile = self.specdir + '/' + '_'.join(name_arr) + '.expo'

        filtstr = re.search('(px[0-9]+)', os.path.basename(evl)).group(1)
        pixgtifile = glob.glob(self.ufdir + '/%srsl_%s_exp.gti*' % (self.stem, filtstr))[0]

        args = ['xaexpmap',
                'ehkfile=%s' % self.ehk,
                'gtifile=%s' % evl,
                'instrume=RESOLVE',
                'badimgfile=NONE',
                'pixgtifile=%s' % pixgtifile,
                'outfile=%s' % outfile,
                'outmaptype=EXPOSURE',
                'delta=20.0',
                'numphi=1']

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'xaexpmap']).wait()
        proc = subprocess.Popen(args).wait()

    def make_region_file(self, outfile=None, pixels='0:11,13:26,28:35'):
        if outfile is None:
            outfile = self.scratchdir + '/' + self.stem + '_DET.reg'

        with open(outfile, 'w') as f:
            for pixrange in pixels.split(','):
                a, b = pixrange.split(':')
                for box in rsl_pixel_region[int(a):int(b)+1]:
                    f.write(box + '\n')

    def make_arf_pointsource(self, outfile=None, xrtevtfile=None, ra=None, dec=None, expomap=None, regionfile=None, rmffile=None, erange='1.5 18.0 0 0', numphoton=300000):
        if ra is None:
            ra = self.ra_nom
        if dec is None:
            dec = self.dec_nom

        if xrtevtfile is None:
            xrtevtfile = self.scratchdir + '/' + os.path.basename(outfile).replace('.arf', '_raytrace_pt.evt')

        if os.path.exists(outfile):
            os.remove(outfile)
        if os.path.exists(xrtevtfile):
            os.remove(xrtevtfile)

        if expomap is None or not os.path.exists(expomap):
            raise AssertionError('Exposure map does not exist')
        if regionfile is None or not os.path.exists(regionfile):
            raise AssertionError('Region file does not exist')
        if rmffile is None or not os.path.exists(rmffile):
            raise AssertionError('RMF does not exist')

        args = ['xaarfgen',
                'xrtevtfile=%s' % xrtevtfile,
                'source_ra=%0.5f' % ra,
                'source_dec=%0.5f' % dec,
                'telescop=XRISM',
                'instrume=RESOLVE',
                'emapfile=%s' % expomap,
                'regmode=DET',
                'regionfile=%s' % regionfile,
                'sourcetype=POINT',
                'rmffile=%s' % rmffile,
                'erange="%s"' % erange,
                'outfile=%s' % outfile,
                'numphoton=%d' % numphoton,
                'qefile=CALDB',
                'contamifile=CALDB',
                'gatevalvefile=CALDB',
                'onaxisffile=CALDB',
                'onaxiscfile=CALDB',
                'mirrorfile=CALDB',
                'obstructfile=CALDB',
                'frontreffile=CALDB',
                'backreffile=CALDB',
                'pcolreffile=CALDB',
                'scatterfile=CALDB',
                'imgfile=NONE']

        print(' '.join(args))

        proc = subprocess.Popen(['punlearn', 'xaarfgen']).wait()
        proc = subprocess.Popen(args).wait()

    def get_spectrum(self, grade=['Hp'], pixels='0:11,13:26,28:35', extract_evl=False, whichrmf='X', split_rmf=True, ra=None, dec=None, suffix=None, extract_spectrum=True, make_rmf=True, make_arf=True):
        for evl in self.evls:
            name_arr = ['%srsl' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.path.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            if suffix is not None:
                name_arr.append(suffix)

            print(name_arr)

            spec_filename = '_'.join(name_arr) + '_%s' % ''.join(grade) + '.pha'
            spec_file = self.specdir + '/' + spec_filename

            print(spec_file)

            if not os.path.exists(self.specdir):
                os.mkdir(self.specdir)

            rmf_root = self.specdir + '/' + '_'.join(name_arr) + '_%s' % ''.join(grade) + '_%s' % whichrmf
            rmf_file = rmf_root + '_comb.rmf' if split_rmf else rmf_root + '.rmf'
            expomap_file = self.scratchdir + '/' + '_'.join(name_arr) + '.expo'
            region_file = self.scratchdir + '/' + self.stem + '_DET.reg'
            arf_file = self.specdir + '/' + '_'.join(name_arr) + '_src.arf'

            if extract_spectrum:
                self.extract_spectrum(evl, spec_file, grade, pixels, extract_evl, suffix)
            if make_rmf:
                self.make_rmf(rmf_root, whichrmf, grade, pixels.replace(':','-'), split_rmf)
            if make_arf:
                self.exposure_map(expomap_file, evl)
                self.make_region_file(region_file, pixels)
                self.make_arf_pointsource(arf_file, ra=ra, dec=dec, expomap=expomap_file, regionfile=region_file, rmffile=rmf_file)

