"""
XRISM Resolve extractor
"""
import glob
import subprocess
from .xselect import Xselect
import os

rsl_grade_def = {'Hp': 0, 'Mp': 1, 'Ms': 2, 'Lp': 3, 'Ls': 4}

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

        self.evls = self.find_evls(filt_level=filt_level)
        if len(self.evls) == 0:
            print('Filter level %s event list not available, falling back to level 1')
            self.evls = self.find_evls(filt_level='')
            self.filt_level = 1

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
            outfile = evl.replace('_cl.evt', '_cl%d.evt' % filt_level)
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
            spec_file = self.specdir + '/' spec_filename

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

    def mkrmf(self, outroot=None, whichrmf='X', grade=['Hp'], pixels='0-11,13-26,28-35', split=True):
        resolist = ','.join([rsl_grade_def[g] for g in grade])

        args = ['rslmkrmf',
                'infile=%s' % evl,
                'outfileroot=' % outroot,
                'regmode=DET',
                'whichrmf=%s' % whichrmf,
                'resolist=%s',
                'regionfile=None',
                'pixlist=%s']

        if split:
            args.append('splitrmf=yes')
            args.append('splitcomb=yes')

        proc = subprocess.Popen(args).wait()

    def get_spectrum(self, grade=['Hp'], pixels='0:11,13:26,28:35', extract_evl=False, suffix=None):
        for evl in self.evls:
            name_arr = ['%srsl' % self.stem]
            if len(self.evls) > 1:
                evl_name = os.basename(evl).split('_')[1]
                name_arr.append(evl_name)
            name_arr.append(''.join(grade))
            if suffix is not None:
                name_arr.append(suffix)

            spec_filename = '_'.join(name_arr) + '.pha'
            spec_file = self.specdir + '/'

            rmf_filename = '_'.join(name_arr) + '_%s' % ''.join(grade) + '.rmf'

            self.extract_spectrum(evl, spec_file, grade, pixels, extract_evl, suffix)

