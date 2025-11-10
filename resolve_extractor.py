"""
XRISM Resolve extractor
"""
import glob
import subprocess

class ResolveExtractor(object):
    #
    # class to extract data products from XRISM Resolve observations
    #

    def __init__(self, obsdir, filt_level='', run_reduction=False, suffix=None):
        self.obsdir = obsdir

        self.stem = 'xa%s' % obsdir

        self.evlsdir = self.obsdir + '/event_cl' + ('_%s' % suffix if suffix is not None else '')
        self.extractdir = self.obsdir + '/products' + ('_%s' % suffix if suffix is not None else '') + ('_%s' % region_suffix if region_suffix is not None else '')
        self.lcdir = self.obsdir + '/lightcurves' + ('_%s' % suffix if suffix is not None else '') + ('_%s' % region_suffix if region_suffix is not None else '')
        self.regiondir = self.obsdir + '/regions'
        self.gtidir = self.obsdir + '/gti'

    def find_evls(self, filt_level=''):
        self.evls = [evl for evl in sorted(glob.glob(self.evlsdir + '/*_cl%s.evt*' % filt_level)) if 'px5000' not in evl and 'px0000' not in evl]

    def filter(self, rise_time=True, anomalous_Ls=True, frame_events=True, adj_channels=False):
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

