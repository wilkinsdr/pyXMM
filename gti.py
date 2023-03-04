import os
import astropy.io.fits as pyfits
import numpy as np

class GTI(object):
    #
    # class to create GTI files for XMM-Newton observations
    #

    def __init__(self):
        self.start_time = []
        self.stop_time = []
        self.interval = []

        # the offset is added to each start and stop time
        # set this to enter times relative to observation start time
        # if the GTI requires satellite times (e.g. for XMM)
        self.offset = 0

    def add_row(self, start, stop):
        #
        # add an 'on' interval to the GTI file
        #
        self.start_time.append(start + self.offset)
        self.stop_time.append(stop + self.offset)

    def sort(self):
        #
        # sort the GTI intervals in ascending order of their start time
        # NOTE: there is no check that the intervals do not overlap!
        #
        tab = sorted(zip(self.start_time, self.stop_time))
        self.start_time = [start for (start, stop) in tab]
        self.stop_time = [stop for (start, stop) in tab]

    def calculate_intervals(self):
        #
        # calculate the on-time of each GTI interval
        #
        self.interval = [(stop - start) for (start, stop)
                         in zip(self.start_time, self.stop_time)]

    def on_time(self):
        #
        # calculate the total on time of all intervals in this GTI
        #
        self.calculate_intervals()
        return sum(self.interval)

    def primary_header(self):
        #
        # create the primary header for the GTI FITS file
        #
        hdu = pyfits.PrimaryHDU()
        hdu.header['CREATOR'] = 'pyGTI'
        hdu.header['LONGSTRN'] = 'OGIP 1.0'
        return hdu

    def gti_hdu(self):
        #
        # create the GTI table extension for the FITS file
        #
        self.sort()

        start_col = pyfits.Column(
            name='START', format='1D', array=self.start_time)
        stop_col = pyfits.Column(
            name='STOP', format='1D', array=self.stop_time)
        cols = pyfits.ColDefs([start_col, stop_col])

        hdu = pyfits.BinTableHDU.from_columns(cols)

        # write the appropriate header keywords
        hdu.header['TFORM1'] = ('D', 'data format of field: 8-byte DOUBLE')
        hdu.header['TFORM2'] = ('D', 'data format of field: 8-byte DOUBLE')
        hdu.header['TUNIT1'] = ('s', 'physical unit of field')
        hdu.header['TUNIT2'] = ('s', 'physical unit of field')
        hdu.header['EXTNAME'] = ('STDGTI', 'The name of this table')
        hdu.header['HDUCLASS'] = ('OGIP', 'format conforms to OGIP standard')
        hdu.header['HDUCLAS1'] = ('GTI', 'table contains Good Time Intervals')
        hdu.header['HDUCLAS2'] = (
            'STANDARD', 'standrad Good Time Interval table')
        hdu.header['ONTIME'] = (
            self.on_time(), '[s] sum of all Good Time Intervals')
        hdu.header['TSTART'] = (min(self.start_time),
                                '[s] Lower bound of first GTI')
        hdu.header['TSTOP'] = (max(self.stop_time),
                               '[s] Upper bound of last GTI')
        hdu.header['TIMEUNIT'] = (
            's', 'All times in s unless specified otherwise')
        hdu.header['TIMESYS'] = ('TT', 'XMM time will be TT (Terrestial Time)')
        hdu.header['MJDREF'] = (5.08140000000000E+04,
                                '1998-01-01T00:00:00 (TT) expressed in MJD')
        hdu.header['TIMEREF'] = (
            'LOCAL', 'Reference location of photon arrival times')
        hdu.header['TASSIGN'] = ('SATELLITE', 'Location of time assignment')
        hdu.header['TIMEZERO'] = (0, 'Clock correction (if not zero)')
        hdu.header['CLOCKAPP'] = (True, 'Clock correction applied?')
        hdu.header['MINGTISZ'] = (
            0.00000000000000E+00, '[s] minimum allowed GTI length')
        hdu.header['SUMSHORT'] = (
            0.00000000000000E+00, '[s] sum of all GTIs shorter than MINGTISZ')
        return hdu

    def write(self, filename):
        #
        # write the GTI FITS file
        #
        if(os.path.exists(filename)):
            os.remove(filename)
        prihdr = self.primary_header()
        gtiext = self.gti_hdu()
        hdulist = pyfits.HDUList([prihdr, gtiext])
        hdulist.writeto(filename)

    @staticmethod
    def fromtxt(filename, **kwargs):
        dat = np.genfromtxt(filename, **kwargs)
        gti = GTI()
        for start, stop in dat:
            gti.add_row(start, stop)
        return gti


