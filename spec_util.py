import subprocess
try:
    import astropy.io.fits as pyfits
except:
    import pyfits
import os


def grppha(grpfile, specfile, bkgfile=None, rmffile=None, arffile=None, grpmin=20):
    #
    # command string for grppha
    #
    grp_commands = ['reset grouping',
                    'group min %d' % grpmin]

    if bkgfile is not None:
        grp_commands += ['chkey backfile %s' % bkgfile]
    if rmffile is not None:
        grp_commands += ['chkey respfile %s' % rmffile]
    if arffile is not None:
        grp_commands += ['chkey ancrfile %s' % arffile]

    grp_commands += ['exit']

    grpcomm = ' & '.join(grp_commands)

    #
    # Argument array for Popen
    #
    args = ['grppha',
            'infile=' + specfile,
            'outfile=!' + grpfile,
            'comm=' + grpcomm]
    #
    # and execute it
    #
    proc = subprocess.Popen(args).wait()


def group_spec(grpfile, specfile, bkgfile=None, rmffile=None, grptype='opt', grpscale=None, usebkg=False):
    #
    # command string for ftgrouppha
    #
    specdir = os.path.dirname(specfile)

    args = ['ftgrouppha',
            os.path.basename(specfile),
            'outfile=%s' % os.path.relpath(grpfile, specdir),
            'grouptype=%s' % grptype,
            'clobber=yes']

    if grpscale is not None:
        args += ['groupscale=%g' % grpscale]

    if rmffile is not None:
        args += ['respfile=%s' % os.path.relpath(rmffile, specdir)]
    elif grptype=='opt':
        raise ValueError('Optimal binning requires the RMF to be provided.')

    if bkgfile is not None:
        args += ['backfile=%s' % bkgfile]

    proc = subprocess.Popen(args, cwd=specdir).wait()


def link_spectra(specfile, bkg=None, rmf=None, arf=None):
    specdir = os.path.dirname(specfile)

    specfits = pyfits.open(specfile, mode='update')

    if bkg is not None:
        specfits['SPECTRUM'].header['BACKFILE'] = os.path.relpath(bkg, specdir)
    if rmf is not None:
        specfits['SPECTRUM'].header['RESPFILE'] = os.path.relpath(rmf, specdir)
    if arf is not None:
        specfits['SPECTRUM'].header['ANCRFILE'] = os.path.relpath(arf, specdir)

    specfits.flush()
    specfits.close()


def update_exposure(specfile, exp):
    srcfits = pyfits.open(specfile, mode='update')
    srcfits[1].header['EXPOSURE'] = exp
    srcfits.flush()
    srcfits.close()

