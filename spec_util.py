import subprocess
try:
    import astropy.io.fits as pyfits
except:
    import pyfits


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


def update_exposure(specfile, exp):
    srcfits = pyfits.open(specfile, mode='update')
    srcfits[1].header['EXPOSURE'] = exp
    srcfits.flush()
    srcfits.close()
