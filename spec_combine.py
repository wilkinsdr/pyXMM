import subprocess
try:
    import astropy.io.fits as pyfits
except:
    import pyfits
import numpy as np
import glob

from pyxmm.spec_util import *


def combine_spectra(srcfiles='*src*.pha', bkgfiles=None, rmffiles=None, arffiles=None, comb_spec='src_comb.pha', comb_bkg='bkg_comb.pha', comb_rmf='src_comb.rmf', comb_arf='src_comb.arf', comb_grp='src_comb.grp', grpmin=20, exposure_calc='sum'):
    src_list = sorted(glob.glob(srcfiles))

    # if we're not passed search expressions for background, RMF and ARF, try
    # to work them out by replacing bits of the filename
    if bkgfiles is None:
        bkg_list = []
        for src in src_list:
            bkg_list.append( src.replace('src','bkg') )
    else:
        bkg_list = sorted(glob.glob(bkgfiles))
    if rmffiles is None:
        rmf_list = []
        for src in src_list:
            rmf_list.append( src.replace('.pha','.rmf') )
    else:
        rmf_list = sorted(glob.glob(rmffiles))
    if arffiles is None:
        arf_list = []
        for src in src_list:
            arf_list.append( src.replace('.pha','.arf') )
    else:
        arf_list = sorted(glob.glob(arffiles))

    # the weighting of each RMF and ARF in the combined files is equal to the
    # fraction of the total exposure in the corresponding spectrum
    src_exp = []
    for spec in src_list:
        f = pyfits.open(spec)
        hdr = f[1].header
        try:
            src_exp.append( hdr['EXPOSURE'] )
        except:
            src_exp.append(0.)
        f.close()
    sum_srcexp = np.sum(src_exp)
    weight = []
    for exp in src_exp:
        weight.append(exp/sum_srcexp)

    sum_spec(src_list, comb_spec, exposure_calc)
    if comb_bkg is not None:
        sum_spec(bkg_list, comb_bkg, exposure_calc)
    if comb_rmf is not None:
        print("combining rmf")
        combine_rmf(rmf_list, weight, comb_rmf)
    if comb_arf is not None:
        combine_arf(arf_list, weight, comb_arf)
    if comb_grp is not None:
        grppha(comb_grp, comb_spec, comb_bkg, comb_rmf, comb_arf, grpmin)


def sum_spec(spec_list, comb_spec, exposure_calc='sum'):
    exposure_list = []
    backscale_list = []
    for spec in spec_list:
        f = pyfits.open(spec)
        hdr = f[1].header
        try:
            exposure_list.append( hdr['EXPOSURE'] )
            backscale_list.append( hdr['BACKSCAL'] )
        except:
            exposure_list.append(0.)
            backscale_list.append(0.)
        f.close()

    if exposure_calc == 'sum':
        exposure = np.sum(exposure_list)
    if exposure_calc == 'mean':
        exposure = np.mean(exposure_list)

    backscale = np.mean(backscale_list)

    sum_expr = '+'.join(spec_list)

    args = ['mathpha',
            'expr='+sum_expr,
            'errmeth=Gauss',
            'units=COUNTS',
            'outfil='+comb_spec,
            'exposure='+str(exposure),
            'backscal='+str(backscale),
            'areascal=%',
            'ncomments=0',
            'clobber=yes']
    print(' '.join(args))
    subprocess.Popen(args).wait()


def combine_rmf(rmf_list, weights, comb_rmf):
    weightstr_list = []
    for w in weights:
        weightstr_list.append(str(w))

    rmflist = ' '.join(rmf_list)
    weightlist = ' '.join(weightstr_list)

    print("running addrmf")
    print(rmf_list)
    print(weights)

    args = ['addrmf',
            'list='+rmflist,
            'weights='+weightlist,
            'rmffile=!'+comb_rmf,
            'clobber=yes']
    subprocess.Popen(args).wait()

    print(' '.join(args))


def combine_arf(arf_list, weights, comb_arf):
    weightstr_list = []
    for w in weights:
        weightstr_list.append(str(w))

    arflist = ' '.join(arf_list)
    weightlist = ' '.join(weightstr_list)

    args = ['addarf',
            'list='+arflist,
            'weights='+weightlist,
            'out_ARF=!'+comb_arf,
            'clobber=yes']
    subprocess.Popen(args).wait()
