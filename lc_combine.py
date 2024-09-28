import re
import glob
import os

try:
    from pylag import *
except ImportError:
    raise ImportError("pylag is required for light curve manipulation")

def sum_mos_pn_energy_lightcurves(pn_lcpath='.', mos_lcpath='.', output_path='.'):
    pn_lc = sorted(glob.glob(pn_lcpath + '/*src_pn*.lc_corr'))

    # regex to pull information out of light curve filenames
    lc_filename_re = re.compile('([0-9]{10})_([A-Z][0-9]{3})_src_(pn|mos)_tbin([0-9]+)_(en[0-9\-]+)\.lc_corr')
    # and some lambdas to extract useful information
    lc_obsid = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(1)
    lc_expid = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(2)
    lc_instr = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(3)
    lc_tbin = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(4)
    lc_en = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(5)

    # find out which obsids, time binnings and energy bands we have
    obsids = set([lc_obsid(l) for l in pn_lc])
    tbins = set([lc_tbin(l) for l in pn_lc])
    enbands = set([lc_en(l) for l in pn_lc])

    # and iterate through them to produce the combibed light curves
    for obsid in obsids:
        for tbin in tbins:
            dt = int(tbin)
            for enband in enbands:
                lcl_pn = get_lclist(pn_lcpath + '/%s_*src_pn_tbin%s_%s.lc_corr' % (obsid, tbin, enband))
                try:
                    lcl_mos1 = get_lclist(mos_lcpath + '/%s_EMOS1_*src_mos_tbin%s_%s.lc_corr' % (obsid, tbin, enband))
                except:
                    lcl_mos1 = []
                try:
                    lcl_mos2 = get_lclist(mos_lcpath + '/%s_EMOS2_*src_mos_tbin%s_%s.lc_corr' % (obsid, tbin, enband))
                except:
                    lcl_mos2 = []

                # if we have more than one light curve from each instrument within an obsid
                # we concatenate them into a single light curve, interpolating over the gaps
                # otherwise, if we just have one light curve, use it
                if len(lcl_pn) == 1:
                    pn_lc = lcl_pn[0].interp_gaps()
                elif len(lcl_pn) > 1:
                    pn_lc = LightCurve().concatenate(lcl_pn).remove_gaps().fill_time(dt=dt).interp_gaps()

                if len(lcl_mos1) == 1:
                    mos1_lc = lcl_mos1[0].interp_gaps()
                elif len(lcl_mos1) > 1:
                    mos1_lc = LightCurve().concatenate(lcl_mos1).remove_gaps().fill_time(dt=dt).interp_gaps()
                else:
                    print('MOS1 missing for', obsid, tbin, enband)
                    mos1_lc = None

                if len(lcl_mos2) == 1:
                    mos2_lc = lcl_mos2[0].interp_gaps()
                elif len(lcl_mos2) > 1:
                    mos2_lc = LightCurve().concatenate(lcl_mos2).remove_gaps().fill_time(dt=dt).interp_gaps()
                else:
                    print('MOS2 missing for', obsid, tbin, enband)
                    mos2_lc = None

                # if we have both MOS1 and MOS2, add the light curves together
                if mos1_lc is not None and mos2_lc is not None:
                    try:
                        mos1_lc, mos2_lc = extract_sim_lightcurves(mos1_lc, mos2_lc)
                        mos_sum_lc = mos1_lc + mos2_lc
                    except:
                        mos_sum_lc = None
                elif mos1_lc is not None:
                    mos_sum_lc = mos1_lc
                elif mos2_lc is not None:
                    mos_sum_lc = mos2_lc
                else:
                    mos_sum_lc = None

                # if we have a MOS light curve, add it to the pn
                if mos_sum_lc is not None:
                    mos_sum_lc, pn_lc = extract_sim_lightcurves(mos_sum_lc, pn_lc)
                    sum_lc = mos_sum_lc + pn_lc
                else:
                    sum_lc = pn_lc

                outfile = output_path + '/%s_src_mospn_tbin%s_%s.lc' % (obsid, tbin, enband)
                sum_lc.write_fits(outfile)


def sum_mos_pn_lightcurves(pn_lcpath='.', mos_lcpath='.', output_path='.'):
    pn_lc = sorted(glob.glob(pn_lcpath + '/*src_pn*.lc_corr'))

    # regex to pull information out of light curve filenames
    lc_filename_re = re.compile('([0-9]{10})_([A-Z][0-9]{3})_src_(pn|mos)_tbin([0-9]+)\.lc_corr')
    # and some lambdas to extract useful information
    lc_obsid = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(1)
    lc_expid = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(2)
    lc_instr = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(3)
    lc_tbin = lambda lc: lc_filename_re.match(os.path.basename(lc)).group(4)

    # find out which obsids, time binnings and energy bands we have
    obsids = set([lc_obsid(l) for l in pn_lc])
    tbins = set([lc_tbin(l) for l in pn_lc])

    # and iterate through them to produce the combibed light curves
    for obsid in obsids:
        for tbin in tbins:
            dt = int(tbin)
            lcl_pn = get_lclist(pn_lcpath + '/%s_*src_pn_tbin%s.lc_corr' % (obsid, tbin))
            try:
                lcl_mos1 = get_lclist(mos_lcpath + '/%s_EMOS1_*src_mos_tbin%s.lc_corr' % (obsid, tbin))
            except:
                lcl_mos1 = []
            try:
                lcl_mos2 = get_lclist(mos_lcpath + '/%s_EMOS2_*src_mos_tbin%s.lc_corr' % (obsid, tbin))
            except:
                lcl_mos2 = []

            # if we have more than one light curve from each instrument within an obsid
            # we concatenate them into a single light curve, interpolating over the gaps
            # otherwise, if we just have one light curve, use it
            if len(lcl_pn) == 1:
                pn_lc = lcl_pn[0].interp_gaps()
            elif len(lcl_pn) > 1:
                pn_lc = LightCurve().concatenate(lcl_pn).remove_gaps().fill_time(dt=dt).interp_gaps()

            if len(lcl_mos1) == 1:
                mos1_lc = lcl_mos1[0].interp_gaps()
            elif len(lcl_mos1) > 1:
                mos1_lc = LightCurve().concatenate(lcl_mos1).remove_gaps().fill_time(dt=dt).interp_gaps()
            else:
                print('MOS1 missing for', obsid, tbin, enband)
                mos1_lc = None

            if len(lcl_mos2) == 1:
                mos2_lc = lcl_mos2[0].interp_gaps()
            elif len(lcl_mos2) > 1:
                mos2_lc = LightCurve().concatenate(lcl_mos2).remove_gaps().fill_time(dt=dt).interp_gaps()
            else:
                print('MOS2 missing for', obsid, tbin, enband)
                mos2_lc = None

            # if we have both MOS1 and MOS2, add the light curves together
            if mos1_lc is not None and mos2_lc is not None:
                mos1_lc, mos2_lc = extract_sim_lightcurves(mos1_lc, mos2_lc)
                mos_sum_lc = mos1_lc + mos2_lc
            elif mos1_lc is not None:
                mos_sum_lc = mos1_lc
            elif mos2_lc is not None:
                mos_sum_lc = mos2_lc
            else:
                mos_sum_lc = None

            # if we have a MOS light curve, add it to the pn
            if mos_sum_lc is not None:
                mos_sum_lc, pn_lc = extract_sim_lightcurves(mos_sum_lc, pn_lc)
                sum_lc = mos_sum_lc + pn_lc
            else:
                sum_lc = pn_lc

            outfile = output_path + '/%s_src_mospn_tbin%s.lc' % (obsid, tbin)
            sum_lc.write_fits(outfile)
