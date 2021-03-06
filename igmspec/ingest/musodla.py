""" Module to ingest MUSoDLA survey

Jorgensen et al. 2013
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os, json, glob, imp

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, Column
from astropy.time import Time
from astropy import units as u

from linetools import utils as ltu
from linetools.spectra import io as lsio

from specdb.build.utils import chk_meta
from specdb.build.utils import init_data

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Ingest supplied meta table
    Returns
    -------
    meta : Table
    """
    # Cut down to unique QSOs
    musodla_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/MUSoDLA/datatab_v2.dat', format='ascii')
    mdict = {1:'MagE', 2:'XSHOOTER', 3:'UVES', 4:'HIRES'}
    mRdict = {'MagE':71., 'XSHOOTER':59., 'UVES':7., 'HIRES':7.}
    gdict = {'MagE':'N/A', 'XSHOOTER':'ALL', 'UVES':'BOTH', 'HIRES':'RED'}
    tdict = {'MagE':'Magellan', 'XSHOOTER':'VLT', 'UVES':'VLT', 'HIRES':'Keck-I'}
    spfdict = {'MagE':'MagE', 'XSHOOTER':'XShooter', 'UVES':'UVES', 'HIRES':'HIRES'}
    coords = []
    zems = []
    instrs = []
    names = []
    dinfos = []
    dates = []
    sfiles = []
    Rs = []
    gratings = []
    telescopes = []
    for row in musodla_meta:
        # Coord
        coord = ltu.radec_to_coord((row['RA(J2000)'],row['DEC(J2000)']))
        # Instruments
        insts = row['I'].split(',')
        dinfo = row['date_info'].split(';')
        for jj,inst in enumerate(insts):
            coords.append(coord)
            zems.append(row['z_em'])
            instr = mdict[int(inst)]
            instrs.append(instr)
            Rs.append(3e5/mRdict[instr])
            gratings.append(gdict[instr])
            telescopes.append(tdict[instr])
            names.append(row['QSOname'])
            sfiles.append(row['QSOname']+'_{:s}.ascii'.format(spfdict[instr]))
            # Date
            assert dinfo[jj][0] == inst
            dinfos.append(str(dinfo[jj]))
            dsplit = dinfo[jj].split(',')
            dates.append(dsplit[-1])
            # Special case
            if 'J1201+0116' in row['QSOname']:
                sfiles[-1] = 'J1201+0116_HIRES_1.3kms.ascii'
                sfiles.append('J1201+0116_HIRES_2.6kms.ascii')
                instrs.append(instr)
                Rs.append(3e5/mRdict[instr])
                gratings.append(gdict[instr])
                telescopes.append(tdict[instr])
                names.append(row['QSOname'])
                dinfos.append(str(dinfo[jj]))
                dates.append(dsplit[-1])
                coords.append(coord)
                zems.append(row['z_em'])
    # Generate
    meta = Table()
    meta['RA_GROUP'] = [coord.ra.deg for coord in coords]
    meta['DEC_GROUP'] = [coord.dec.deg for coord in coords]
    meta['NAME'] = names
    meta['zem_GROUP'] = zems
    meta['INSTR'] = instrs
    meta['DISPERSER'] = gratings
    meta['DATE-INFO'] = dinfos
    meta['SPEC_FILE'] = sfiles
    meta['TELESCOPE'] = telescopes
    meta['R'] = Rs
    t = Time(dates, out_subfmt='date')  # Fixes to YYYY-MM-DD
    meta.add_column(Column(t.iso, name='DATE-OBS'))
    #
    meta['sig_zem'] = 0.
    meta['flag_zem'] = str('SDSS')
    meta['STYPE'] = str('QSO')
    # Check
    assert chk_meta(meta, chk_cat_only=True)
    # Return
    return meta


def hdf5_adddata(hdf, sname, musodla_meta, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append MUSoDLA data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGM_ID values in mainDB
    sname : str
      Survey name
    chk_meta_only : bool, optional
      Only check meta file;  will not write
    mk_test_file : bool, optional
      Generate the debug test file for Travis??

    Returns
    -------

    """
    from specdb.build.utils import chk_meta
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    hdlls_grp = hdf.create_group(sname)
    # Checks
    if sname != 'MUSoDLA':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    max_npix = 230000  # Just needs to be large enough
    data = init_data(max_npix, include_co=False)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    nspec = len(musodla_meta)
    spec_set.resize((nspec,))
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    # Loop
    for jj,row in enumerate(musodla_meta):
        kk = jj
        # Extract
        f = os.getenv('RAW_IGMSPEC')+'/MUSoDLA/data/'+row['SPEC_FILE']
        try:
            spec = lsio.readspec(f, masking='edges')
        except:
            pdb.set_trace()
        # Parse name
        fname = f.split('/')[-1]
        # npix
        head = spec.header
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        if 'MagE' in f:
            if fname in ['J2122-0014_MagE.ascii','J0011+1446_MagE.ascii']:
                data['sig'][0][:npix] = spec.sig.value  # Special cases..
            else:
                data['sig'][0][:npix] = 1./np.sqrt(spec.sig.value)  # IVAR
        else:
            data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        # Meta
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[kk] = data

    # Add columns
    nmeta = len(musodla_meta)
    musodla_meta.add_column(Column([2000.]*nmeta, name='EPOCH'))
    musodla_meta.add_column(Column(npixlist, name='NPIX'))
    musodla_meta.add_column(Column(wvminlist, name='WV_MIN'))
    musodla_meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    musodla_meta.add_column(Column(np.arange(nmeta,dtype=int),name='GROUP_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(musodla_meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = musodla_meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2013MNRAS.435..482J',
                 bib='regina+13'),
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


def add_ssa(hdf, dset):
    """  Add SSA info to meta dataset
    Parameters
    ----------
    hdf
    dset : str
    """
    from specdb.ssa import default_fields
    Title = '{:s}: The Magellan uniform survey of damped Lya systems'.format(dset)
    ssa_dict = default_fields(Title, flux='normalized')
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))


def chk_mage_flux():
    """ At least one MagE file has a mix of fluxed and normalized spectra
    And a 'normal' error array
    """
    mage_files = glob.glob(os.getenv('RAW_IGMSPEC')+'/MUSoDLA/data/*MagE.ascii')
    for mfile in mage_files:
        spec = lsio.readspec(mfile)
        nhigh = np.sum(spec.flux.value > 1e10)
        print("File={:s} with nhigh={:d}".format(mfile, nhigh))
        if nhigh > 1:
            pdb.set_trace()
