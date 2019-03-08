""" Module to ingest UVES SQUAD DR1 data

Murphy et al. 2018
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os, json, glob, imp
import datetime

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io import fits
from astropy.time import Time

from linetools import utils as ltu
from linetools.spectra import io as lsio

from specdb.build.utils import chk_meta
from specdb.build.utils import init_data
from specdb.build.utils import set_resolution

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Generates the meta data needed for the IGMSpec build

    Returns
    -------
    meta : Table
    spec_files : list
      List of spec_file names
    """
    # Load the summary table
    path = os.path.join(os.getenv('RAW_IGMSPEC'), 'UVES_SQUAD_DR1')
    squad_meta = Table.read(os.path.join(path, 'DR1_quasars_master.csv'))
    # Limit to those with spectra
    keep = np.array([True]*len(squad_meta))
    for ii in [3,4,5]:
        keep = keep & np.invert(squad_meta['Spec_status'] == str(ii))
    squad_meta = squad_meta[keep]
    # Cut down to unique QSOs
    spec_files = []
    for row in squad_meta:
        # SPEC_FILE
        fname = row['Name_Adopt']+'.fits'
        spec_files.append(fname)
    nqso = len(squad_meta)
    # Coord me
    coord = SkyCoord(ra=squad_meta['RA_Adopt'], dec=squad_meta['Dec_Adopt'], unit=(units.hour, units.deg))
    squad_meta['RA_GROUP'] = coord.ra.value
    squad_meta['DEC_GROUP'] = coord.dec.value
    # Rename a few
    squad_meta.rename_column('zem_Adopt', 'zem_GROUP')
    squad_meta['sig_zem'] = 0.
    #
    zsource = [None]*nqso
    gd_sdss = squad_meta['zem_SDSS'] > 0.
    for ii in np.where(gd_sdss)[0]:
        zsource[ii] = str('SDSS')
    gd_NED = squad_meta['zem_SDSS'].mask & (squad_meta['zem_NED'] > 0.)
    for ii in np.where(gd_NED)[0]:
        zsource[ii] = str('NED')
    gd_SIMBAD = squad_meta['zem_SDSS'].mask & squad_meta['zem_NED'].mask
    for ii in np.where(gd_SIMBAD)[0]:
        zsource[ii] = str('SIMBAD')
    squad_meta['flag_zem'] = zsource
    #
    squad_meta.rename_column('WavStart', 'WV_MIN')
    squad_meta.rename_column('WavEnd', 'WV_MAX')
    squad_meta['INSTR'] = 'UVES'
    squad_meta['TELESCOPE'] = 'VLT'
    squad_meta['DISPERSER'] = 'BOTH'
    squad_meta['EPOCH'] = 2000.
    squad_meta['STYPE'] = str('QSO')
    squad_meta['SPEC_FILE'] = spec_files
    # Check
    assert chk_meta(squad_meta, chk_cat_only=True)
    return squad_meta


def hdf5_adddata(hdf, sname, squad_meta, debug=False, chk_meta_only=False):
    """ Append UVES SQUAD data to the h5 file

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
    from specdb import defs
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    hdlls_grp = hdf.create_group(sname)
    # Load up
    Rdicts = defs.get_res_dicts()
    # Checks
    if sname != 'SQUAD_DR1':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    #if mk_test_file:
    #    hdla100_full = hdlls_full[0:3]
    max_npix = 300000  # Just needs to be large enough
    data = init_data(max_npix, include_co=True)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    nspec = len(squad_meta)
    spec_set.resize((nspec,))
    Rlist = []
    dateobslist = []
    npixlist = []
    gratinglist = []
    # Loop
    for jj,row in enumerate(squad_meta):
        kk = jj
        # Extract
        f = os.path.join(os.getenv('RAW_IGMSPEC'),'UVES_SQUAD_DR1', 'spectra',
                         row['SPEC_FILE'][:-5], row['SPEC_FILE'])
        spec = lsio.readspec(f)
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
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        data['co'][0][:npix] = spec.co.value
        # Meta
        hdu = fits.open(f)
        # Date
        exp_tbl = hdu[3].data
        items = []
        for item in exp_tbl['UTDate']:
            items.append(item.strip())
        times = Time(items, format='iso')
        dateobslist.append(times.min().value[0:10])
        # R
        comb_tbl = hdu[1].data
        Rlist.append(int(np.mean(comb_tbl['NomResolPower'][comb_tbl['NomResolPower'] > 0.])))
        # npix
        npixlist.append(npix)
        # Done
        spec_set[kk] = data

    # Add columns
    squad_meta['GROUP_ID'] = np.arange(nspec, dtype=int)
    squad_meta['R'] = Rlist
    squad_meta['NPIX'] = npixlist
    squad_meta['DATE-OBS'] = dateobslist

    # Add HDLLS meta to hdf5
    if chk_meta(squad_meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = squad_meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2019MNRAS.482.3458M',
                 bib='murphy+19'),
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
    Title = '{:s}: The Keck/HIRES Survey of 100 Damped Lya Systems'.format(dset)
    ssa_dict = default_fields(Title, flux='normalized')
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))
