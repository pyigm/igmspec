""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, warnings

import h5py
import json
import datetime
import pdb
from collections import OrderedDict

from specdb import defs
from specdb.build import utils as sdbbu

import igmspec
from igmspec.ingest import boss, hdlls, kodiaq, ggg, sdss, hst_z2, myers, twodf, xq100
from igmspec.ingest import hdla100
from igmspec.ingest import esidla
from igmspec.ingest import cos_halos
from igmspec.ingest import hst_qso
from igmspec.ingest import hst_cooksey as hst_c
from igmspec.ingest import cos_dwarfs
from igmspec.ingest import musodla
from igmspec.ingest import uves_dall
from igmspec.ingest import boss_dr14
from igmspec.ingest import esi_z6
from igmspec.ingest import kodiaq_two
from igmspec.ingest import uves_squad

from astropy.table import Table


def ver01(test=False, clobber=False, publisher='J.X. Prochaska', **kwargs):
    """ Build version 1.0

    Parameters
    ----------
    test : bool, optional
      Run test only

    Returns
    -------

    """
    pdb.set_trace()  # THIS VERSION IS NOW FROZEN
    raise IOError("THIS VERSION IS NOW FROZEN")
    version = 'v01'
    # HDF5 file
    outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_{:s}.hdf5'.format(version)
    # Chk clobber
    if os.path.isfile(outfil):
        if clobber:
            warnings.warn("Overwriting previous DB file {:s}".format(outfil))
        else:
            warnings.warn("Not overwiting previous DB file.  Use clobber=True to do so")
            return
    # Begin
    hdf = h5py.File(outfil,'w')

    ''' Myers QSOs '''
    myers.orig_add_to_hdf(hdf)

    # Main DB Table
    idkey = 'IGM_ID'
    maindb, tkeys = sdbbu.start_maindb(idkey)

    # Group dict
    group_dict = {}

    # Organize for main loop
    groups = get_build_groups(version)

    pair_groups = ['SDSS_DR7']

    meta_only = False
    # Loop over the groups
    for gname in groups:
        # Meta
        if gname == 'SDSS_DR7':
            meta = groups[gname].grab_meta(hdf)
        else:
            meta = groups[gname].grab_meta()
        # Survey flag
        flag_g = sdbbu.add_to_group_dict(gname, group_dict)
        # IDs
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups))
        # Spectra
        if not meta_only:
            groups[gname].hdf5_adddata(hdf, gname, meta, idkey)

    # Check for duplicates -- There is 1 pair in SDSS (i.e. 2 duplicates)
    if not sdbbu.chk_for_duplicates(maindb, dup_lim=2):
        raise ValueError("Failed duplicates")

    # Check for junk
    zpri = defs.z_priority()

    # Finish
    sdbbu.write_hdf(hdf, str('igmspec'), maindb, zpri,
                    group_dict, version, Publisher=str(publisher))
    print("Wrote {:s} DB file".format(outfil))
    print("Update DB info in specdb.defs.dbase_info !!")


def ver02(test=False, skip_copy=False, publisher='J.X. Prochaska', clobber=False,
          version='v02', out_path=None):
    """ Build version 2.X

    Reads previous datasets from v1.X

    Parameters
    ----------
    test : bool, optional
      Run test only
    skip_copy : bool, optional
      Skip copying the data from v01

    Returns
    -------
    """
    import os
    from specdb.specdb import IgmSpec
    # Read v01
    v01file = os.getenv('SPECDB')+'/IGMspec_DB_v01.hdf5'
    #v01file_debug = igmspec.__path__[0]+'/tests/files/IGMspec_DB_v01_debug.hdf5'
    print("Loading v01")
    igmsp_v01 = IgmSpec(db_file=v01file)
    v01hdf = igmsp_v01.hdf
    maindb = igmsp_v01.cat.copy()

    # Start new file
    if out_path is None:
        out_path = igmspec.__path__[0]+'/../DB/'
    outfil = out_path + 'IGMspec_DB_{:s}.hdf5'.format(version)
    # Clobber?
    if not chk_clobber(outfil, clobber=clobber):
        return
    # Begin
    hdf = h5py.File(outfil,'w')


    # Copy over the old stuff
    redo_groups = ['HD-LLS_DR1']
    skip_groups = []#'BOSS_DR12', 'SDSS_DR7'] #warnings.warn("NEED TO PUT BACK SDSS AND BOSS!")
    skip_copy = False
    if (not test) and (not skip_copy):
        old_groups = get_build_groups('v01')
        for key in v01hdf.keys():
            if key in ['catalog','quasars']+redo_groups+skip_groups:
                continue
            else:
                #v01hdf.copy(key, hdf)  # ONE STOP SHOPPING
                grp = hdf.create_group(key)
                # Copy spectra
                v01hdf.copy(key+'/spec', hdf[key])
                # Modify v01 meta and add
                if key == 'BOSS_DR12':
                    meta = boss.add_coflag(v01hdf)
                else:
                    meta = Table(v01hdf[key+'/meta'].value)
                meta.rename_column('GRATING', 'DISPERSER')
                hdf[key+'/meta'] = meta
                for akey in v01hdf[key+'/meta'].attrs.keys():
                    hdf[key+'/meta'].attrs[akey] = v01hdf[key+'/meta'].attrs[akey]
                # SSA info
                old_groups[key].add_ssa(hdf, key)
    skip_myers = False
    if skip_myers:
        warnings.warn("NEED TO INCLUDE MYERS!")
    else:
        myers.add_to_hdf(hdf)

    # Setup groups
    old_groups = get_build_groups('v01')
    pair_groups = []
    group_dict = igmsp_v01.qcat.group_dict
    # Set/Check keys (and set idkey internally for other checks)
    idkey = 'IGM_ID'
    _, tkeys = sdbbu.start_maindb(idkey)
    mkeys = list(maindb.keys())
    for key in tkeys:
        assert key in mkeys

    # Loop over the old groups to update (as needed)
    new_IDs = False
    for gname in redo_groups:
        print("Working to replace meta/spec for group: {:s}".format(gname))
        # Meta
        meta = old_groups[gname].grab_meta()
        # Group flag
        flag_g = group_dict[gname]
        # IDs
        if new_IDs:
            pdb.set_trace()  # NOT READY FOR THIS
            #maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
            #                   first=(flag_g==1), close_pairs=(gname in pair_groups))
        else:
            _, _, ids = sdbbu.set_new_ids(maindb, meta, idkey)
        # Spectra
        old_groups[gname].hdf5_adddata(hdf, gname, meta)
        old_groups[gname].add_ssa(hdf, gname)

    meta_only = False
    new_groups = get_build_groups(version)
    # Loop over the new groups
    for gname in new_groups:
        print("Working on group: {:s}".format(gname))
        # Meta
        meta = new_groups[gname].grab_meta()
        # Survey flag
        flag_g = sdbbu.add_to_group_dict(gname, group_dict, skip_for_debug=True)
        # IDs
        debug= False
        #if gname == 'XQ-100':
        #    debug = True
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups),
                               debug=debug)
        # Spectra
        if not meta_only:
            new_groups[gname].hdf5_adddata(hdf, gname, meta)
            new_groups[gname].add_ssa(hdf, gname)

    # Check for duplicates -- There is 1 pair in SDSS (i.e. 2 duplicates)
    if not sdbbu.chk_for_duplicates(maindb, dup_lim=2):
        raise ValueError("Failed duplicates")

    # Check stacking
    if not sdbbu.chk_vstack(hdf):
        print("Meta data will not stack using specdb.utils.clean_vstack")
        print("Proceed to write at your own risk..")
        pdb.set_trace()

    # Finish
    zpri = v01hdf['catalog'].attrs['Z_PRIORITY']
    sdbbuwrite_hdf(hdf, str('igmspec'), maindb, zpri,
                    group_dict, version, Publisher=str(publisher))

    print("Wrote {:s} DB file".format(outfil))
    print("Update DB info in specdb.defs.dbase_info !!")


def ver03(test=False, skip_copy=False, publisher='J.X. Prochaska', clobber=False,
          version='v03.1', out_path=None, redo_dr14=False):
    """ Build version 3.X
    Reads several previous datasets from v1.X
    Remakes the maindb using BOSS DR14 as the main driver

    v3.0 will be BOSS DR14 only to speed up generation of the rest

    Parameters
    ----------
    test : bool, optional
      Run test only
    skip_copy : bool, optional
      Skip copying the data from v01

    Returns
    -------
    """
    import os
    from specdb.specdb import IgmSpec
    # Read v02
    v02file = os.getenv('SPECDB')+'/IGMspec_DB_v02.1.hdf5'
    igmsp_v02 = IgmSpec(db_file=v02file)
    v02hdf = igmsp_v02.hdf

    # Start new file
    if out_path is None:
        out_path = '/scratch/IGMSpec/'
    outfil = out_path + 'IGMspec_DB_{:s}.hdf5'.format(version)
    # Clobber?
    if not chk_clobber(outfil, clobber=clobber):
        return

    # Other bits
    pair_groups = ['SDSS_DR7']

    # Begin
    hdf = h5py.File(outfil,'w')

    # Set/Check keys (and set idkey internally for other checks)
    idkey = 'IGM_ID'
    maindb, tkeys = sdbbu.start_maindb(idkey)
    group_dict = {}

    # BOSS DR14
    new_groups = get_build_groups('v03')
    gname = 'BOSS_DR14'
    # Survey flag
    flag_g = sdbbu.add_to_group_dict(gname, group_dict, skip_for_debug=True)

    if not redo_dr14:
        v030file = os.getenv('SPECDB')+'/IGMspec_DB_v03.0.hdf5'
        igmsp_v030 = IgmSpec(db_file=v030file)

        grp = hdf.create_group(gname)
        # Copy spectra
        igmsp_v030.hdf.copy(gname+'/spec', hdf[gname])
        # Copy meta
        igmsp_v030.hdf.copy(gname+'/meta', hdf[gname])
        # Meta for maindb (a little risky as Meta needs to be aligned to the spectra but they should be)
        meta = igmsp_v030['BOSS_DR14'].meta
        meta.remove_column('IGM_ID')
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups),
                               debug=False)
    else:
        # BOSS DR14
        print("Working on group: {:s}".format(gname))
        # Meta
        meta = new_groups[gname].grab_meta()
        # IDs
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups),
                               debug=False)
        # Spectra
        new_groups[gname].hdf5_adddata(hdf, gname, meta)
        new_groups[gname].add_ssa(hdf, gname)

    # Pop me
    new_groups.pop('BOSS_DR14')

    # Loop on new v3 groups before copying in the others
    for gname in new_groups.keys():
        print("Working on group: {:s}".format(gname))
        # Meta
        meta = new_groups[gname].grab_meta()
        # Survey flag
        flag_g = sdbbu.add_to_group_dict(gname, group_dict, skip_for_debug=True)
        # IDs
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups),
                               debug=False)
        # Spectra
        new_groups[gname].hdf5_adddata(hdf, gname, meta)
        new_groups[gname].add_ssa(hdf, gname)

    # Copy over all the old stuff
    redo_groups = []#'HD-LLS_DR1']
    skip_groups = ['BOSS_DR12']# 'SDSS_DR7'] #warnings.warn("NEED TO PUT BACK SDSS AND BOSS!")
    skip_copy = False
    if (not test) and (not skip_copy):
        old1 = get_build_groups('v01')
        old2 = get_build_groups('v02')
        # Add v02 to v01 list
        for key,item in old2.items():
            old1[key] = item
        # Loop on the combined
        for key in old1.keys():
            if key in ['catalog']+redo_groups+skip_groups:
                continue
            print("Working on: {:s}".format(key))
            grp = hdf.create_group(key)
            # Meta
            meta = Table(v02hdf[key+'/meta'].value)
            meta.remove_column('IGM_ID')
            # Survey flag
            flag_g = sdbbu.add_to_group_dict(key, group_dict, skip_for_debug=True)
            # IDs
            maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                                   first=(flag_g==1), close_pairs=(key in pair_groups),
                                   debug=False)
            # Add meta to HDF5
            #meta.rename_column('GRATING', 'DISPERSER')
            hdf[key+'/meta'] = meta
            for akey in v02hdf[key+'/meta'].attrs.keys():
                hdf[key+'/meta'].attrs[akey] = v02hdf[key+'/meta'].attrs[akey]
            # SSA info
            old1[key].add_ssa(hdf, key)
            # Copy spectra
            v02hdf.copy(key+'/spec', hdf[key])


    skip_myers = False
    if skip_myers:
        warnings.warn("NEED TO INCLUDE MYERS!")
    else:
        # Copy from v02
        _ = hdf.create_group('quasars')
        v02hdf.copy('quasars', hdf['quasars'])
        #myers.add_to_hdf(hdf)

    # Setup groups
    pair_groups = []

    # Check for duplicates -- There is 1 pair in SDSS (i.e. 2 duplicates)
    if not sdbbu.chk_for_duplicates(maindb, dup_lim=2):
        raise ValueError("Failed duplicates")

    # Check stacking
    if not sdbbu.chk_vstack(hdf):
        print("Meta data will not stack using specdb.utils.clean_vstack")
        print("Proceed to write at your own risk..")
        pdb.set_trace()

    # Finish
    zpri = v02hdf['catalog'].attrs['Z_PRIORITY']
    sdbbu.write_hdf(hdf, str('igmspec'), maindb, zpri,
                    group_dict, version, Publisher=str(publisher))

    print("Wrote {:s} DB file".format(outfil))
    print("Update DB info in specdb.defs.dbase_info !!")


def chk_clobber(outfil, clobber=False):
    """ Simple clobber check
    outfil : str
    clobber : bool, optional
    """
    # Chk clobber
    if os.path.isfile(outfil):
        if clobber:
            warnings.warn("Overwriting previous DB file {:s}".format(outfil))
            return True
        else:
            warnings.warn("Not overwiting previous DB file.  Set clobber=True to do so")
            return False
    else:
        return True


def get_build_groups(version):
    """
    Parameters
    ----------
    version : str

    Returns
    -------
    build_groups : dict

    """

    groups = OrderedDict()
    if version == 'v01':
        groups['BOSS_DR12'] = boss
        groups['SDSS_DR7'] = sdss
        groups['KODIAQ_DR1'] = kodiaq
        groups['HD-LLS_DR1'] = hdlls
        groups['GGG'] = ggg
    elif version[0:3] == 'v02':
        groups['HST_z2'] = hst_z2         # O'Meara et al. 2011
        groups['XQ-100'] = xq100          # Lopez et al. 2016
        groups['HDLA100'] = hdla100       # Neeleman et al. 2013
        groups['2QZ'] = twodf             # Croom et al.
        groups['ESI_DLA'] = esidla        # Rafelski et al. 2012, 2014
        groups['COS-Halos'] = cos_halos   # Tumlinson et al. 2013
        groups['COS-Dwarfs'] = cos_dwarfs # Bordoloi et al. 2014
        groups['HSTQSO'] = hst_qso        # Ribaudo et al. 2011; Neeleman et al. 2016
        groups['MUSoDLA'] = musodla       # Jorgensen et al. 2013
        groups['UVES_Dall'] = uves_dall   # Dall'Aglio et al. 2008
        groups['UVpSM4'] = hst_c          # Cooksey et al. 2010, 2011
    elif version == 'v03':
        groups['BOSS_DR14'] = boss_dr14   # Paris et al. 2018
        groups['SQUAD_DR1'] = uves_squad  # Murphy et al. 2018
        groups['ESI_z6'] = esi_z6         # Eiler et al. 2018
        groups['KODIAQ_DR2'] = kodiaq_two
    else:
        raise IOError("Not ready for this version")
    # Return
    return groups
