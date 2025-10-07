# AGN-VAST crossmatch
# this script reads in my AGN catalog coordinates, prepares VAST samples, and does a coordinate crossmatch
# Chloe Klare

# note! this code is public, but the VAST data is not, so this script is more for transparency/documentation of my methodology

import numpy as np
import pandas as pd
import pickle
import os
import astropy.units as u # working with unit conversions
from astropy.coordinates import SkyCoord, search_around_sky

# vast packages 
from vasttools.moc import VASTMOCS
from vasttools import source
from vasttools.query import Query
from vasttools.pipeline import Pipeline


# function to actual perform the crossmatch
# sources is the vast-source df and agn_cat is the catalog of AGN to crossmatch
def agn_vast_crossmatch(sources,agn_cat):
   
   # how many sources to start?
    print(f'total number of sources: {len(sources)}')
    
    # first, remove sources with 'relations' (if any sources are mapped to the same source in a different epoch as any other vast source, remove)
    mask = sources['n_relations']==0
    sources = sources[mask]
    print(f'number of sources with no relations: {len(sources)}')

    # keep only compact sources (total flux/peak flux<1.4)
    mask = sources['avg_compactness']<1.4
    print(f'sources with average compactness < 1.4: {sum(mask)}')
    sources = sources[mask]

    # remove sources with any siblings (epochs with another source's beam overlapping its beam)
    mask = (sources['n_siblings']==0)
    print(f'sources with no siblings: {sum(mask)}')
    sources = sources[mask]
   
    print(f'total VAST radio sources to crossmatch: {len(sources)}')
    
    # assign them their index as an additional column
    sources = sources.assign(vast_id=sources.index.values)
    
    # get the coordinates to crossmatch
    ras = agn_cat['ra'].values
    decs = agn_cat['dec'].values
    agn_coords = SkyCoord(ras,decs,unit=u.deg)
    vast_coords = SkyCoord(sources['wavg_ra'].values,sources['wavg_dec'].values,unit=u.deg)
    
    # perform the crossmatch
    print('crossmatching catalogs!')
    idx_matches_agn, idx_matches_vast, sep, dist = search_around_sky(agn_coords,vast_coords,2.5*u.arcsec)
    print(f'unique agn matches: {len(np.unique(idx_matches_agn))}')
    print(f'length of idx_matches_agn: {len(idx_matches_agn)}')
    print(f'unique vast matches: {len(np.unique(idx_matches_vast))}')
    print(f'length of idx_matches_vast: {len(idx_matches_vast)}')
    
    # select all the matches
    agn_matches = agn_cat.iloc[idx_matches_agn]
    vast_matches = sources.iloc[idx_matches_vast]
    # get the identifiers from each
    agn_ids = agn_matches['ObjID'].astype(str).values
    vast_ids = vast_matches['vast_id'].values

    print(f'unique agn matches: {len(np.unique(agn_ids))}')
    print(f'length of idx_matches_agn: {len(agn_ids)}')
    print(f'unique vast matches: {len(np.unique(vast_ids))}')
    print(f'length of idx_matches_vast: {len(vast_ids)}')
    print(f'number of indices in agn cat: {len(agn_cat)} {len(np.unique(agn_cat.index))}')

    # assign the identifiers to the crossmatched catalog to tie them together
    vast_matches = vast_matches.assign(ObjID=agn_ids)
    agn_matches = agn_matches.assign(vast_id=vast_ids)
    sep = (sep*u.arcsec).value
    vast_matches = vast_matches.assign(sep=sep)
    
    # dropping duplicate sources - any sources not uniquely mapped get dropped
    vast_matches = vast_matches.drop_duplicates(subset='ObjID',keep=False)
    agn_matches = agn_matches.drop_duplicates(subset='ObjID',keep=False)
    print(f'number of matches after dropping AGN which were not uniquely mapped: {len(vast_matches)} {len(agn_matches)}')
    vast_matches = vast_matches.drop_duplicates(subset='vast_id',keep=False)
    agn_matches = agn_matches.drop_duplicates(subset='vast_id',keep=False)
    print(f'number of matches after dropping vast sources which were not uniquely mapped: {len(vast_matches)} {len(agn_matches)}')

    # we just want to return the vast matches df (since it has the agn identifiers added)
    print(f'total number of vast-agn matches: {len(vast_matches)}')
    return vast_matches

# define the vast pipeline object 
pipe=Pipeline()

# pull the vast survey (extragalactic fields) data (master source catalog)
print('loading full survey sources...')
full_run = pipe.load_run('workshop2025_extragalactic_skyregfixed')
full_sources = full_run.sources

# read in the agn catalog to crossmatch
agn_cat = pd.read_pickle('final_paper/agn_coords')
print(f'total number of agn: {len(agn_cat)}')

# perform the crossmatch
agn_full_matches = agn_vast_crossmatch(full_sources,agn_cat)

# repeat for the pilot data
pilot_run = pipe.load_run('pilot_p1_redux')
pilot_sources = pilot_run.sources

# perform the pilot matches
agn_pilot_matches = agn_vast_crossmatch(pilot_sources,agn_cat)

# save the matches
agn_pilot_matches.to_pickle('final_paper/agn_pilot_matches')
agn_full_matches.to_pickle('final_paper/agn_full_matches')

# retrieve the radio-AGN measurements for individual epochs
pilot_measurements = pilot_run.measurements
# keep only measurements for sources matched to an agn
mask = (pilot_measurements['source'].isin(agn_pilot_matches['vast_id']))
pilot_measurements = pilot_measurements[mask]
pilot_measurements = pilot_measurements.to_pandas_df()
pilot_measurements.to_pickle('final_paper/pilot_measurements')
del(pilot_measurements,pilot_run)
# repeat with full data

full_measurements = full_run.measurements
mask = (full_measurements['source'].isin(agn_full_matches['vast_id']))
full_measurements = full_measurements[mask]
full_measurements = full_measurements.to_pandas_df()
full_measurements.to_pickle('final_paper/full_measurements')