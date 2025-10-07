# this script queries online AGN catalogs as a starting point for our crossmatch
# chloe klare

import numpy as np
import pandas as pd
import pickle
import os
import astropy.units as u # working with unit conversions
from astropy.coordinates import SkyCoord, search_around_sky
from astroquery.sdss import SDSS
from astroquery.vizier import Vizier
from astroquery.heasarc import Heasarc

# function for removing suspected blazars
def remove_blazars(agn_cat, blazar_coords):
    n_agn = len(agn_cat)
    # create a coordinate object
    agn_coords = SkyCoord(agn_cat['ra'].values, agn_cat['dec'].values, unit=u.deg)
    # perform the catalog crossmatch
    bz_idx, bz_cat_idx, sep, dist = search_around_sky(agn_coords,blazar_coords, 1.5*u.arcsec)
    bz_ids = agn_cat['ObjID'].iloc[bz_idx]
    # drop all possible blazar contaminants
    agn_cat_mask = ~agn_cat['ObjID'].isin(bz_ids)
    agn_cat = agn_cat[agn_cat_mask]
    # print results
    print(f'dropped {n_agn-len(agn_cat)} blazars')
    del(sep, dist, bz_ids, bz_idx, bz_cat_idx)
    # return the decontaminated catalog of agn
    return agn_cat

# first, query SDSS catalog
max_sep = 2.5*u.arcsec
#************************************** querying the initial catalogs of AGN *****************************

# first, query the sdss database
print('querying sdss...')
# define the SDSS sql query
SDSS_query = "SELECT top 500000 p.ObjID, p.ra, p.dec, p.type, p.PSFMag_g, s.z, s.zwarning, s.class FROM \
PhotoObj AS p JOIN SpecObj AS s ON s.bestobjid = p.objid WHERE s.class = 'qso' AND \
s.zwarning = 0 AND p.dec BETWEEN -71.3 AND 30.0"
# execute the query
sdss_data = SDSS.query_sql(SDSS_query)
del(SDSS_query)
# output the size of the query
print('size of sdss query: %d' %len(sdss_data))

# convert to pandas dataframe, which is easier to work with
sdss_data = sdss_data.to_pandas()
# now, repeat for wise catalog
# wise database query
print('querying wise...')
agn_search = Vizier(columns = ['WISEA','RAJ2000', 'DEJ2000', 'W1mag'],row_limit=-1)
wise_data = agn_search.query_constraints(catalog='J/ApJS/234/23/r90cat',DEJ2000= '<30.0' )
wise_data = wise_data[0]
print('size of wise query: %d' %len(wise_data))
# convert to pandas dataframe, which is easier to work with
wise_data = wise_data.to_pandas()
# rename this to the sdss conventions, which is easier
column_mapping = {'WISEA':'ObjID','RAJ2000':'ra', 'DEJ2000':'dec'}
wise_data = wise_data.rename(columns=column_mapping)

# first, remove any sources near a blazar
print('beginning blazar-ectomy...')
# blazar-ectomy - removing blazars from my sample
# load the heasarc query class
heasarc = Heasarc() # class object
# retrieve the blazar catalog
# define query
request_payload = heasarc._args_to_payload(entry = '', mission = 'romabzcat',fields = 'ra, dec', resultmax = 4000)
# execute the query
blazar_cat = heasarc.query(request_payload)
del(request_payload)
blazar_cat = blazar_cat.to_pandas()

print(f'number of blazars in heasarc catalog: {len(blazar_cat)}')

# define our query - we only need the coordinates
blazar_search = Vizier(columns = ['RAJ2000', 'DEJ2000'],row_limit=-1)
# get the catalog we want from the paper - this gets us a list of catalogs
blazar_candidates = blazar_search.query_constraints(catalog= 'J/ApJS/215/14/table4')
# get the actual data
blazar_candidates = blazar_candidates[0]
blazar_candidates = blazar_candidates.to_pandas()
print(f'number of blazars in the wise blazar catalog: {len(blazar_candidates)}')
blazar_candidates = blazar_candidates.rename(columns=column_mapping)
blazar_cat = blazar_cat.rename(columns={'RA':'ra','DEC':'dec'})
# merge both blazar catalogs into one df
blazar_cat = pd.concat([blazar_cat, blazar_candidates], axis=0)
del(blazar_candidates)
print(f'number of blazars in catalog to be xmatched: {len(blazar_cat)}')
blazar_coords = SkyCoord(blazar_cat['ra'].values, blazar_cat['dec'].values,unit=u.deg)
# remove the blazar suspects from both agn samples
sdss_data = remove_blazars(sdss_data, blazar_coords)
print(f'sdss sources after removing blazars: {len(sdss_data)}')
wise_data = remove_blazars(wise_data, blazar_coords)
print(f'wise sources after removing blazars: {len(wise_data)}')


# convert to numpy arrays, which makes skycoord run much faster
wise_coords = SkyCoord(wise_data['ra'].values, wise_data['dec'].values, unit =u.deg)
sdss_coords = SkyCoord(sdss_data['ra'].values, sdss_data['dec'].values, unit=u.deg)

# split the two dataframes into 3 dataframes - sdss only sources, sdss/wise sources, and wise sources
sdss_sources_in_wise, wise_sources_in_sdss, sep, dist =search_around_sky(sdss_coords, wise_coords, 1.5*u.arcsec)
del(sdss_coords, wise_coords)

print(f'number of AGN matches between SDSS and WISE catalogs: {len(sdss_sources_in_wise)} {len(wise_sources_in_sdss)}')
print(f'unique sdss sources matched to a wise source: {len(np.unique(sdss_sources_in_wise))}')
print(f'unique wise sources matched to an sdss source: {len(np.unique(wise_sources_in_sdss))}')


# let's try dropping any wise sources mapped to an sdss source, since the sdss source has more accurate data
overlapping_sources = wise_data.iloc[wise_sources_in_sdss]
mask = wise_data['ObjID'].isin(overlapping_sources['ObjID'])
wise_only = wise_data[~mask]
wise_coords = wise_only[['ObjID','ra','dec']]
sdss_coords = sdss_data[['ObjID','ra','dec']]
agn_coords = pd.concat([sdss_coords,wise_coords],axis=0, ignore_index=True)


print(f'total agn in catalog:{len(agn_coords)}')
print(agn_coords)

agn_data_folder = '../../final_paper/'
if not os.path.isdir(agn_data_folder):
    os.mkdir(agn_data_folder)

# save the datafiles
agn_coords.to_pickle(f'{agn_data_folder}agn_coords')
wise_data.to_pickle(f'{agn_data_folder}wise_data')
sdss_data.to_pickle(f'{agn_data_folder}sdss_data')
