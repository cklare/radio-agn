# plot ztf light curves
# chloe klare 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import os
import juliandate as jd
from datetime import datetime

# first, get candidate coordinates so I know which light curve goes with which
agn_coords = pd.read_pickle('../../final_paper/agn_coords')
# get the ids of all the candidates
main_folder = '../../final_paper/no_relations/lightcurves/'
good_sources = os.scandir(main_folder)
candidate_ids = []
for file in good_sources:
    objid = file.name.replace(".png","")
    candidate_ids += [objid]

print(f'number of sources: {len(candidate_ids)}')

# function for plotting one of the ztf bands
def plot_one_color(df, axs, color,label):
    # first - get dates of measurements and convert to elapsed time
    # any optical flares will not line up with the radio flares anyway, so the actual dates are less important
    julian_date = df['jd']
    first_epoch = julian_date.iloc[0]
    days_elapsed = julian_date-(np.ones((len(julian_date),),)*first_epoch)
    mag = df['mag']
  
    # plot just the magnitudes since this is the most meaningful
    ax.scatter(days_elapsed,mag,color=color,label=label)
  
    return axs

# get coordinates of all candidates
mask = agn_coords['ObjID'].astype(str).isin(candidate_ids)
agn_coords = agn_coords[mask]

# folder with all the ztf data
ztf_data = '../../final_paper/ZTF/'
ztf_data = os.scandir(ztf_data)

for file in ztf_data:
    file_name = file.name
    file_name = file_name.replace('.txt','')
    # read in the file
    source_data = pd.read_csv(file,sep=' ',names=['index', 'field', 'ccdid', 'qid', 'filter', 'pid', 'infobitssci', 'sciinpseeing',\
    'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms', 'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2',\
    'diffmaglim', 'zpdiff', 'programid', 'jd', 'rfid', 'forcediffimflux', 'forcediffimfluxunc', 'forcediffimsnr', 'forcediffimchisq', 'forcediffimfluxap',\
    'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr', 'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi', 'nearestrefsharp',\
    'refjdstart', 'refjdend', 'procstatus'],skiprows=58)


    source_data = source_data.set_index('index')
    # the last line is a footer 
    source_data = source_data[:-1]
    # only keep measurements with no warnings     
    mask = (source_data['procstatus']==0) 
    source_data = source_data[mask]
   
    # first, remove any bad epochs as suggested in the ZTF forced photometry paper
    mask = source_data['infobitssci']>0
    source_data = source_data[~mask]
    
    mask = source_data['scisigpix']<25
    source_data = source_data[mask]
    
    mask = source_data['sciinpseeing']<4
    source_data = source_data[mask]
    mask = source_data['forcediffimchisq']<2
    source_data = source_data[mask]

    # if we have no epochs left, move to the next file
    if len(source_data)==0:
        continue

    # now, calculations
    # the ab magnitude of the reference source
    m_ref = source_data['nearestrefmag'].values
    # the ab magnitude a source detected with f=1 DN by ZTF would have
    zpdiff = source_data['zpdiff'].values
    
    # convert the magnitude to a flux (in counts)
    flux_ref = 10**(0.4*(zpdiff-m_ref)) # in counts
    # the forced difference flux of my source
    flux_diffs = source_data['forcediffimflux'].values

    # determine the error in the reference flux
    flux_ref_err = np.sqrt((10**(0.4*zpdiff))**2*(10**(-0.4*m_ref))**2*(-0.4*np.log(10)*source_data['nearestrefmagunc'])**2)

    # add the difference flux to the flux it was substracted from
    flux = flux_ref+source_data['forcediffimflux'].values
    flux_err = np.sqrt(flux_ref_err**2+source_data['forcediffimfluxunc']**2)

    # add new columns to the df
    source_data = source_data.assign(total_flux=flux)
    source_data = source_data.assign(total_flux_err=flux_err)
    # remove the sources with negative total flux counts bc idk what is going on there
    mask = source_data['total_flux']>0
    print(f'total negative sources: {sum(~mask)}/{len(mask)}')
    source_data = source_data[mask]

    flux = source_data['total_flux']
    flux_err = source_data['total_flux_err']
    zpdiff = source_data['zpdiff']
    # now, calculate the AB magnitude of our sources
    mag = zpdiff-2.5*np.log10(flux)
    # convert this to an analog flux value (in Jy)
    flux_analog = 3631*10**(-mag/2.5)
    # determine the error in AB magnitude
    mag_err = np.sqrt((-2.5*flux_err/(flux*np.log(10)))**2)
    # determine the error in the analog flux (Jy)
    flux_analog_err = np.sqrt(flux_analog**2*(-1/2.5*np.log(10)*mag_err)**2)
    # add these to my df
    source_data = source_data.assign(mag=mag)
    source_data = source_data.assign(flux_analog=flux_analog)
    source_data = source_data.assign(flux_analog_err=flux_analog_err)
    
    # only consider sources with snr>3 as detections (also following ZTF paper)
    mask = source_data['flux_analog']/source_data['flux_analog_err']>3
    source_data = source_data[mask]

    # if nothing left to plot, continue
    if len(source_data)==0:
        continue

    # plot the light curves
    f,ax = plt.subplots(1)
    f.set_size_inches(12,8)
    mask = source_data['filter']=='ZTF_g'
    if sum(mask)>0:
        ax = plot_one_color(source_data[mask],ax,'g','g')

    mask = source_data['filter']=='ZTF_r'
    if sum(mask)>0:
        axs= plot_one_color(source_data[mask],ax,'r','r')

    mask = source_data['filter']=='ZTF_i'
    if sum(mask)>0:
        ax = plot_one_color(source_data[mask],ax,'orange', 'i')
    
    ax.set_xlabel('Days Elapsed')
    ax.set_ylabel('AB Magnitude')
    ax.yaxis.set_inverted(True)
    plt.legend(loc='2')
  
 
    f.suptitle('Optical Lightcurve [ZTF]')
    plt.tight_layout()
    plt.savefig(f'../../final_paper/ZTF_lightcurves/{file_name}.png')
    plt.close()
              



                