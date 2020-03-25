#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os

from tqdm.auto import tqdm
from multiprocessing import Pool

import numpy  as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.cm        import ScalarMappable
from matplotlib.colors    import Normalize, LogNorm
from matplotlib.offsetbox import AnchoredText


import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader


# In[2]:


# https://gist.github.com/JeffPaine/3083347
# https://gist.github.com/tlancon/9794920a0c3a9990279de704f936050c
US_STATES = {
	'Alabama': 'AL',
	'Alaska': 'AK',
	'Arizona': 'AZ',
	'Arkansas': 'AR',
	'California': 'CA',
	'Colorado': 'CO',
	'Connecticut': 'CT',
	'Delaware': 'DE',
	'District of Columbia': 'DC',
	'Florida': 'FL',
	'Georgia': 'GA',
	'Hawaii': 'HI',
	'Idaho': 'ID',
	'Illinois': 'IL',
	'Indiana': 'IN',
	'Iowa': 'IA',
	'Kansas': 'KS',
	'Kentucky': 'KY',
	'Louisiana': 'LA',
	'Maine': 'ME',
	'Maryland': 'MD',
	'Massachusetts': 'MA',
	'Michigan': 'MI',
	'Minnesota': 'MN',
	'Mississippi': 'MS',
	'Missouri': 'MO',
	'Montana': 'MT',
	'Nebraska': 'NE',
	'Nevada': 'NV',
	'New Hampshire': 'NH',
	'New Jersey': 'NJ',
	'New Mexico': 'NM',
	'New York': 'NY',
	'North Carolina': 'NC',
	'North Dakota': 'ND',
	'Ohio': 'OH',
	'Oklahoma': 'OK',
	'Oregon': 'OR',
	'Pennsylvania': 'PA',
	'Rhode Island': 'RI',
	'South Carolina': 'SC',
	'South Dakota': 'SD',
	'Tennessee': 'TN',
	'Texas': 'TX',
	'Utah': 'UT',
	'Vermont': 'VT',
	'Virginia': 'VA',
	'Washington': 'WA',
	'West Virginia': 'WV',
	'Wisconsin': 'WI',
	'Wyoming': 'WY'
}


# In[3]:


name2map = {
    'Inner Mongolia':'Nei Mongol',
    'Ningxia':'Ningxia Hui',
    'Xinjiang':'Xinjiang Uygur',
    'Macau':'Macao',
    'Tibet':'Xizang'
}


# In[4]:


china_geo_data = {}
usa_geo_data = {}

# 中国大陆
for record in shpreader.Reader('./data_GADM/gadm36_CHN_1.shp').records():
    name = record.attributes['NAME_1']
    geo  = record.geometry
    china_geo_data[name] = geo

# 香港、澳门、台湾
for sp in ['HKG','MAC','TWN']:
    record = list(shpreader.Reader('./data_GADM/gadm36_{:s}_0.shp'.format(sp)).records())[0]
    name = record.attributes['NAME_0']
    geo  = record.geometry
    china_geo_data[name] = geo

# USA
for record in shpreader.Reader(shpreader.natural_earth(resolution='110m',category='cultural', name='admin_1_states_provinces')).records():
    name = record.attributes['name']
    geo  = record.geometry
    usa_geo_data[name] = geo


# In[5]:


color_mapper = ScalarMappable(norm=Normalize(0,5,clip=True), cmap='Reds')
# color_mapper = ScalarMappable(norm=LogNorm(1,5,clip=True), cmap='Reds')
# plt.scatter(np.arange(1,5e4,50),np.arange(1,5e4,50),c=color_mapper.to_rgba(np.log10(np.arange(1,5e4,50))))


# In[6]:


usa_data = []
china_data = []

for file_type in ['Confirmed','Deaths','Recovered']:
    # Load csv
    total_data_df = pd.read_csv('time_series_19-covid-{:s}.csv'.format(file_type)).set_index('Country/Region')

    # Save TimeSeries as strings
    date_idx = total_data_df.columns.drop(['Province/State','Lat','Long'])
    date_str = pd.to_datetime(date_idx.tolist()).strftime('%Y-%m-%d')
    data_dict = dict(zip(date_idx, date_str))

    usa_data_df   = (
        total_data_df.loc['US']
        .set_index('Province/State')
        .loc[US_STATES.keys()]
        .fillna(0)
        .rename(columns=data_dict)
    )

    china_data_df = (
        total_data_df.loc['China']
        .append(
            total_data_df.loc['Taiwan*'].fillna('Taiwan')
        )
        .set_index('Province/State')
        # Rename for some provinces
        .rename(index=name2map)
        .loc[china_geo_data.keys()]
        .fillna(0)
        .rename(columns=data_dict)
    )

    # Convert to int
    usa_data_df.loc[  :,date_str] = usa_data_df.loc[  :,date_str].astype(int)
    china_data_df.loc[:,date_str] = china_data_df.loc[:,date_str].astype(int)

    usa_data.append(usa_data_df)
    china_data.append(china_data_df)
    
# Define `existed` = `confirmed` - `cured` - `dead`
usa_data_df   = usa_data[0] - usa_data[1] - usa_data[2]
china_data_df = china_data[0] - china_data[1] - china_data[2]

usa_data_df.loc[  :,['Lat','Long']] *= -1
china_data_df.loc[:,['Lat','Long']] *= -1


# ---

# In[7]:


def plot_main_land(ax):
    '''
    Plot the main land, ocean, coastline and borders.
    '''        
    ax.add_feature(cfeature.LAND.with_scale('110m'), facecolor='white', alpha=0.5)
    ax.add_feature(cfeature.OCEAN.with_scale('110m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), zorder=100)
    ax.add_feature(cfeature.BORDERS.with_scale('110m'), zorder=100)

    return 

def plot_states(ax, df, geo_data):
    '''
    Plot provinces/states.
    '''
    for k in geo_data.keys():
        if df[k] == 0:
            gcolor = 'white'
        else:
            gcolor = color_mapper.to_rgba( np.log10(df[k]) )
    
        ax.add_geometries(
            geo_data[k],
            crs=ccrs.PlateCarree(),
            facecolor=gcolor,
            lw=0.1,
            edgecolor='k',
            zorder=0
        )

    cax = ax.inset_axes([0.9, 0.1, 0.02, 0.35])
    plt.colorbar( color_mapper, cax=cax, extend='max', ticks=np.arange(0,6) )
    cax.set_yticklabels(['$10^{:d}$'.format(x) for x in np.arange(0,6)], fontsize=10, ha='left',va='center')

    return


# In[8]:


def two_countries_plot(df1, df2):
    fig = plt.figure(figsize=(14,5))
    ax1 = fig.add_subplot(121, projection=ccrs.LambertConformal(central_latitude=90,central_longitude=105))
    ax2 = fig.add_subplot(122, projection=ccrs.LambertConformal())

    ax1.set_title('China', fontsize=24)
    ax2.set_title('US', fontsize=24)

    ax1.set_extent([80, 130, 15, 53])
    ax2.set_extent([-120, -70, 22, 53])
    
    plot_main_land(ax1)
    plot_main_land(ax2)

    plot_states(ax1, df1, china_geo_data)
    plot_states(ax2, df2, usa_geo_data)
    
    text = AnchoredText(
        'Visualization by ch', 
        loc='lower left', 
        prop={'size': 8, 'alpha':0.25}, frameon=False, 
    )
    ax1.add_artist(text)

    text = AnchoredText(
        'Data from CSSE @ Johns Hopkins University', 
        loc='lower right', 
        bbox_to_anchor=(1.015, -0.15),
        bbox_transform=ax2.transAxes,
        prop={'size': 10}, frameon=True, 
    )
    ax2.add_artist(text)

    text = AnchoredText(
        '@chAwater', 
        loc='lower left', 
        prop={'size': 8, 'alpha':0.25}, frameon=False, 
    )
    ax2.add_artist(text)

    return fig


# In[9]:


def worker(idx):
    day_str = date_str[idx]

    file_name = 'frames/frame_{:02d}.jpg'.format(idx)
    
    ndf1 = china_data_df[day_str]
    ndf2 = usa_data_df[day_str]
    
    fig = two_countries_plot(ndf1, ndf2)
    fig.suptitle('Existing COVID-19\n'+day_str,y=1.1,fontsize=28,ha='center',va='top')
    
    
    fig.savefig(file_name, dpi=150, bbox_inches='tight', facecolor=None)
    
    plt.close(fig)
    
    return 1


# In[ ]:





# In[10]:


pool = Pool(4)
_ = list(
    tqdm( 
        pool.imap(worker, np.arange(date_str.shape[0])), total=date_str.shape[0] 
    )
)
pool.close()
pool.join()


# In[11]:


get_ipython().system(' ffmpeg -f image2 -framerate 8 -y -i ./frames/frame_%002d.jpg China_vs_US.gif')


# In[12]:


get_ipython().system(' ffmpeg -f image2 -framerate 8 -y -i ./frames/frame_%002d.jpg -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" China_vs_US.mp4')


# In[ ]:




