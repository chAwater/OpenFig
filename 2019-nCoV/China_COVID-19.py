#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os

from tqdm.auto import tqdm

import numpy  as np
import pandas as pd

import imageio
import matplotlib.pyplot as plt

from matplotlib.colors import Normalize

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader


# ## Data for the Map
# 
# - Download data from [GADM](https://gadm.org)

# In[2]:


name2map = {
    'Neimenggu':'Nei Mongol',
    'Ningxia':'Ningxia Hui',
    'Xinjiang':'Xinjiang Uygur',
    'Macau':'Macao',
}

geo_data = {}

# 中国大陆
for record in shpreader.Reader('./data_GADM/gadm36_CHN_1.shp').records():
    name = record.attributes['NAME_1']
    geo  = record.geometry
    geo_data[name] = geo

# 香港、澳门、台湾
for sp in ['HKG','MAC','TWN']:
    record = list(shpreader.Reader('./data_GADM/gadm36_{:s}_0.shp'.format(sp)).records())[0]
    name = record.attributes['NAME_0']
    geo = record.geometry
    geo_data[name] = geo


# ## Load data
# 
# - Download data from https://github.com/BlankerL/DXY-COVID-19-Data
# - Keep the max value for a day
# - Define `existed` = `confirmed` - `cured` - `dead`

# In[3]:


date_str = pd.date_range(start='2020-01-22',end='2020-03-23').strftime('%Y-%m-%d')

china_df = (
    pd.read_csv('./DXYArea.csv')
    .query('countryEnglishName=="China"')
    .query('province_confirmedCount!=0')
    .loc[:,['provinceEnglishName','province_confirmedCount','province_curedCount','province_deadCount','updateTime']]
    .rename(columns={'provinceEnglishName':'Name','province_confirmedCount':'confirmed','province_curedCount':'cured','province_deadCount':'dead'})
    .copy()
)

china_df['date_str'] = china_df['updateTime'].apply(lambda x:pd.to_datetime(x).date().strftime('%Y-%m-%d'))

china_df = (
    china_df
    # Keep the max value for a day
    .groupby(['Name','date_str'],sort=False).max()
    .reset_index(level=1)
    # Drop summary
    .drop('China')
    # Rename for some provinces
    .rename(index=name2map)
    # Sort by date
    .set_index('date_str',append=True)
)

# Define existed = confirmed - cured - dead
china_df['existed'] = (china_df['confirmed'] - china_df['cured'] - china_df['dead'])


# ### Clean data
# - linear interpolate
# - fillna(0)
# - +0.1, log10

# In[4]:


# Clean data
clean_df = pd.DataFrame(
    index=pd.MultiIndex.from_product([geo_data.keys(),date_str]),
    columns=china_df.columns,
)

# For each province
for g,df in china_df.reset_index().groupby('Name',sort=False):
    ndf = (
        df.set_index('date_str')
        # Sort by date
        .reindex(index=date_str)
        # Fill the name
        .assign(Name=g)
        .set_index('Name',append=True)
        # Interpolate missing values
        .interpolate()
        # Fillna by 0
        .fillna(0)
    )
    clean_df.loc[ndf.index.swaplevel(0,1)] = ndf.values

clean_df['color'] = plt.cm.Reds(
    Normalize(0, 5, clip=True)(
        np.log10( clean_df['existed'].astype(int)+0.1 ).values
    )
).tolist()


# ---

# In[5]:


files = []
for n,d in tqdm(enumerate(date_str), total=date_str.shape[0]):
    
    file_name = 'frames/frame_{:02d}.jpg'.format(n)
    files.append(file_name)
    
    if os.path.exists(file_name):
        continue
    
    ndf = clean_df.xs(d,level=1)

    # Plot the main land
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection=ccrs.LambertConformal(central_latitude=90,central_longitude=105))
    ax.set_extent([80, 130, 13, 55])
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='white', alpha=0.5)
    ax.add_feature(cfeature.OCEAN.with_scale('110m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), zorder=1000)
    ax.add_feature(cfeature.BORDERS.with_scale('110m'), zorder=1000)

    # Plot provinces
    for k in geo_data.keys():
        if k in ndf.index:
            gcolor = ndf.loc[k,'color']
        else:
            gcolor = 'white'
    
        ax.add_geometries(
            geo_data[k],
            crs=ccrs.PlateCarree(),
            facecolor=gcolor,
            lw=0.1,
            edgecolor='k',
            zorder=0
        )
    
    ax.set_title(d, fontsize=24)
    
    # Add color bar
    cax = fig.add_axes([0.825, 0.2, 0.02, 0.2])
    fig.colorbar(
        plt.cm.ScalarMappable(norm=Normalize(0, 5, clip=True), cmap='Reds'), 
        cax=cax,
#         extend='max',
    )
    cax.set_yticklabels(['$10^{:d}$'.format(x) for x in np.arange(0,6)], fontsize=12, ha='left',va='center')
    
    fig.savefig(file_name, dpi=150, facecolor=None)
    
    plt.close(fig)
    


# ---

# ## Merge jpg to gif

# In[7]:


get_ipython().system(' ffmpeg -f image2 -framerate 6 -y -i ./frames/frame_%002d.jpg China.gif')


# In[6]:


# with imageio.get_writer('./China.gif', mode='I') as writer:
#     for file_name in files:
#         image = imageio.imread(file_name)
#         writer.append_data(image)
        
# from pygifsicle import optimize
# optimize('./China.gif')


# ---
# ---
# ---

# In[8]:


# fig = plt.figure(figsize=(8,8))
# ax = fig.add_subplot(projection=ccrs.LambertConformal(central_latitude=90,central_longitude=105))

# ax.set_extent([80, 130, 13, 55])
# ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='white', alpha=0.5)
# ax.add_feature(cfeature.OCEAN.with_scale('110m'))
# ax.add_feature(cfeature.COASTLINE.with_scale('110m'), zorder=1000)
# ax.add_feature(cfeature.BORDERS.with_scale('110m'), zorder=1000)

# for k in geo_data.keys():
#     if k in ndf.index:
#         c = ndf.loc[k,'color']
#     else:
#         c = 'white'
#     ax.add_geometries(
#         geo_data[k],
#         crs=ccrs.PlateCarree(),
#         facecolor=c,
#         lw=0.1,
#         edgecolor='k',
# #         hatch='//',
#         zorder=0
        
#     )
    
#     ax.set_title(d, fontsize=24)
    
#     cax = fig.add_axes([0.825, 0.2, 0.02, 0.2])
#     fig.colorbar(
#         plt.cm.ScalarMappable(norm=Normalize(0, 5), cmap='Reds'), 
#         cax=cax
#     )
    
#     cax.set_yticklabels(['$10^0$','$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'], fontsize=12, ha='left',va='center')
    
#     plt.show()


# In[9]:


# # https://gmt-china.org/
# with open('CN-border-La.dat') as src:
#     context = src.read()
#     blocks = [cnt for cnt in context.split('>') if len(cnt) > 0]
#     borders = [np.fromstring(block, dtype=float, sep=' ') for block in blocks]

# for line in borders:
#     ax.plot(
#         line[0::2], line[1::2],
#         '-', lw=1, color='k', 
#         transform=ccrs.Geodetic(),
#     )


# In[ ]:





# In[ ]:




