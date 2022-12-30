import calendar

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from bs4 import BeautifulSoup

data_from_web_txt = open('./Kaggle_Activity_2022.txt', 'r').read()

bsObj = BeautifulSoup(data_from_web_txt, 'html')

fig_df = pd.DataFrame()

n = 0

qdict = {
    'color-empty' : 0,
    'q1' : 1,
    'q2' : 2,
    'q3' : 3,
}

for col in bsObj.find_all("g", {"class":"react-calendar-heatmap-week"}):
    x = int(col['transform'].replace('translate(', '').split(',')[0])
    
    for row in col.find_all('rect'):
        name = row['class'][0]
        time = pd.to_datetime(row['data-tip'][:24])
        y = int(row['y'])
        fig_df = pd.concat(
            [
                fig_df,
                pd.DataFrame.from_dict(
                    {n:{
                        'name' : name,
                        'time' : time,
                        'value' : qdict[name],
                        'x' : x,
                        'y' : y,
                        'month' : calendar.month_abbr[time.month],
                    }}
                ).T
            ]
        )
        n+=1


fig, ax = plt.subplots(figsize=(15,2))
cbar_ax = fig.add_axes([0.15, 0.0, 0.04, 0.08]) # [left, bottom, width, height]

month_idx = fig_df.pivot(columns='x', index='y', values='month').loc[0].drop_duplicates()
heatmap_df = fig_df.pivot(columns='x', index='y', values='value').astype(int)

g = sns.heatmap(
    heatmap_df, 
    vmax=4,
    vmin=-0.5,
    cmap='Blues',
    yticklabels=False,
    linewidths=2.5,
    cbar_kws={ 
        'ticks':[], 
        'boundaries' : np.arange(0,5,1),
        'orientation':'horizontal',
    },
    cbar_ax = cbar_ax,
    ax = ax
)

cbar_ax.get_xaxis().set_ticklabels([]);

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks(heatmap_df.columns.get_indexer(month_idx.index[1:]))
ax.set_xticklabels(month_idx.values[1:])

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

ax.set_title('Activity @ '+pd.Timestamp.today().strftime('%Y-%m-%d'))
ax.tick_params(top=False)

fig.savefig('Kaggle_Activity_' + pd.Timestamp.today().strftime('%Y') + '.png', bbox_inches='tight')
fig.savefig('Kaggle_Activity_' + pd.Timestamp.today().strftime('%Y') + '.svg', bbox_inches='tight')

plt.show()