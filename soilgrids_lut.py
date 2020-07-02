#!/usr/bin/env python

import pandas as pd
import numpy as np
from simpledbf import Dbf5

'''
- this script creates a look-up table for mHM soil classes derived by k-mean clustering soilgrids database
- should work for an arbitrary amount of cluster
- by Marco Hannemann UFZ ENINV, marco.hannemann@ufz.de'''


# dbf derived from k-mean clustering by Erik Nixdorf
dbf = Dbf5('data/soilgrids/soil_example.dbf')
df = dbf.to_dataframe()

# get number of cluster and create empty data frame
n_cluster = len([x for x in df['k_means'].unique() if ~np.isnan(x)])
df_cluster = pd.DataFrame(columns=df.columns[1:-1].to_list())

# calculate mean of parameters and fill data frame
for cluster in range(0, n_cluster):
    values = []
    for parameter in df.columns[1:-1]:
        values.append(df[parameter].loc[df['k_means'] == float(cluster)].mean())
    df_cluster.loc[cluster] = values

# convert units: sand/clay in g kg-1 -> %; bulk density in Mg m-3 -> g cm-3
bd_names = [name for name in df_cluster.columns if 'bdod' in name]
sand_names = [name for name in df_cluster.columns if 'sand' in name]
clay_names = [name for name in df_cluster.columns if 'clay' in name]

df_cluster[df_cluster.columns.difference(bd_names)] = df_cluster[df_cluster.columns.difference(bd_names)] / 10
df_cluster[bd_names] = df_cluster[bd_names] / 100

# add ID starting with 1 because soil class can not be zero
df_cluster['ID'] = list(range(1, len(df_cluster) + 1))

# generate columns for LUT
horizon_id = list(range(1, len(df_cluster) + 1)) * 6
soil_id = sorted(horizon_id)
upper_depth = [0, 5, 15, 30, 60, 100] * 6
lower_depth = [5, 15, 30, 60, 100, 200] * 6
clay = [x for depth in df_cluster[clay_names].values.tolist() for x in depth]
sand = [x for depth in df_cluster[sand_names].values.tolist() for x in depth]
bd = [x for depth in df_cluster[bd_names].values.tolist() for x in depth]

# create LUT for mHM and output txt
data = {'LBA10': soil_id,
        'HORIZON': horizon_id,
        'UD[mm]': upper_depth,
        'LD[mm]': lower_depth,
        'CLAY[%]': [int(round(x)) for x in clay],
        'SAND[%]': [int(round(x)) for x in sand],
        'Bd[gcm-3]': [round(x, 1) for x in bd]}

lut = pd.DataFrame(data)

dfStyler = lut.style.set_properties(**{'text-align': 'left'})
dfStyler.set_table_styles([dict(selector='th', props=[('text-align', 'left')])])

with open('soil_classdefinition.txt', "w+") as f:
    f.write(f'\tnSoil_Types  {len(df_cluster)}\n')

lut.to_csv('soil_classdefinition.txt', index=False, sep='\t', mode='a')
