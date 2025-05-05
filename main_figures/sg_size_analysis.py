# Script for processing schoolgroup size checkpoint file output ('fl_schoolgroups.csv')
# Output file captures for each student/teacher: census tract, ethnicity, schoolgroup size 

import pandas as pd

import numpy as np

 

fl = pd.read_csv('fl_schoolgroups.csv')

ethnicity = []

census_tracts = []

sg_size = []

for index, row in fl.iterrows():

    if row['n_student'] > 0:

        sg_size_tmp = row['n_student'] + 1

        for i in range(0,row['n_non_hisp']):

            ethnicity.append(0)

        for i in range(0,row['n_hisp']):

            ethnicity.append(1)

        sg_size_tmp_eth = row['n_non_hisp'] + row['n_hisp']

        for i in range(0,sg_size_tmp_eth):

            census_tracts.append(row['tract_fips'])

            sg_size.append(sg_size_tmp)

 

output = pd.DataFrame()

output['ethnicity'] = ethnicity

output['census_tracts'] = census_tracts

output['sg_size'] = sg_size

output.to_csv('output.csv')
