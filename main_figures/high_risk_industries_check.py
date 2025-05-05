# Script to print high risk industries (NAICS codes)

import pandas as pd
import numpy as np

work_contact = pd.read_csv('worker_contact_scalers_expanded.csv')
workgroups = pd.read_csv('workgroups.csv')
workgroups = workgroups.groupby('NAICS').mean()

work = workgroups.merge(work_contact,on='NAICS')

work['combined'] = np.log(work['workgroup size']) * work['public_exposure'] * work['proximity']

high_risk = np.quantile(work['combined'],0.9)

print([int(h) for h in list(work['NAICS'][work['combined'] > high_risk])])
