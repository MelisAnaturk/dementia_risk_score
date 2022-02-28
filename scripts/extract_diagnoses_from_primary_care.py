
import csv
from typing import List, Any, Union
import os
import pandas as pd
#import seaborn as sns

# ------- 1. Specify pathways and read in csv
datapath = "../../raw_data/"
ICD910_df = pd.read_csv(datapath + 'gp_clinical.csv')
meds_list = pd.read_csv(datapath +'read_codes.csv')

ICD910_df.head(5)
col_names = pd.DataFrame(list(ICD910_df))
col_names = pd.DataFrame(list(meds_list))

# ------- 2. Create a loop to subset over all individuals
# create a loop
variable_list = ['Depression', 'Stroke', 'TBI','Diabetes', 'Diabetes_II', 'TIA', 'Atrial_fibrillation']

#list_to_extract = meds_list["%s_med_codes" % m].tolist()
for m in variable_list:
    ICD910_df = pd.read_csv(datapath + 'gp_clinical.csv')
    meds_list = pd.read_csv(datapath + 'read_codes.csv')
    list_to_extract = meds_list['%s_read_codes' % m].tolist()
    idx = pd.concat((ICD910_df[col].astype(str).str.startswith(tuple(list_to_extract)) for col in ['read_2', 'read_3']), axis=1).any(axis=1)
    ICD910_df = ICD910_df.loc[idx]
    print ('saving file containing participants with %s' %m)
    ICD910_df['idx'] = ICD910_df.groupby('eid').cumcount() + 1
    ICD910_df = ICD910_df.pivot_table(index=['eid'], columns='idx',
                        values=['event_dt'], aggfunc='first')
    ICD910_df = ICD910_df.sort_index(axis=1, level=1)
    ICD910_df.columns = [f'{x}_{y}' for x, y in ICD910_df.columns]
    ICD910_df = ICD910_df.reset_index()
    ICD910_df.to_csv(datapath+'participants_with_%s.csv' %m, sep=',',index=None)

