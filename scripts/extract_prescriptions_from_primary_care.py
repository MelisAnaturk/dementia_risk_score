# script written by Melis Anaturk (June 2020)
# contact: melis.anaturk@psych.ox.ac.uk

import csv
from typing import List, Any, Union
import os
import pandas as pd
#import seaborn as sns

# ------- 1. Specify pathways and read in csv
datapath = "../../raw_data"

#meds_df.head(5)
#col_names = pd.DataFrame(list(meds_df))
#col_names = pd.DataFrame(list(meds_list))

# ------- 2. Create a loop to subset over all individuals
# create a loop
variable_list = ['dementia', 'HRTs', 'statins', 'NSAIDs', 'antidepressants', 'Diabetes_II']
#list_to_extract = meds_list["%s_med_codes" % m].tolist()
for m in variable_list:
    meds_df = pd.read_csv(datapath + 'gp_scripts.csv')
    meds_list = pd.read_csv(datapath + 'prescription_list.csv')
    list_to_extract = meds_list['%s_med_codes' % m].tolist()
    idx = pd.concat((meds_df[col].astype(str).str.startswith(tuple(list_to_extract)) for col in ['read_2', 'bnf_code', 'dmd_code', 'drug_name']), axis=1).any(axis=1)
    meds_df = meds_df.loc[idx]
    print ('saving file containing participants with %s' %m)
    meds_df['idx'] = meds_df.groupby('eid').cumcount() + 1
    meds_df = meds_df.pivot_table(index=['eid'], columns='idx',
                        values=['issue_date','drug_name'], aggfunc='first')
    meds_df = meds_df.sort_index(axis=1, level=1)
    meds_df.columns = [f'{x}_{y}' for x, y in meds_df.columns]
    meds_df = meds_df.reset_index()
    meds_df.to_csv(datapath+'participants_%s.csv' %m, sep=',',index=None)


