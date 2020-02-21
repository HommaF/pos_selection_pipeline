#!/usr/bin/env python


import os
import sys
import numpy as np
import pandas as pd
import re

data = pd.read_csv(sys.argv[1], sep=',', engine='python', index_col=0)

print(data)

for i in data.index:

    print('prot ID:',i)
    print('gene name:',data.loc[i,'gene'])

    protein_file = data.loc[i, 'gene'] + '_prot_linear.fasta'
    os.system('grep -A 1 ' + i + ' ' + sys.argv[2] + ' > ' + protein_file) 


