#!/usr/bin/env python


import os
import sys
import numpy as np
import pandas as pd
import re

data = pd.read_csv(sys.argv[1], sep=',', engine='python', index_col=0)

dic_chrom = {
        'NC_015438.2':1,
        'NC_015439.2':2,
        'NC_015440.2':3,
        'NC_015441.2':4,
        'NC_015442.2':5,
        'NC_015443.2':6,
        'NC_015444.2':7,
        'NC_015445.2':8,
        'NC_015446.2':9,
        'NC_015447.2':10,
        'NC_015448.2':11,
        'NC_015449.2':12
        }

new = list()


counter = 0
for i in range(int(len(data))):
    counter += 1

row_counter = 0

for i in data.index:

    row_counter += 1

    if row_counter%10 == 0:
        print(row_counter, '/', counter)


    prot = i.split(';')
#    fold_change = data.loc[i, 'difference_log2']

    for j in range(int(len(prot))):
        cand = prot[j].strip()

        ID_temp = list()
        gene_temp = list()
        product_temp = list()
        chrom_temp = list()
        pos_temp = list()
        
        os.system('grep ' + cand + ' ' + sys.argv[2] + ' > temp_hits.txt')

        with open('temp_hits.txt') as f:
            for lines in f:
                gene = re.search(';gene=(.*);product',lines)
                product = re.search(';product=(.*);protein_id=', lines)
                chrom = re.search('^[A-Za-z0-9_\.]*', lines)
                pos_start = re.search('CDS\\t([0-9]*)', lines)
                pos_end = re.search('CDS\\t[0-9]*\\t([0-9]*)', lines)

                gene_temp.append(gene.group(1))
                product_temp.append(product.group(1))
                chrom_temp.append(chrom.group(0))
                pos_temp.append(pos_start.group(1))
                pos_temp.append(pos_end.group(1))

        if (len(set(gene_temp)) + len(set(product_temp))) > 2:
            print('ERROR: too many gene and product names found')
            break

        gene = list(set(gene_temp))[0]
        product = list(set(product_temp))[0]
        chrom = list(set(chrom_temp))[0]

        pos_temp = [int(x) for x in pos_temp]

        pos_start = min(pos_temp)
        pos_end = max(pos_temp)

        if chrom in dic_chrom:
            chrom = int(dic_chrom[chrom])

        else:
            chrom = chrom


        comb = list()

        comb.append(cand)
        comb.append(gene)
        comb.append(product)
        comb.append(chrom)
        comb.append(pos_start)
        comb.append(pos_end)
#        comb.append(fold_change)

        new.append(comb)

#df = pd.DataFrame(new, columns=['prot_id', 'gene', 'product_name', 'chromosome', 'low_bound', 'upp_bound', 'bion_fold2_change']) 
df = pd.DataFrame(new, columns=['prot_id', 'gene', 'product_name', 'chromosome', 'low_bound', 'upp_bound']) 

df.low_bound = df.low_bound.astype('int')
df.upp_bound = df.upp_bound.astype('int')


###################################################################################################
## Loop that looks duplicated loci in the data frame and keeps the locus with the longest product #
###################################################################################################

total = list()
curr = list()
curr.append(list())
curr.append(list())
curr.append(list())

prevID = None
counter = 0


for i in df.index:
    locus = df.iloc[i]['gene']
#    print(locus)

    if locus == prevID:
        print('1')
        curr[0].append(prevID)
        curr[1].append(prevL)
        curr[2].append(prevP)

        counter = 1

    elif locus != prevID:
        if counter == 1:
#            print(curr)
            counter = 0
            curr[0].append(prevID)
            curr[1].append(prevL)
            curr[2].append(prevP)
            print(curr)
            pos = curr[1].index(max(curr[1]))
            del curr[2][pos]

            for j in curr[2]:
                print(j)
                total.append(j)

            curr = list()
            curr.append(list())
            curr.append(list())
            curr.append(list())

    prevID = locus
    prevL = df.iloc[i]['upp_bound'] - df.iloc[i]['low_bound']
    prevP = df.iloc[i]['prot_id']


df = df.set_index('prot_id')
df2 = df.drop(total)


locus_info = df2[df2.chromosome.isin(list(range(12)))]

chrom_name = sys.argv[1][:-4] + '_chrom_total_info'
locus_info.to_csv(chrom_name)


locus_info = locus_info[['gene', 'chromosome', 'low_bound', 'upp_bound']]
#print(locus_info)
#locus_info = locus_info.drop_duplicates(subset=['gene'], keep='first')
#print(locus_info)


df2_name = sys.argv[1][:-4] + '_locus_info.csv'
df2.to_csv(df2_name)

locus_name = sys.argv[1][:-4] + '_locus_pipeline_selection.csv'
locus_info.to_csv(locus_name, index=False)

#print(df)


