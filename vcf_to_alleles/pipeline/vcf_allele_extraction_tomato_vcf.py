#!/usr/bin/env python


# python script to extract genome region of interest with inserted SNPs or indels from a vcf file

import sys
import glob
import pandas as pd
import numpy as np
import random
import time


## include QC filter and how to deal with several alleles

# sys.argv[1]: file with genes of interest
# sys.argv[2]: genome of interest
# sys.argv[3]: path to directory with vcf files

data = pd.read_csv(sys.argv[1], sep=",", engine="python")

# nb of vcf files:

counter = 0
for nb_vcf in glob.glob(sys.argv[3] + "*vcf"):
    counter += 1 

print(counter)

# looping through all genes of interest

for i in range(len(data)):

    finish_counter = 0

    chrom_counter = 0
    print("gene of interest: " + str(data.iloc[i][0]))
    chrom_oi = int(data.iloc[i][1])
    lower_bound = int(data.iloc[i][2])
    upper_bound = int(data.iloc[i][3])
    gene_name = data.iloc[i][0]

    gene_variants = list()

# genome of interest is opened and region of genome with gene of interest extracted

    with open(sys.argv[2]) as f:
        for line in f:
            line = line.strip()
            if line[:1] == ">":

                if chrom_counter < 10:

                    chrom_counter += 1

                elif chrom_counter > 9:

                    chrom_counter += 1

# matching the chromosome of interest

            elif chrom_counter == chrom_oi:

                print('current chromosome: ' + str(chrom_counter))
                
# defining the name of the extracted reference sequence

                if chrom_counter < 10:
                    target_chrom = ">" + "Sl_" + data.iloc[i][0] + "_Chromosome_0" + str(data.iloc[i][1])

                else:
                    target_chrom = ">" + "Sl_" + data.iloc[i][0] + "_Chromosome_" + str(data.iloc[i][1])


                region = line[lower_bound-500:upper_bound+500]
                backup = region

                print(target_chrom)
                print(region)

                gene_variants.append(target_chrom)
                gene_variants.append(region)


# opening all vcf files to extract SNP info for region of interest

    for vcf in glob.glob(sys.argv[3] + '*.vcf'):

        region = backup
        genome_name = vcf[14:-4]
        seq_name = ">" + gene_name + "_" + genome_name
        gene_variants.append(seq_name)

        linear_line = None
        snp_stats = list()
        lin = list()

        chrom_counter = 0
        prog_counter = 0
        pos_counter = 0
        line_counter = 0

# determine number of comments lines to skip in vcf file

        with open(vcf) as v:

            lines_of_comments = 0

            for line in v:
                if line[:6] != '#CHROM':
                    lines_of_comments += 1

                elif line[:6] == '#CHROM':
                    break

# reading in the current vcf file
        if chrom_oi > 9:
            clip_chrom = 'SL2.50ch' + str(chrom_oi)

        else:
            clip_chrom = 'SL2.50ch0' + str(chrom_oi)


        snps = pd.read_csv(vcf, sep="\\t", skiprows = lines_of_comments, engine = "python")
#        snps = snps.loc[(snps['#CHROM'] == clip_chrom) & (snps['POS'].between(lower_bound-500, upper_bound+500)) & (snps['QUAL'] > 80)]
        snps = snps.loc[snps['#CHROM'] == clip_chrom]
        snps = snps.loc[snps['POS'] > lower_bound-500]
        snps = snps.loc[snps['POS'] < upper_bound+500]
        snps = snps.loc[snps['QUAL'] > 80]
        snps = snps.reset_index()
        current_genome = vcf[:-4]

        seq = upper_bound - lower_bound + 1000
        snp_range = range(0,seq+1)
        nb_snps = 0
        nb_counting_snps = 0
        indel_correction = 0


# adjusting the pointer in the vcf file to the right chromosome

        if len(snps) == 0:
            gene_variants = gene_variants[:-1]
            finish_counter += 1

        else:
            
            for j in range(len(snps)):

                ref = str(snps['REF'][j])
                len_ref = len(str(ref))
                before = snps['POS'][j] - (lower_bound-500) - 1 + indel_correction
                after = before + len(str(ref))

                alt = random.choice(str(snps['ALT'][j]).split(","))
                len_alt = len(alt)
                
                if region[before:after].lower() == ref.lower():
                            
                    region = region[:before] + alt + region[before+len_ref:]
                    indel_correction = indel_correction + len_alt - len_ref
                    nb_counting_snps += 1
                    nb_snps += 1

                   
                if j == len(snps) - 1:
                    gene_variants.append(region)
                    finish_counter += 1

                    if finish_counter%10 ==0:
                        print('finished vcf', finish_counter, '/', counter)
                        print(time.process_time())

    thefile = open(gene_name + "_genomic.fasta", 'w')
    for item in gene_variants:
        thefile.write("%s\n" % item)
    thefile.close()
    print(time.process_time())


        






#                for vcf in glob.glob(sys.argv[3] + "*.vcf"):
#                    snps = pd.read_csv(vcf, sep="/\t", engine="python")
#                    print(snps)



#                hits.append(">" + target_chrom)
#                hits.append(line[lower_bound-500:upper_bound+500])

# genome_directory = "/" + sys.argv[1] + "*.fasta"
#print(genome_directory)

'''

chrom_oi = int(sys.argv[2])
lower_bound = int(sys.argv[3])
upper_bound = int(sys.argv[4])

for genome in glob.glob("./genomes/*fasta"):
    print(genome)

    with open(genome) as f:
        for line in f:
            line = line.strip()
            if line[:1] == ">":

                if chrom_counter < 10:

                    chrom_counter += 1

                    if chrom_counter == chrom_oi:
                        target_chrom = line


                elif chrom_counter > 9:

                    chrom_counter += 1

                    if chrom_counter == chrom_oi:
                        target_chrom = line


            elif chrom_counter == chrom_oi:

                print('current chromosome: ' + str(chrom_counter-1))
    #            print(target_chrom)
                
                print(len(line))
                print('found it')
                print(line[lower_bound:upper_bound])

                hits.append(">" + target_chrom)
                hits.append(line[lower_bound-500:upper_bound+500])

#print(hits)

file_name = "mfasta_" + sys.argv[1] + "_nt.fasta"

thefile = open(file_name, 'w')
for item in hits:
    thefile.write("%s\n" % item)




'''
# 
#                        print(alt)
#                        for j in range(len(alt)):
#                            region = region[:before] + alt[j] + region[before+len_ref:]

