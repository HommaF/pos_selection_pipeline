#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import sys
import glob
import re


# Script to identify residues under selection in a protein of interest and look up the alternative residues

dic_codons = {
        'GCT' : 'A',
        'GCC' : 'A',
        'GCA' : 'A',
        'GCG' : 'A',
        'CGT' : 'R',
        'CGC' : 'R', 
        'CGA' : 'R', 
        'CGG' : 'R', 
        'AGA' : 'R',
        'AGG' : 'R',
        'AAT' : 'N',
        'AAC' : 'N',
        'GAT' : 'D',
        'GAC' : 'D',
        'TGT' : 'C',
        'TGC' : 'C',
        'GAA' : 'E',
        'GAG' : 'E',
        'CAA' : 'Q',
        'CAG' : 'Q',
        'GGT' : 'G',
        'GGC' : 'G',
        'GGA' : 'G',
        'GGG' : 'G',
        'CAT' : 'H',
        'CAC' : 'H',
        'ATT' : 'I',
        'ATC' : 'I',
        'ATA' : 'I',
        'CTT' : 'L',
        'CTC' : 'L',
        'CTA' : 'L',
        'CTG' : 'L',
        'TTA' : 'L',
        'TTG' : 'L',
        'AAA' : 'K',
        'AAG' : 'K',
        'ATG' : 'M',
        'TTT' : 'F',
        'TTC' : 'F',
        'CCT' : 'P',
        'CCC' : 'P',
        'CCA' : 'P',
        'CCG' : 'P',
        'TCT' : 'S',
        'TCC' : 'S',
        'TCA' : 'S',
        'TCG' : 'S',
        'AGC' : 'S',
        'AGT' : 'S',
        'ACT' : 'T',
        'ACC' : 'T',
        'ACA' : 'T',
        'ACG' : 'T',
        'TGG' : 'W',
        'TAT' : 'Y',
        'TAC' : 'Y',
        'GTT' : 'V',
        'GTC' : 'V',
        'GTA' : 'V',
        'GTG' : 'V',
        'TAA' : 'STOP',
        'TAG' : 'STOP',
        'TGA' : 'STOP',
        '---' : 'GAP'
        }

mlc_files = list()
alignments = list()

for files in glob.glob(sys.argv[1] + '*'):
    if files[-3:] == 'txt':
        mlc_files.append(files)

    elif files[-5:] == 'fasta':
        alignments.append(files)

mlc_files.sort()
alignments.sort()

print(mlc_files)
print(alignments)
    
new_file = list()

for h in range(int(len(alignments))):

    start = 0
    new = list()
    loci = list()

    name = re.match('.+\/(.+)_prot', alignments[h]).group(1)
    print(name)

    with open(mlc_files[h]) as f:
        for line in f:

            begin = re.match('Bayes Empirical Bayes \(BEB\)', line)
            end = re.match('The grid', line)

            if begin != None:
                start = 1

            if end != None:
                start = 0

            if start == 1:

#                sign = re.match(' +[0-9]+ [A-Z]{1} +[0-9\.]{5}(\*+) ', line)
                sign = re.match(' +[0-9]+ [A-Z]{1} +([0-9\.]{5})', line)
                if sign != None:
#                    print('its not none:', sign.group(1))

                    if float(sign.group(1)) >= 0.8:

                        print('its larger than 0.8:', sign.group(1))
#                        print('it might have an asteriks:', sign.group(2))

                        entry = list()

                        pos_aa = int(re.match(' +([0-9]+) +[A-Z]{1}?', line).group(1))
                        ref_aa = re.match(' +[0-9]+ {1}?([A-Z]{1}?)', line).group(1)
                        sign_lv = float(sign.group(1))

                        loci.append(pos_aa-1)

                        entry.append(pos_aa)
                        entry.append(ref_aa)
                        entry.append(sign_lv)

                        new.append(entry)


    print(new)
    print(loci)

    alt_seq = list()
    for i in range(int(len(loci))):
        alt_seq.append(list())

    codon_list = list()
    for i in range(int(len(loci))):
        codon_list.append(list())
        
    header = list()
    ref_seq = list()

    ref_hit = 0

    with open(alignments[h]) as f:
        for line in f:
            line = line.strip()

            if ref_hit < 2:
                ref_hit += 1

            if ref_hit == 2:
                ref_hit = 3

                gap_counter = 0



                for i in range(int(len(line)/3)):

                    codon = line[0+3*i:3+3*i]

                    if codon  == '---':
                        gap_counter += 1

                    elif i in loci:
    #                    print(str(i+1) + ' (-' + str(gap_counter) + ')')
    #                    print(dic_codons[codon.upper()])
                        
                        if gap_counter == 0:
                            header.append(str(i+1))

                        elif gap_counter != 0:
                            header.append(str(i+1) + '(-' + str(gap_counter) + ')')


                        ref_seq.append(dic_codons[codon.upper()])

            elif (ref_hit == 3) & (line[:1] != '>'):


                counter = 0

                for i in range(int(len(line)/3)):

                    codon = line[0+3*i:3+3*i]

                    if i in loci:
#                        print('this is i:', i)
#                        print('codon:', codon)
#                        print('counter', counter)

                        if dic_codons[codon.upper()] not in alt_seq[counter]:
                            if dic_codons[codon.upper()] != ref_seq[counter]:
                                alt_seq[counter].append(dic_codons[codon.upper()])

                        elif dic_codons[codon.upper()] in alt_seq[counter]:
                            if codon.upper() not in codon_list[counter]:
                                codon_list[counter].append('*')


                        if codon.upper() not in codon_list[counter]:
                            if dic_codons[codon.upper()] != ref_seq[counter]:
                                codon_list[counter].append(codon.upper())

                        counter += 1

    print(header)
    print(ref_seq)
    print(alt_seq)
    print(codon_list)

    final = list()

    for i in range(int(len(header))):

        temp = list()
        ref_aa_temp = None

        temp.append(ref_seq[i])
        temp.append(header[i])
        temp.append(''.join(alt_seq[i]))

        final.append(''.join(temp))

    final = '; '.join(final)
    temp = list()

    temp.append(name)
    temp.append(final)

    new_file.append(temp)

data = pd.DataFrame(new_file, columns=['seq_name', 'residues under positive selection'])

print(data)

data.to_csv(sys.argv[2], index=None)


