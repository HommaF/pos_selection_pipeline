#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import glob
import os
import subprocess
import re
from scipy import stats
import string



tab = str.maketrans("ACTG", "TGAC")

def rev_compl_table(seq):
    return seq.translate(tab)[::-1]


#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna


# sys.argv[1]: directory with fasta files - all sequences must be within one line
# sys.argv[2]: yes or no: if yes, a phylogenetic tree will be made using raxml-ng


#####################################################################################
### read in files and run exonerate on them using cds and genomic regions ###########
#####################################################################################

proteins = list()
loci = list()

for files in glob.glob(sys.argv[1] + '/*prot_linear.fasta'):
    proteins.append(files)
proteins = sorted(proteins)

for files in glob.glob(sys.argv[1] + '/*genomic.fasta'):
    loci.append(files)
loci = sorted(loci)

if len(loci) != len(proteins):
    print('ERROR - NOT THE SAME NUMBER OF PROTEIN AND GENOMIC FILES FOUND')

lenght = len(proteins)

print(proteins)
print(loci)

for i in range(int(lenght)):
    genomic = loci[i]
    prt = proteins[i]

    exonerate_temp_file = prt[:-6] + '_exonerate_temp.fasta'
    exonerate_file = prt[:-6] + '_exonerate.fasta'
#    print('exonerate_file', exonerate_file)

    line_counter = 0

    with open(genomic) as f:
        new = list()

        for lines in f:
            lines = lines.strip()

#            print(lines)

            line_counter += 1

            if line_counter%2 == 1:
                seq_name = lines

            elif line_counter%2 == 0:
                seq_code = lines
                thefile = open('temp_seq.fasta', 'w')
                thefile.write('%s\n' % seq_name)
                thefile.write('%s\n' % seq_code)
                thefile.close()

                call = 'exonerate --model protein2genome --refine full --score 250 --bestn 1 --showvulgar no --showalignment FALSE --showsugar FALSE --showtargetgff TRUE ' + prt + ' temp_seq.fasta > ' + exonerate_temp_file
                os.system(call)

                os.system('cat ' + exonerate_temp_file + ' >> ' + exonerate_file)


#                call = 'exonerate --model protein2genome --refine full --score 250 --bestn 1 --showvulgar no --showalignment True --showsugar FALSE --showtargetgff TRUE ' + prt + ' temp_seq.fasta >> test2.txt'
#                os.system(call)

                counter = 2


                with open(exonerate_temp_file) as f:
                    for line in f:

                        line = line.strip('\r\n')
                        line = line.strip('\n')

                        if line == '# --- START OF GFF DUMP ---':
                            counter = 1
                            exons = list()

                        elif line == '# --- END OF GFF DUMP ---':
                            counter = 0

                        if counter == 1:
                            s_exon = re.search('.+exon\t([0-9]+)\t+([0-9]+)', line)
                            s_name = re.match('(.+)(\texonerate)', line)

                            if s_exon != None:
#                                print(s_name.group(1))
#                                print(s_exon.group(2))
                                exons.append(s_exon.group(1))
                                exons.append(s_exon.group(2))
                                s_rev = re.match('.+exon\t([0-9\t\.]+)([-,+]+)', line)

                        elif counter == 0:

                            counter = 2

                            exons = [int(x) for x in exons]
#                            print(exons)

                            if s_rev.group(2) == '+':

                                for j in range(int(len(exons)/2)):
#                                    print('here we go!!')
#                                    print(s_rev.group(2))
#                                    print('seq_code_exon:',seq_code[exons[j*2]-1:exons[j*2+1]])

                                    if j == 0:
                                        seq_exon = seq_code[exons[j*2]-1:exons[j*2+1]]
#                                        print(seq_exon)

                                    else:
                                        seq_exon += seq_code[exons[j*2]-1:exons[j*2+1]]
#                                        print(seq_exon)

                                new.append(seq_name)
                                new.append(seq_exon)

                            elif s_rev.group(2) == '-':
                                for j in range(int(len(exons)/2)):
                                    rev_seq = seq_code[exons[j*2]-1:exons[j*2+1]].upper()
                                    rev_compl_seq = rev_compl_table(rev_seq)

#                                    print('rev seq:',rev_seq)
#                                    print('rev compl sequence:',rev_compl_seq)

                                    if j == 0:
                                        seq_exon = rev_compl_seq

                                    else:
                                        seq_exon += rev_compl_seq

                                new.append(seq_name)
                                new.append(seq_exon)


        call = 'rm ' + exonerate_temp_file
        os.system(call)

        name = exonerate_file[:-15] + 'multi_sl25_prt.fasta'

        thefile = open(name, 'w')
        for item in new:
            thefile.write('%s\n' % item)
        thefile.close()


#                            for j in range(int(len(exons)/2)):
#                                print('low_bound', exons[j*2])
#                                print('upp_bound', exons[j*2+1])
#                                print('seq_name: ',seq_name)
#                                print('seq_code_total:',seq_code)
#                                print('exons', exons)
#                                print('seq_code_exon:',seq_code[exons[j*2]-1:exons[j*2+1]])


#####################################################################################
### extracts all the cds from the regions generated with a previous pipeline ########
### removes all end of line characters, looks for the name of the sequence, #########
### reads in the indicator line 'indic' to know which nucleotides of the subject ####
### alignment to integrate in the new sequence representing the accession ###########
#####################################################################################

    skip = 0
    query = 0

    counter = 0



'''

# this should be in comments


    with open(exonerate_file) as f:

        for line in f:
            line = line.strip('\r\n')
            line = line.strip('\n')


#            print(line)


            if line[8:16] == 'Target: ':

                query = 0
                indic_h = 0

                seq_name = '>' + line[16:]
                sequence = ''

                counter += 1

#                print(seq_name)

            elif query == 0:
#                if line[6:7] == ':':
                if re.match('[0-9 ]+: [a-zA-Z]*', line):
                    spacer = len(re.match('([0-9 ]+: )[a-zA-Z]*', line).group(1))
                    query = 1

            elif query == 1:

                if line[:15] != '         Query:':

#                    if re.match('[ \:\!]{5,}[\|]+', line):
                    if re.match('[ ]{5,}', line[:spacer]):
                        if indic_h == 0:

#                        indic = str(re.match('[ ]+([\|\!\:\. ]*)', line).group(1))
                            indic = line[spacer:]

                            indic_h = 1

                        elif indic_h == 1:
                            helper = line[spacer:]
                            indic_h = 0

                    elif re.match('[0-9 ]+: [a-zA-Z]*', line):
                        query = 0

                        target = re.match('[ ]+[0-9]+ : (.*)', line)
                        target = target.group(1)

#                        print(indic)
#                        print(helper)
#                        print(target)

                        for i in range(len(indic)):
                            if indic[i:i+1] == '|':
#                                if indic[i+1:i+2] != '}':
                                if '}' not in target[i:i+3]:
                                    if '{' not in target[i-2:i+1]:
                                        sequence += target[i:i+1]

                            elif indic[i:i+1] == ' ':
                                if target[i:i+1] != '.':
                                    if helper[i:i+1] != '+':
                                        if helper[i:i+1] != '-':
#                                            if '}' not in target[i-2:i+1]:
#                                                if '{' not in target[i:+3]:
                                            if target[i:i+1] != '-':
                                                sequence += target[i:i+1]

                            elif indic[i:i+1] == '!':
                                sequence += target[i:i+1]

                            elif indic[i:i+1] == ':':
                                sequence += target[i:i+1]

                            elif indic[i:i+1] == '.':
                                sequence += target[i:i+1]


            if line[:2] == '--':

                if skip == 0:
                    skip = 1

                elif skip == 1:

#                    if reverse == 1:
#                        sequence = Seq(sequence, generic_dna)
#                        print(sequence)
#                        sequence = sequence.complement()
#                        print(sequence)
    #                print(sequence)
                    new.append(seq_name)
                    new.append(sequence)
                    seq_name = ''
                    sequence = ''

                    skip = 0

    name = exonerate_file[:-15] + 'multi_sl25_prt.fasta'

    print('name after exonerate analysis:',name)

    thefile = open(name, 'w')
    for item in new:
        thefile.write('%s\n' % item)
    thefile.close()

'''

#####################################################################################
### reorganise all files and move the new cds files into the new 'cds' directory ####
### define all accepted nucleotides and all stop codons. Stop codons are all in #####
### small letters since the strings later will also be in small letters #############
#####################################################################################



os.system('rm temp_seq.fasta')
os.system('mkdir ./prt')
#os.system('ll ' + sys.argv[1] + '*multi_sl25_prt.fasta')
os.system('mv ' + sys.argv[1] + '*multi_sl25_prt.fasta ' + './prt') 


deleted = list()

for files in glob.glob('./prt/*fasta'):

    statinfo = os.stat(files)

    if statinfo.st_size == 0:
        print('the following file ' + str(files) + ' was removed due to his size of: ' + str(statinfo.st_size))
        deleted.append(files)
        os.remove(files)


nucleotides = set('AGCTagct_-')
#stop_codons = ("TAG", "tag", "TAA", "taa", "TGA", "tga")
stop_codons = ("tag", "taa", "tga")

log_removed = list() 
log_wrong = list()

os.system("mkdir no_frameshift")


print("removing pseudogenes")

hit = 0

w_dir = './prt'

#####################################################################################
### removes all sequences w/ intenal stop codons, frameshifts and sequences #########
### containing non-normal nucleotides ###############################################
### all terminal stop codons are removed for downstream analysis ####################
#####################################################################################


for files in glob.glob(w_dir + '/*.fasta'):

    print(files)

    helper = files[5:-6]

    helper = helper + "_noFrameshift.fasta"
    print('helper', helper)

    log_removed.append(helper)
    log_wrong.append(helper)   

    genes = list()

    with open(files) as f:
        for line in f:
#            print('line read in:\n',line)
            line = line.strip()


            if line[:1] == '>':
#                print(line)
                genes.append(line)
                log_wrong.append(line)
                seq_name = line

                hit = 0

            elif line[:1] != '>':
#                print('overview:\n', line)

# check for non-regular bases in nt sequence

                if set(line).issubset(nucleotides):
                    log_wrong = log_wrong[:-1]
#                    print('characters_correct:\n', line)

# remove terminal stop codons
                    
                    end = line[-3:]
                    if end.lower() in stop_codons:
                        line = line[:-3]
#                        print('term stop removed:\n:', line)
                    

                    counter = 0
                    valid = 2

                    for i in range(len(line)):
                        trimmed = line[i:]

                        a = trimmed[:1]
                        b = trimmed[1:2]
                
                        if (i == len(line)-1):
                            if (a == '_'):
                                counter += 1

                                if (counter%3 != 0):
                                    genes = genes[:-1]        
                                    valid = 0

                                else:
                                    valid = 1
 



                        elif (a == '_' and b == '_'):
                            counter += 1
                            
                        elif (a == '_' and b != '_'):
                            counter += 1
                            
                            if (counter%3 != 0):
                                valid = 0
                                genes = genes[:-1]

                                break

                            else:
                                valid = 1

                        
# remove all sequences that have internal stop codons

                    nbCodons = int(len(line)/3)
                    internalCounter = 1

                    for i in range(nbCodons):
                        trimmed = line[i*3:]
                        a = trimmed[:3]

                        if a.lower() in stop_codons:
                            internalCounter = 0
                            i = nbCodons-1


                        if (i == nbCodons-1):
#                            if hit ==1:
#                                print('internalCounter')
#                                print(internalCounter)
#                                print('counter')
#                                print(counter)
#                                print('valid')
#                                print(valid)
#                            if (internalCounter == 1 and counter%3 == 0):
#                                valid = 1

                            if (valid > 0 and internalCounter == 1):
                                genes.append(line)
#                                if hit == 1:
#                                    print('made it')
#                                    print(genes)

                            elif (valid > 0 and internalCounter == 0):
                                genes = genes[:-1]
                                log_removed.append(seq_name)
                                log_removed.append(line)
                            
                            break

                else: 
                    log_wrong.append(line)
                    genes = genes[:-1]

# writing log file with removed sequences and the cleaned fasta files without terminal stop codons

    thefile = open("./no_frameshift/" + helper, 'w')
    for item in genes:
        thefile.write("%s\n" % item)
    thefile.close()

thefile = open("./no_frameshift/log_stops.txt", 'w')
for item in log_removed:
    thefile.write("%s\n" % item)
thefile.close()

thefile = open("./no_frameshift/log_wrong_characters.txt", 'w')
for item in log_wrong:
    thefile.write("%s\n" % item)
thefile.close()




# moving files into respective directories

os.system("mkdir translated")

print("translating sequences for msa")

# translation of sequences using transeq
#####################################################################################
### running transeq to translate the curated cds ####################################
#####################################################################################


os.chdir("./no_frameshift")
 
for files in glob.glob("*.fasta"):
    print(files)
    helper = files[:-6]
    helper = helper + "_translated.fasta"

    os.system("transeq " + files + " ../translated/" + helper)



# making msa with amino acid sequences
#####################################################################################
### running mafft to generate an aa alignment #######################################
#####################################################################################

print("making msa alignments")

os.chdir("../translated")
os.system("mkdir ../aa_alignment")

for files in glob.glob("*.fasta"):
    helper = files[:-6]
    helper = helper + "_aa_alignment.fasta"

    os.system("mafft --quiet " + files + " > ../aa_alignment/" + helper)

os.chdir("../aa_alignment")


for files in glob.glob("../no_frameshift/*.fasta"):
    os.rename(files, files[17:])
    
# make protein-guided codon alignment 


print("generate codon alignemnt")

os.system("mkdir ../codon_alignments")

# read fasta files into a list

array = list()

for files in sorted(glob.glob("*.fasta")):
    array.append(files)
#    print(array)


# check even number of fasta files for QC

print(array)
len_array = len(array)
checkpoint = len_array%2
rounds = int(len_array/2)

if checkpoint != 0:
    print("ERROR - not the same amount of files indicated for protein guided alignment!")
    


# run pal2nal.pl script for all pairs of arrays - it is assumed that the order of files is always both files that belong together are sequentially and that the nt sequence comes firts - should be ensured due to the systematic naming approach



# call command in os.systems using a variable that contains everything - modify so that it works every time

#####################################################################################
### run pal2nal to make codon-guided msa to combine cds w/ aa alignment #############
#####################################################################################


nt = 0
aa = 1

for i in range(rounds):
    output = array[nt][:-19]
    output = output + "_codonAlignment.fasta"

    
    call = "pal2nal.pl " + array[aa] + " " + array[nt] + " -output fasta > ../codon_alignments/" + output

    os.system(call) 
    nt = nt + 2
    aa = aa + 2


# remove duplicates using raxml

#####################################################################################
### run raxml-ng with codon-guided msa for later paml analysis ######################
#####################################################################################

for files in glob.glob('../codon_alignments/*fasta'):

    statinfo = os.stat(files)

    if statinfo.st_size == 0:
        print('the following file ' + str(files) + ' was removed due to his size of: ' + str(statinfo.st_size))
        deleted.append(files)
        os.remove(files)




if sys.argv[2] == 'yes':

    print("remove duplicates")
    os.system("mkdir ../raxml-ng")
    os.system("mkdir ../raxml-ng/total_length")
    os.system("cp ../codon_alignments/*.fasta ../raxml-ng")
    os.chdir("../raxml-ng")

    for files in glob.glob("*.fasta"):
        cropped = files[:-6]
        call = "raxml-ng --msa " + files + " --prefix " + cropped + " --model GTR+G --seed 123 --bs-tree 100 --force" 
        os.system(call)
        os.system("mv " + files + " ./total_length")


#####################################################################################
### tidy up all files in a summary called 'codon_guided_msa' ########################
#####################################################################################


os.chdir("..")
os.system("mkdir codon_guided_msa")
os.system("mv aa_alignment ./codon_guided_msa")
os.system("mv codon_alignments ./codon_guided_msa")
os.system("mv no_frameshift ./codon_guided_msa")
os.system("mv translated ./codon_guided_msa")
os.system("mv raxml-ng ./codon_guided_msa")
os.system("mv prt ./codon_guided_msa")

os.system("mkdir exonerate_output")
os.system("mv *" + sys.argv[1] + "*exonerate* " + "./exonerate_output")


# extract information of which subsitution model is best suited - count occurence


for files in glob.glob('./codon_guided_msa/codon_alignments/*'):
    print(files)
    name_msa = re.match('./codon_guided_msa/codon_alignments/([a-zA-Z0-9_]*)_prot_linear_multi', files)
    print(name_msa.group(1))

    hit = 0

    for trees in glob.glob('./codon_guided_msa/raxml-ng/*bestTree'):
        name_tree = re.match('./codon_guided_msa/raxml-ng/([a-zA-Z0-9_]*)prot_linear_multi', trees)
        print(trees)
        print(name_tree.group(1))

        if name_msa.group(1) in name_tree.group(1):
            hit = 1

    if hit == 0:
        print('the following file', str(files), 'was removed due to no matching phylogenetic tree')
        deleted.append(files)
        os.remove(files)


thefile = open('./codon_guided_msa/no_frameshift/empty_files_removed.log', 'w')
for item in deleted:
    thefile.write('%s\n' % item)
thefile.close()


#####################################################################################
### read in codon-guided alignemts and trees to make control files for paml #########
#####################################################################################

msa_list = list()
tree_list = list()
print('here we go')

os.system('mkdir model_testing_CODEML')

# change so that they're read in separatey as fasta and tree files

for files in glob.glob('./codon_guided_msa/codon_alignments/*fasta'):
    msa_list.append(files[36:])

for files in glob.glob('./codon_guided_msa/raxml-ng/*bestTree'):
    tree_list.append(files[28:])

msa_list = sorted(msa_list)
tree_list = sorted(tree_list)
length = int(len(msa_list))

#print('msa_list',msa_list)
#print('tree_list', tree_list)


for i in range(length):

    msa = './codon_guided_msa/codon_alignments/' + msa_list[0+i]
    tree = './codon_guided_msa/raxml-ng/' + tree_list[0+i]

    print('msa: ', msa)
    print('msa[36:-47]: ', msa[36:-47])
    print('tree ', tree)
#    print('tree ', tree)


    os.system('mkdir model_testing_CODEML/' + msa[36:-47])

#####################################################################################
### read in files and run exonerate on them using cds and genomic regions ###########
#####################################################################################


    models = (1, 2, 7, 8)
#    models = (1, 2, 0)
    for model in models:

        print('model to be used:',model)

        codeml_ctl = list()

        with open("codeml_default.ctl") as f:
            for line in f:
                line = line.strip()

                if line[:7] == "seqfile":
                    line = "seqfile = " + msa
                    codeml_ctl.append(line)

                elif line[:8] == "treefile":
                    line = "treefile = " + tree
                    codeml_ctl.append(line)

                elif line[:7] == "outfile":
                    line = "outfile = mlc_" + str(model) + "_" + msa_list[0+i][:-5] + "txt"
                    codeml_ctl.append(line)

                elif line[:7] == "NSsites":
                    line = "NSsites = " + str(model)
                    codeml_ctl.append(line)

                else:
                    codeml_ctl.append(line)

        name = str(msa[36:-47]) + '_' + str(model) + '_codeml.ctl'
        thefile = open(name, "w")
        for item in codeml_ctl:
            thefile.write("%s\n" % item)
        thefile.close()

#        os.system("cat codeml.ctl")
#        subprocess.call(["codeml"])

#    os.system("mv " + data[0+i*2][:-5] + 'txt' + "*" + " model_testing_CODEML/" + msa[36:-47])
#    os.system("mv mlc* model_testing_CODEML/" + msa[36:-47])

#os.system('ls > zzzz_9999_codeml.ctl')
#os.system('rm *_0_codeml.ctl')

for files in glob.glob('*codeml.ctl'):
    print('file to be used for analysis:',files)
    os.system('codeml ' + files)

#os.system('rm zzzz_9999_codeml.ctl')
os.system('rm 2*')
os.system('rm lnf')
os.system('rm rst*')
os.system('rm rub')


os.system('grep lnL mlc* > modeltesting.txt')

print(models)




info = list()
new = list()

with open('modeltesting.txt') as f:

    for lines in f:
        model = re.search('^mlc_([0-9])_', lines)
        name = re.search('^mlc_[0-9]_(.*)_prot_linear_', lines)
        df = re.search('np: ([0-9]*)\):', lines)
        ln = re.search('\): +([0-9\.\-]*) ', lines)
        test = re.search(":\t([0-9\.\-]*)\t", lines)


        print('name', name)
        print('model', model)
        print('df', df)
        print('ln', ln)


        if name != None:

            print('name:', name.group(1))
            print('model:', model.group(1))
            print('df:', df)
            print('ln:', ln.group(1))


            info.append(str(name.group(1)))
            info.append(int(model.group(1)))
            info.append(int(df.group(1)))
            info.append(float(ln.group(1)))
#            name_list.append(str(name.group(1)))
#            model_list.append(int(model.group(1)))


for i in range(int(len(info)/4)):
    new.append(info[0+4*i:4+4*i])

df_sum = pd.DataFrame(new, columns=['name', 'model', 'df', 'logL'])
df_sum = df_sum.sort_values(by=['name', 'model'])
print(df_sum)

comp = list()
new = list()

for i in df_sum.name.unique():

    df_temp = df_sum[df_sum.name == i]

    ln1 = df_temp.logL[df_temp.model == 1].item()
    ln2 = df_temp.logL[df_temp.model == 2].item()
    ln7 = df_temp.logL[df_temp.model == 7].item()
    ln8 = df_temp.logL[df_temp.model == 8].item()
#    print('ln:',ln1, ln2, ln7, ln8)

    df1 = df_temp.df[df_temp.model == 1].item()
    df2 = df_temp.df[df_temp.model == 2].item()
    df7 = df_temp.df[df_temp.model == 7].item()
    df8 = df_temp.df[df_temp.model == 8].item()
#    print('df:',df1, df2, df7, df8)

    dln1_2 = ln2 - ln1
    dln7_8 = ln8 - ln7

#    print('dln:', dln1_2, dln7_8)

    ddf1_2 = df2 - df1
    ddf7_8 = df8 - df7

#    print('ddf:', ddf1_2, ddf7_8)

#    print(dln1_2)
#    print(ddf1_2)

    pchisq1_2 = stats.chi2.sf(2*dln1_2, ddf1_2)
    pchisq7_8 = stats.chi2.sf(2*dln7_8, ddf7_8)

#    print(pchisq1_2)
#    print(pchisq7_8)

    comp.append(i)
    comp.append(pchisq1_2)
    comp.append(pchisq7_8)

#    print(comp)

for i in range(int(len(comp)/3)):
    new.append(comp[0+3*i:3+3*i])
#    print('test')

df_ana = pd.DataFrame(new, columns=['name', 'M1-M2', 'M7-M8'])
print(df_ana)

df_sum.to_csv('summary_selection_models.csv', index=None)
df_ana.to_csv('summary_model_testing.csv', index=None)

os.system('mv *codeml.ctl ./model_testing_CODEML/')
os.system('mv *codonAlignment.txt ./model_testing_CODEML/')
