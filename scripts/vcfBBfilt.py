#! /usr/bin/python
#
# This script filters a vcf file based on overall sequence coverage, number of
# non-reference reads, number of alleles, and reverse orientation reads.
# See below for default values, and to change them, if necessary. Additionally,

# Usage: ./vcfBBffilt.py <variant file> <typls> > outfile.vcf
# NOTE: that typls has to be written in double quotes and in python list() format (e.g. "['DEL', 'INS', 'SUB']")

from sys import argv
import re
import shutil
from ast import literal_eval
#import tempfile

### stringency variables, edit as desired
minCoverage = 128 # minimum number of seqs; DP
minAltRds = 4 # minimum number of sequences with the alternative allele; AD
fixed = [0.0, 1.0] # removes loci fixed for alt; RAF
mapQual = 50 # minimum mapping quality
#typls = ["DEL", "INS", "SUB"]
typls = literal_eval(argv[2])


with open(argv[1], 'rb') as file:
    for line in file:
        line = bytes.decode(line).strip('\n')
        if line[0] == '#':
            #continue
            print(line)
        else:
            lspl = line.split('\t')
            test = lspl[6]
            altAll = lspl[4]
            if test == 'PASS':
                dp = int(re.findall('DP=[0-9]+', line)[0].split('=')[1])
                ac = int(re.findall('AD=[0-9]+', line)[0].split('=')[1])
                raf = float(re.findall('RAF=[0.0-9.0]+', line)[0].split('=')[1])
                mqs = int(re.findall('MQS=[0-9]+', line)[0].split('=')[1])
                mqm = int(re.findall('MQM=[0-9]+', line)[0].split('=')[1])
                ls = int(re.findall('LS=[0-9]+', line)[0].split('=')[1])
                typ = str(re.findall('TYP=[A-Z]+', line)[0].split('=')[1])
                ppc = int(re.findall('PPC=[0-9]+', line)[0].split('=')[1])
                #
                if (raf not in fixed and dp >= minCoverage and typ in typls and mqs >= mapQual and ac >= minAltRds):
                    print(line.strip('\n'))

    file.close()


