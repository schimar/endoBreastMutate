#! /usr/bin/python

# This script reads through a vcf file and writes a bimbam file for all individuals and the given genotypes, with '0' for ref homozygote, '1' for heterozygote, and '2' for alt homozygote.

# Usage: ~/vcf2bimbam.py file.vcf > file.mg


from sys import argv
import numpy as np
import pandas as pd


def get012(x):
    if x == '0/0':
        return '0'
    elif x == '0/1':
        return '1'
    else:
        return '2'

with open(argv[1], 'rb') as vcf:
    for line in vcf:
        line = bytes.decode(line).strip('\n')
        if line[0] == '#':
            continue
        else:
            lspl = line.split('\t')
            scafbp = ':'.join(lspl[0:2])
            ref, alt = lspl[3:5]
            gts = pd.Series(lspl[9:len(lspl)])
            gts = gts.map(lambda x:x.split(':', 1)[0])
            newline = str('.:' + scafbp + ' ' + ref + ' ' + alt + ' ' + ' '.join(list(gts.map(get012))))
            print(newline)

