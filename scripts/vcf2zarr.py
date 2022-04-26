#! /usr/bin/python

# This script reads a vcf file and converts it into zarr format, for faster access from scikit-allel.

# Usage: ./vcf2zarr.py <file.vcf>

import os
import numexpr
os.environ["NUMEXPR_MAX_THREADS"]="272"
#if 'NUMEXPR_MAX_THREADS' in os.environ: os.environ.pop('NUMEXPR_MAX_THREADS')
#import numexpr
#print('NumExpr.nthreads =' + str(numexpr.nthreads))
import allel as al
import zarr
import numpy as np
import sys


vcfPath = sys.argv[1]

variants = al.read_vcf(vcfPath, numbers= {'GT': 2, 'ALT': 1}, fields= '*')	# with all (*) fields read

zarrPath = str(vcfPath.split('Bypop.')[0] + '.zarr')
al.vcf_to_zarr(vcfPath, zarrPath, fields='*', log=sys.stdout, overwrite=True)



