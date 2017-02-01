#!/usr/bin/env python
"""

Date : March 18, 2016

Author : Heather Landry

Remove reads in bam file resulting from PCR duplicates. This script records all barcodes and coordinates
at a specific position. For every bam line, if the barcode and coordinate has not been seen previously, it
will print; if the barcode and position has been seen previously, it will not print to a new file.

use : python removePCRdupsFromBAM.py    iBAM (input BAM file with only unique alignments) [1]
                                        oBAM (output BAM file containing only non duplicated reads) [2]
                                            
"""
import sys, pysam, os, numpy, re

iBAM = pysam.Samfile(sys.argv[1], 'rb')
oBAM = pysam.Samfile(sys.argv[2], 'wb', template=iBAM)

MB = set()

# read through starting bam file
for read in iBAM:
    mb = read.qname.split('_MolecularBarcode:')[1]
    chrom = iBAM.getrname(read.tid)
    
    # selecting the 3' position for pos strand 
    if read.is_reverse:
        start = read.aend
        std='pos'
    
    # selecting the 3' position for neg strand 
    if not read.is_reverse:
        start = read.pos
        std='neg'

    key = str(chrom)+"_"+str(start)+"_"+str(std)+"_"+str(mb)
    
    # output 1 read per molecular barcode
    if key not in MB:
        MB.add(key)
        oBAM.write(read)

iBAM.close()
oBAM.close()

