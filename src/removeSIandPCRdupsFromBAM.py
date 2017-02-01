#!/usr/bin/env python
"""

Date : March 18, 2016

Author : Heather Landry

Remove reads in bam file resulting from PCR duplicates. This script records all barcodes and coordinates
at a specific position. For every bam line, if the barcode and coordinate has not been seen previously, it
will print; if the barcode and position has been seen previously, it will not print to a new file. This
script will also remove reads that could be generated from splice intermediates that map to either the
5'SS or 3'SS of an intron.

use : python removeSIandPCRdupsFromBAM.py   iBAM (input BAM file with only unique alignments) [1]
                                            iSI (input file containing splicing intermediates (SI) positions with the folowing
                                                    nomenclature : chrX_strand_Start (where Start position are 1 based) [2]
                                            oBAM1 (output BAM file containing no duplicated reads) [3]
                                            oBAM2 (output BAM file containing no splicing intermediates and no duplicated reads) [4]
                                            
"""
import sys, pysam, os, numpy, re

iBAM = pysam.Samfile(sys.argv[1], 'rb')
iSI = set(open(sys.argv[2], 'r').readlines())
oBAM1 = pysam.Samfile(sys.argv[3], 'wb', template=iBAM)
oBAM2 = pysam.Samfile(sys.argv[4], 'wb', template=iBAM)

MB = set()

# read through starting bam file
for read in iBAM:
    mb = read.qname.split('_MolecularBarcode:')[1]
    chrom = iBAM.getrname(read.tid)
    
    # selecting the 3' position for pos strand 
    if read.is_reverse:
        start = read.aend
        std='pos'
        coord = str(chrom)+"_"+str(std)+"_"+str(start)+"\n"
        
    # selecting the 3' position for neg strand 
    if not read.is_reverse:
        start = read.pos
        std='neg'
        coord = str(chrom)+"_"+str(std)+"_"+str(start+1)+"\n"
    
    key = str(chrom)+"_"+str(start)+"_"+str(std)+"_"+str(mb)
    
    # output 1 read per molecular barcode
    if key not in MB:
        MB.add(key)
        oBAM1.write(read)
        # output 1 read per molecular barcode and
        # no splicing intermediates
        if coord not in iSI:
            oBAM2.write(read)
            
iBAM.close()
oBAM1.close()
oBAM2.close()

