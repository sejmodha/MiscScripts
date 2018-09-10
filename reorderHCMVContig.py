#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:56:17 2018

@author: sejalmodha
"""
import argparse
from Bio import Seq,SeqIO
from pyfaidx import Fasta
from pymummer import nucmer

#Define all commandline arguments
parser = argparse.ArgumentParser(description='This script runs mummer using python')
parser.add_argument('-q','--query', help='Input query genome fasta file name with one sequence',required=True)
parser.add_argument('-r','--reference',help='Input reference fasta file', required=True)
parser.add_argument('-c','--coordinates',help='Output coordinates file name',required=True)
parser.add_argument('-i','--minidentity',help='Minimum percent identity for filtering results',required=True)
parser.add_argument('-o','--output',help='Output file in fasta format',required=True)
parser.add_argument('-n','--header',help='Header name for output fasta format',required=False)
args = parser.parse_args()

#function to check input file format - fasta
def checkFormat(file):
    #print(file)
    inputfile=open(file,'r')
    fasta=SeqIO.parse(inputfile,'fasta')
    return any(fasta)
    
if checkFormat(args.query) == False:
    print('Error: Please ensure that QUERY (-q) file is in fasta format')
    exit()
if checkFormat(args.reference) == False:
    print('Error: Please ensure that REFERENCE (-r) file is in fasta format')
    exit()
#check if valid parameter is provided for identity value
if args.minidentity.isdigit() ==False:
    print('Error: Identity (-i) threshold must be a numeric value')
    exit()
#Set a default header if user-defined header name is absent
if args.header == None:
    args.header='reordered_contig'
print(args.header)

#Save query file in faidx index    
contigfile=Fasta(args.query)    

#run nucmer program and filter the coords output file
runner = nucmer.Runner(args.reference,args.query,args.coordinates,min_id=args.minidentity,coords_header=False)
runner.run()

#open output files
coordsfile=open(args.coordinates)
outfile=open(args.output,'w')

#reorder the sequences based on the reference genome coordinates
reordered=''
#print(contigfile[0].name)
for line in coordsfile:
    fields=line.replace('[BEGIN]','').rstrip('\n').split('\t')[:-1]
    #print(fields)
    start=int(fields[2])
    end=int(fields[3])
    contigid=''.join((fields[-1:]))
    print(str(start)+'\t'+str(end)+'\t'+contigid)
    if start > end:
        segment=Seq.Seq(contigfile[contigid][end:start].seq).reverse_complement()
        #print(segment)
    else:
        segment=contigfile[contigid][start:end].seq
    reordered=reordered+segment
outfile.write('>'+args.header+'\n'+str(reordered))


#close output files
outfile.close()
coordsfile.close()