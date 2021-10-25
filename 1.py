#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import re

for record in SeqIO.parse("spiderman.fasta", 'fasta'):
    seq = record.seq

rna = seq
rna_r = rna.reverse_complement()
rna = str(rna).upper()
rna_r = str(rna_r).upper()

rna1 = rna[1:]
rna1_r = str(Seq(rna1).reverse_complement())

rna2 = rna1[1:]
rna2_r = str(Seq(rna2).reverse_complement())

# split 6 translation frames into codons
codon1 = [rna[i:i+3] for i in range(0, len(rna), 3)]
codon2 = [rna_r[i:i+3] for i in range(0, len(rna_r), 3)]
codon3 = [rna1[i:i+3] for i in range(0, len(rna1), 3)]
codon4 = [rna1_r[i:i+3] for i in range(0, len(rna1_r), 3)]
codon5 = [rna2[i:i+3] for i in range(0, len(rna2), 3)]
codon6 = [rna2_r[i:i+3] for i in range(0, len(rna2_r), 3)]

# extract orfs
def all_orfs(codon):
    orfs = list()
    remaining = codon
    x = 1
    while x == 1:
        if 'ATG' in remaining:
            start = remaining.index('ATG')
            remaining = remaining[start:]
            end1, end2, end3 = 9999, 9999, 9999
            if "TAA" in remaining:
                end1 = remaining.index('TAA')
            if 'TAG' in remaining:
                end2 = remaining.index('TAG')
            if 'TGA' in remaining:
                end3 = remaining.index('TGA')
            end = min(end1, end2, end3)
            orfs.append(remaining[:end+1])
            remaining = remaining[end:]
        else:
            x = 0
    return orfs

# 6 translation frames' orfs
orf1 = all_orfs(codon1)
orf2 = all_orfs(codon2)
orf3 = all_orfs(codon3)
orf4 = all_orfs(codon4)
orf5 = all_orfs(codon5)
orf6 = all_orfs(codon6)

# filter orfs and translate using table2
def newOrf(orf):
    newOrf = []
    for i in orf:
        if len(i) > 20:
            seq = Seq(''.join(i)).translate(table=2)
            newOrf.append(seq)
    return newOrf

n1 = newOrf(orf1)
n2 = newOrf(orf2)
n3 = newOrf(orf3)
n4 = newOrf(orf4)
n5 = newOrf(orf5)
n6 = newOrf(orf6)

# print all AA's.
allorfs = n1 + n2 + n3 + n4 + n5
for i in allorfs:
    print(str(i))



