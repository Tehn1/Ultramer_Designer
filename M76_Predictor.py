
c = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    '_': ['TAA', 'TAG', 'TGA'],
}


#Time system information
import datetime

now = datetime.datetime.now().hour
if now <12: now = "morning"
else: now="afternoon"

#Restriction Enzyme module from BioPython
from Bio.Seq import Seq
from Bio.Restriction import *

import itertools


#Take 300 bp up and downstream of target from Ensembl for input
#Generate all possible mutants
print("Good "+now)
WT ="CTCCAACTTTCTGTAGAAGATACCACCTCTCCAAATACCAAGCCGTGCCCACCTACTCCCACCACCCCAGAAACATCCCCTCCTCCTCCTCCTCCTCCTCCTTCATCTACTCCTTGTTCAGCTCACCTGACCCCCTCCTCCCTGTTCCCTTCCTCCCTGGAATCATCATCGGAACAGAAATTCTATAACTTTGTGATCCTCCACGCCAGGGCAGACGAACACATCGCCCTGCGGGTTCGGGAGAAGCTGGAGGCCCTTGGCGTGCCCGACGGGGCCACCTTCTGCGAGGATTTCCAGGTGCCGGGGCGCGGGGAGCTGAGCTGCCTGCAGGACGCCATAGACCACTCAGCTTTCATCATCCTACTTCTCACCTCCAACTTCGACTGTCGCCTGAGCCTGCACCAGGTGAACCAAGCCATGATGAGCAACCTCACGCGACAGGGGTCGCCAGACTGTGTCATCCCCTTCCTGCCCCTGGAGAGCTCCCCGGCCCAGCTCAGCTCCGACACGGCCAGCCTGCTCTCCGGGCTGGTGCGGCTGGACGAACACTCCCAGATCTTCGCCAGGAAGGTGGCCAACACCTTCAAGCCCCACAGGCTTCAG"
MutPos=101
MutRes="H"
MutList = []
newseq = ""
newseqs = []
insert = ""

for letter in MutRes:
    MutList.append(c[letter])

for p in list(itertools.product(*MutList)):
    for codon in p: insert += codon
    newseq = newseq + WT[:(MutPos-1)*3]
    newseq+=insert
    newseq += WT[((MutPos-1)*3)+(len(MutRes)*3):]
    newseqs.append(newseq)
    print(newseq)
    newseq=""
    insert=""

WTHits = []
MutHits = []
LabBatch = [AciI, AflII, AleI, BamHI, BbsI, BsrGI, BtsI, DpnI, EcoRI, HindIII, KpnI, MluI, NcoI, NdeI, NheI, NotI, PvuI, SbfI, SmaI, SpeI, XmaI]

#WT Restriction Analysis
WT = Seq(WT)

#Change here to switch between in-lab enzymes and common commercial enzymes
rb = RestrictionBatch(CommOnly)
##rb = RestrictionBatch(LabBatch)

WTAna = Analysis(rb, WT, linear=True)
for item in WTAna.full():
    if len(WTAna.full()[item]) == 1: WTHits.append(item)

#Analysis of mutants
#Identify enzymes which cut only once in one sequence and not at all in the other


for seq in newseqs:
    Mut = Seq(seq)
    rb = RestrictionBatch(rb)
    MutAna = Analysis(rb, Mut, linear=True)
    for item in MutAna.full():
        if len(MutAna.full()[item]) == 1: MutHits.append(item)
    for item in WTHits:
        if (item not in MutHits) and (MutAna.full()[item] == 0):
            print("*")
            print(item)
            print(WT)
            print(Mut)
            print(item.site)
    for item in MutHits:
        if item not in WTHits:
            print("*")
            print(item)
            print(WT)
            print(Mut)
            print(item.site)
    MutHits = []
