#!/usr/bin/python
import os
import sys
import glob
import string
import operator 
from itertools import imap
from Bio import SeqIO

def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def list2string(l):
  newl = [] 
  for i in l:
    newl.append(str(i))
  return newl

def readrecord(filename,outfile):
  Rrecords = {}
  outfile = open(outfile,'w')
  for record in SeqIO.parse(filename, "fastq"):
    ID   = record.id
    seq  = record.seq
    qual = list2string(record.letter_annotations["phred_quality"])
    assert(not Rrecords.has_key(ID))
    outfile.write(str(ID)+"\t"+str(seq)+"\t"+'-'.join(qual)+"\n")
  outfile.close() 

def offsetcheck(R_seq,refseqs): # will it be faster if we set F and R set of reference?
  for i in range(0,len(R_seq)-240+1):
    Fseq = R_seq[i:i+240]
    Rseq = rc(Fseq)
    for ref in refseqs.keys():
      refseq  = refseqs[ref]
      Fhdist  = hamming(refseq,Fseq)
      Rhdist  = hamming(refseq,Rseq)
      if Fhdist <= 6: 
        strand = 'F'
        Roffset = i
        return [ref,Roffset,strand]
      if Rhdist <= 6: 
        strand = 'R'
        Roffset = i
        return [ref,Roffset,strand]
  return 'bad'

def MapNPair(R1tmp, R2tmp,mfile,refseqs):
  R1file = open(R1tmp,'r')
  R2file = open(R2tmp,'r')
  mfile  = open(mfile,'w')
  for line in R1file.xreadlines():
    R1record = line.rstrip().rsplit("\t")
    R1_ID    = R1record[0]
    R1_bc    = R1record[1][0:4]
    R1_seq   = R1record[1][4::]
    R1_qual  = R1record[2].rsplit('-')[4::]
    R2record = R2file.readline().rstrip().rsplit("\t") 
    R2_ID    = R2record[0]
    R2_bc    = R2record[1][0:4]
    R2_seq   = R2record[1][4::]
    R2_qual  = R2record[2].rsplit('-')[4::]
    #QUALITY CONTROL#
    assert(R1_ID == R2_ID)
    if R1_bc != R2_bc: continue
    if len(R1_seq) < 240: continue
    if len(R2_seq) < 240: continue
    #END OF QUALITY CONTROL#
    #EXTRACT OFFSET INFO
    R1_info   = offsetcheck(R1_seq,refseqs)
    R2_info   = offsetcheck(R2_seq,refseqs)
    #QUALITY CONTROL#
    if R1_info == 'bad' or R2_info == 'bad': continue
    if R1_info[0] != R2_info[0]: continue
    if R1_info[2] == R2_info[2]: continue # why do you consider this as an issue?
    #END OF QUALITY CONTROL#
    #CALL MUTATION
    WT_Amp    = R1_info[0]
    refseq    = refseqs[WT_Amp]
    R1_offset = R1_info[1]
    R1_strand = R1_info[2]
    R2_offset = R2_info[1]
    R2_strand = R2_info[2]
    R1_seq    = R1_seq[R1_offset:R1_offset+240]
    R1_qual   = R1_qual[R1_offset:R1_offset+240]
    R2_seq    = R2_seq[R2_offset:R2_offset+240]
    R2_qual   = R2_qual[R2_offset:R2_offset+240]
    if R1_strand == 'R': R1_seq = rc(R1_seq); R1_qual.reverse()
    if R2_strand == 'R': R2_seq = rc(R2_seq); R2_qual.reverse()
    Muts = []
    for n in range(0,len(refseq)):
      if R1_seq[n] != refseq[n] and R1_seq[n] == R2_seq[n] and int(R1_qual[n]) >= 30 and int(R2_qual[n]) >= 30:
        Mut = refseq[n]+str(n+1)+R1_seq[n]
        Muts.append(Mut)
    if len(Muts) == 0:
      Muts = ['WT']
    mfile.write(WT_Amp+"\t"+R1_bc[0:3]+"\t"+'-'.join(Muts)+"\n")
  R1file.close()
  R2file.close()
  mfile.close()

#READ IN REFERENCE SEQUENCE
reffile = open('/home/yushendu/RSun/240lib/Flu3/Fasta/flu3amp.fa','r')
refseqs = {}
for line in reffile.xreadlines():
  if '>' in line: 
    ID = line.rstrip().replace('>','')
  else:
    refseqs[ID] = line.rstrip()

#MAIN#
filenames = sorted(glob.glob('/home/yushendu/RSun/Seqdata/121515Miseq/*_R1_*.fastq'))
for filename in filenames:
  print filename
  fileID = filename.rsplit('_')[4].rsplit('.')[0]
  R1file = filename
  R2file = filename.replace('_R1_','_R2_')
  R1tmp  = 'tmp_flu3MC/'+fileID+'_R1'
  R2tmp  = 'tmp_flu3MC/'+fileID+'_R2'
  mfile  = 'tmp_flu3MC/'+fileID+'.m'
  readrecord(R1file,R1tmp)
  readrecord(R2file,R2tmp)
  MapNPair(R1tmp,R2tmp,mfile,refseqs)
  os.system('rm '+R1tmp+' '+R2tmp)

os.system('cat tmp_flu3MC/*.m > tmp_flu3MC/AllM')
