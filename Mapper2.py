#!/usr/bin/python
import os
import sys
import glob

def readframe(aa, pos):
  wtaa  = aa[0]
  mutaa = aa[1]
  pos   = (pos-7)/3+1 # change with different gene
  return wtaa+str(pos)+mutaa

def adjustmutpos(muts, offset):
  if muts == 'WT': return muts
  else:
    muts  = muts.rsplit('-')
    mlist = []
    for m in muts:
      pos = str(int(m[1:-1])+offset-1)
      mlist.append(m[0]+pos+m[-1])
    return '-'.join(mlist)

def list2string(l):
  newl = []
  for i in l:
    newl.append(str(i))
  return newl

#READ IN OFFSET FILE
infile  = open('Fasta/flu2offset','r')
offsetH = {}
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  offsetH[line[0]] = [int(line[1]),int(line[2])]
infile.close()
#READ IN BARCODE FILE
infile   = open('Fasta/BarCode_IFN','r')
barcodes = {}
pops     = []
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  barcodes[line[0]] = line[1]
  pops.append(line[1])
infile.close()

#GENOTYPE AND DEPTH HASH INITIATE
GenotypeH = {}
depthH    = {}
WTcount   = {}
for pop in pops:
  GenotypeH[pop] = {}
  depthH[pop]    = {}
  WTcount[pop]   = {}
  for amp in offsetH.keys():
    depthH[pop][amp] = 0
    WTcount[pop][amp] = 0

##############################MAIN##################################
infile = open('tmp_flu4_IFN/AllM','r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  if line[1] not in barcodes.keys(): continue
  amp    = line[0]
  offset = offsetH[amp][0]
  pop    = barcodes[line[1]]
  muts   = adjustmutpos(line[2],offset)
#  if muts == 'G689T' : print line
  depthH[pop][amp] += 1
  if muts == 'WT':
    WTcount[pop][amp] += 1
  else:
    if GenotypeH[pop].has_key(muts): 
      GenotypeH[pop][muts] += 1
    else: 
      GenotypeH[pop][muts] = 1
infile.close()

#SUMMARIZE GENOTYPE
Genotypes = []
for pop in pops:
  Genotypes.extend(GenotypeH[pop].keys())
Genotypes = list(set(Genotypes))

#for G in Genotypes:
#  if G.count('-')==0:
#    print G
#  elif G.count('-')==1:
#    print G

#OUTPUT
#WTCOUNT AND DEPTH OUTPUT
outfile = open('result_flu4_IFN_101315/depth','w')
outfile.write('pop'+"\t"+'amp'+"\t"+'WTcount'+"\t"+'Depth'+"\n")
for pop in pops:
  for amp in offsetH.keys():
    outfile.write("\t".join([pop, amp, str(WTcount[pop][amp]), str(depthH[pop][amp])])+"\n")
outfile.close()

#GENOTYPES
outfileA = open('result_flu4_IFN_101315/AGenotypes','w')
outfileS = open('result_flu4_IFN_101315/SGenotypes','w')
header   = 'Genotype'+"\t"+"\t".join(pops)
outfileA.write(header+"\n")
outfileS.write(header+"\n")
for G in Genotypes:
  counts = []
  for pop in pops:
    if GenotypeH[pop].has_key(G): 
      counts.append(GenotypeH[pop][G])
    else: 
      counts.append(0)
  out = G+"\t"+"\t".join(list2string(counts))
  outfileA.write(out+"\n")
  if '-' not in G:
    outfileS.write(out+"\n")
outfileA.close()
outfileS.close()
'''
#Double MUTS INFO FILE
outfileC = open('result/DMut','w')
header   = ['Pos1','Pos2','Genotype']
for pop in pops:
  header.extend([pop, pop+'_WT',pop+'_Dep'])
outfileC.write("\t".join(header)+"\n")

'''
#SINGLE MUTS INFO FILE
outfileC = open('result_flu4_IFN_101315/flu4','w')
header   = ['Pos','Genotype','Frame1','Frame2']
for pop in pops:
  header.extend([pop, pop+'_WT',pop+'_Dep'])
outfileC.write("\t".join(header)+"\n")

#READ IN PROTEIN LEVEL DATA
infile = open('/home/yushendu/RSun/Flu_M/Ref/Allflu_flu4','r')
Rframe = {}
for line in infile.xreadlines():
  array = line.rstrip().rsplit("\t")
  ID = ''.join(array[0:3])
  WTaa1 = array[10]
  Mutaa1 = array[11]
  WTaa2 = array[12]
  Mutaa2 = array[13]
#  ID = ID.replace('T291A','C291A') # flu4 BpueI silent mutation
#  ID = ID.replace('T291C','C291T')
#  ID = ID.replace('T291G','C291G')
  Rframe[ID] = [WTaa1+Mutaa1, WTaa2+Mutaa2]
infile.close()

# single mutation
for G in Genotypes:
  if '-' not in G and 'N' not in G:
    pos = int(G[1:-1])
    F1  = readframe(Rframe[G][0],pos)
    F2  = readframe(Rframe[G][1],8)
    out = [str(pos),G,F1,F2]
    for pop in pops:
      wtc = 0
      dep = 0
      for amp in offsetH.keys():
        if pos >= offsetH[amp][0] and pos <= offsetH[amp][1]:
          wtc += WTcount[pop][amp]
          dep += depthH[pop][amp]
      if GenotypeH[pop].has_key(G): out.append(str(GenotypeH[pop][G]))
      else: out.append(str(0))
      out.append(str(wtc))
      out.append(str(dep))
    if len(out) < 5: continue
    else:
      outfileC.write("\t".join(out)+"\n")
outfileC.close()


'''
# double mutation profile #
for G in Genotypes:
  if G.count('-')==1 and 'N' not in G :
    Mut=G.rsplit('-')
    pos1 = int(Mut[0][1:-1])
    pos2 = int(Mut[1][1:-1])
    out = [str(pos1),str(pos2),G]
    for pop in pops:
      wtc = 0
      dep = 0 
      for amp in offsetH.keys():
        if pos1 >= offsetH[amp][0] and pos1 <= offsetH[amp][1] and pos2 >= offsetH[amp][0] and pos2 <= offsetH[amp][1]:
          wtc += WTcount[pop][amp]
          dep += depthH[pop][amp]
      if GenotypeH[pop].has_key(G): out.append(str(GenotypeH[pop][G]))
      else: out.append(str(0))
      out.append(str(wtc))
      out.append(str(dep))
    outfileC.write("\t".join(out)+"\n")
outfileC.close()
'''
