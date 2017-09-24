#!/usr/bin/python

import os
import sys
from string import atof

allflu=open('result/result_flu1/Flu1IFN_Combine_All','r')
outfile=open('result/result_flu1/Flu1IFN_HC','w')
count=0
for line in allflu.xreadlines():
  array=line.rstrip().rsplit('\t')
  if 'Pos' in line: 
      outfile.write('\t'.join(array[0:4])+'\t'+'\t'.join(['DNAfreq','R1Afreq','R1Bfreq','IFNR1Afreq','IFNR1Bfreq','FitR1A','FitR1B','Fit','IFNR1A','IFNR1B','FitIFN','IFNS'])+'\n')
      continue
  else:
      DNAfreq=atof(array[10])/atof(array[11])
      if DNAfreq > 0.0005:
        R1freq=atof(array[4])/atof(array[5])
        R2freq=atof(array[7])/atof(array[8])
        IFNR1freq=atof(array[13])/atof(array[14])
        IFNR2freq=atof(array[16])/atof(array[17])
        fit1=R1freq/DNAfreq
        fit2=R2freq/DNAfreq
        fit3=IFNR1freq/DNAfreq
        fit4=IFNR2freq/DNAfreq
        fit = (atof(array[4])+atof(array[7]))/(atof(array[5])+atof(array[8]))/DNAfreq
        fitIFN = (atof(array[13])+atof(array[16]))/(atof(array[14])+atof(array[17]))/DNAfreq
        if fit>0:
          IFNS = fitIFN/fit
        else: 
          IFNS = 'NA'
        outfile.write('\t'.join(array[0:4])+'\t'+'\t'.join([str(DNAfreq),str(R1freq),str(R2freq),str(IFNR1freq),str(IFNR2freq),str(fit1),str(fit2),str(fit),str(fit3),str(fit4),str(fitIFN),str(IFNS)])+'\n')
    
allflu.close()
outfile.close()
      
  
