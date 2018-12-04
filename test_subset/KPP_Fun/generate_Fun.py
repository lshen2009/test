#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
os.chdir('/Users/lulushen/Documents/GC_speedup/CPU_final/Reduce_mechanism/debug/test/test_subset')


name="KPP_Fun/main.F90"
input_file=open(name,'r')
lines=input_file.readlines()
input_file.close()

for num in range(2,21):
    #if(num==13):
    #    continue
    new_lines=[]
    for line in lines:
        line=line.replace('_3','_'+str(num))        
        new_lines.append(line)
      
    
    newfile='KPP_Fun/main'+'_'+str(num)+'.F90'
    output_file = open(newfile, 'w')        
    for k in range(len(new_lines)):
        ap=new_lines[k]
        output_file.write(ap)
    output_file.close()
