# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 17:20:52 2018

@author: anton
"""
import csv
import numpy as np
import sys

maxInt = sys.maxsize
decrement = True

while decrement:
    decrement = False
    try:
        csv.field_size_limit(maxInt)
    except OverflowError:
        maxInt = int(maxInt/10)
        decrement = True


file_dict = {}

with open('C:\\Downloads\\herit_merged_data_for_heatmap.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        file_dict[row[1]] = row[0].split("; ")
        
def correl(key1,key2):
    return(len(set(file_dict[key1]).intersection(file_dict[key2]))/len(set(file_dict[key1]).union(set(file_dict[key2]))))
    
    
array_to_analyze = np.zeros((len(file_dict), len(file_dict)))
n,m = 0,0


for key1 in file_dict:
    for key2 in file_dict:
        if n>len(file_dict)-1:
            n = 0
        array_to_analyze[m][n] = correl(key1,key2)
        print(m,n)
        n += 1
    m += 1
    
    
np.savetxt('C:\\Downloads\\herit_Jaccard_total_array.out', array_to_analyze, delimiter='\t')    