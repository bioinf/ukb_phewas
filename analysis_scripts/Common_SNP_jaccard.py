# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 22:18:17 2018

@author: anton
"""

import csv
from collections import Counter, defaultdict

Clust_dict = defaultdict(list)


with open('C:\\Downloads\\herit_for_clustering_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        Clust_dict[row[2]].append(row[1])
        


csvwfile = open('C:\\Downloads\\Common_SNP_for_Clusters_union_n2.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')
  
with open('C:\\Downloads\\herit_clust_diags_sort_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        for clust in row[1].split('; '):
#            if  len (Clust_dict[clust])>=3 and len(set(row[2].split('; ')).intersection(set(Clust_dict[clust])))==len(Clust_dict[clust]):
            if  len (Clust_dict[clust])>=3 and len(set(row[2].split('; ')).intersection(set(Clust_dict[clust])))>=2:
               row_to_write = [None]*5 
               row_to_write[0] = clust
               row_to_write[1] = row[0]
               row_to_write[2] = str(('; ').join(list(set(row[2].split('; ')).intersection(set(Clust_dict[clust])))))
               row_to_write[3] = row[1]
               row_to_write[4] = len(row[1].split('; '))
               my_writer.writerow(row_to_write)
               
