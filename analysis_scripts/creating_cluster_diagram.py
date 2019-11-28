# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 17:47:46 2018

@author: anton
"""

import csv
from collections import Counter, defaultdict

d = defaultdict(dict)
Clust_dict = defaultdict(list)
Clust_total = Counter()
Clust_phen_ratio = Counter()


#То, что получается в результате кластеризации(cluster_decoder.tsv)
with open('C:\\Downloads\\herit_for_clustering_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        Clust_dict[row[2]].append(row[1])
  

#Аутпут
csvwfile = open('C:\\Downloads\\Cluster_diagram_jaccard.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')

#То, что сохраняется из Манхеттена кластеров
with open('C:\\Downloads\\herit_clust_diags_sort_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        for clust in row[1].split('; '):
             Clust_total[str(clust)]+=1
        for phen in row[2].split('; '):
             Clust_phen_ratio[str(phen)]+=1
             
         
for cluster in Clust_dict.keys():
    for phen in Clust_dict[cluster]:
        d[cluster][phen] = Clust_phen_ratio[phen]/Clust_total[cluster]


row = [None]*4            
for cluster in d:
    row[0] = cluster
    for phen in d[cluster]:
        row[1] = phen
        row[2] = d[cluster][phen]
        row[3] = len(d[cluster])
        my_writer.writerow(row)
