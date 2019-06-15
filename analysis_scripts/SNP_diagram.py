# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 01:38:41 2018

@author: anton
"""

import csv
from collections import Counter, defaultdict

Clust_dict = defaultdict(list)
Clust_total = Counter()
Clust_snp_percent = defaultdict(list)

with open('C:\\Downloads\\herit_for_clustering_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        Clust_dict[row[2]].append(row[1])

for cluster in Clust_dict:
    Clust_snp_percent[cluster] = [0]*len(Clust_dict[cluster])

csvwfile = open('C:\\Downloads\\SNP_diagram_jaccard.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')


with open('C:\\Downloads\\herit_clust_diags_sort_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        for clust in row[1].split('; '):
             Clust_total[str(clust)]+=1
             Clust_snp_percent[clust][len(set(row[2].split('; ')).intersection(set(Clust_dict[clust])))-1]+=1

for clust in Clust_snp_percent:
    for i in range(len(Clust_snp_percent[clust])):
        Clust_snp_percent[clust][i] = Clust_snp_percent[clust][i]/(Clust_total[clust])

row = [None]*4            
for cluster in Clust_dict:
    row[0] = cluster
    for i in range(len(Clust_snp_percent[cluster])):
        row[1] = Clust_snp_percent[cluster][i]
        row[2] = i+1
        row[3] = len(Clust_snp_percent[cluster])
        my_writer.writerow(row)

        
