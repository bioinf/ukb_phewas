# -*- coding: utf-8 -*-
"""
Created on Wed May  9 11:02:13 2018

@author: anton
"""
import csv
from collections import defaultdict

clust_dict = defaultdict(list)

# Название аутпута
csvwfile = open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')


# То, что получилось после кластеризации
with open('C:\\Downloads\\herit_for_clustering_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        clust_dict[row[2]].append(row[1])

        
with open('C:\\Downloads\\herit_Manhattan_09_no_0.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        row_to_write = row
        clust_list = []
        for phen in row[3].split('; '):
            for key, value in clust_dict.items():
                if phen in value:
                    clust_list.append(key)
        row_to_write[3] = '; '.join(map(str, list(set(clust_list))))  
        row_to_write[4] = len(set(clust_list))
        my_writer.writerow(row_to_write)
