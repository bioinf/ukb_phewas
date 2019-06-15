# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 19:33:31 2018

@author: anton
"""

import csv
from collections import defaultdict

clust_dict = defaultdict()


csvwfile = open('C:\\Downloads\\herit_Clustered_mahattan_pval_jaccard.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')

with open('C:\\Downloads\\herit_pval_SNP.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        clust_dict[row[0]] = row[5]

with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if row[2] in clust_dict:
            row.append(clust_dict[row[2]])
            my_writer.writerow(row)