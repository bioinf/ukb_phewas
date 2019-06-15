# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 23:21:00 2018

@author: anton
"""

import csv
from collections import defaultdict

clust_dict = defaultdict(list)
SNP_dict = defaultdict(list)
cluster_set = set()


with open('C:\\Downloads\\for_vcf.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        SNP_dict[row[2]].append(row[0])
        SNP_dict[row[2]].append(row[1])

with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        clust_dict[row[2]].append(row[0])
        clust_dict[row[2]].append(row[1])
        clust_dict[row[2]].append(row[4])

with open ('C:\\Downloads\\Common_SNP_for_Clusters_union_n2.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
          cluster_set.add(row[0])  
    
cluster_list = list(cluster_set)  
for cluster in cluster_list:
    csvwfile = open('C:\\Downloads\\union_n2\\Bed\\Common_SNP_jaccard_union_n2_{}.bed'.format(cluster), 'w',newline='')
    my_writer = csv.writer(csvwfile, delimiter='\t')
    csvwfile2 = open('C:\\Downloads\\union_n2\\VCF\\Common_SNP_jaccard_union_n2_{}.vcf'.format(cluster), 'w',newline='')
    my_writer2 = csv.writer(csvwfile2, delimiter='\t')  
    with open ('C:\\Downloads\\Common_SNP_for_Clusters_union_n2.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            if int(row[0])==int(cluster):
                row_to_write = [None]*4
                row_to_write[0] = "chr"+str(clust_dict[row[1]][0])
                row_to_write[1] = int(clust_dict[row[1]][1])-1000
                row_to_write[2] = int(clust_dict[row[1]][1])+1000
                row_to_write[3] = len(row[2].split('; '))
                my_writer.writerow(row_to_write)
                
                
                row_to_write = [None]*5
                row_to_write[0] = clust_dict[row[1]][0]
                row_to_write[1] = int(clust_dict[row[1]][1])
                row_to_write[2] = row[1]
                row_to_write[3] = SNP_dict[row[1]][0]
                row_to_write[4] = SNP_dict[row[1]][1]
                for i in range(4):
                    row_to_write.append('.')
                row_to_write[8] = len(row[2].split('; '))
                my_writer2.writerow(row_to_write)
                

                
                