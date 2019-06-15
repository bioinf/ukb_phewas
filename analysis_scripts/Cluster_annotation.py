# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 21:04:20 2018

@author: anton
"""


import csv
Clust_dict = {}

csvwfile = open('C:\\Downloads\\herit_Annotated_clusters_diagram_jaccard.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')

with open('C:\\Downloads\\Codes_ann.csv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter=';')
    for row in my_reader:
        Clust_dict[row[0]]= row[1]


with open('C:\\Downloads\\Cluster_diagram_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        row[1] = Clust_dict[row[1]]
        my_writer.writerow(row)
        
        
        
        
with open('C:\\Downloads\\herit_names_for_map.csv', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    with open('C:\\Downloads\\herit_annotated_names_for_map.csv','w',newline='') as csvfile1:
       my_writer = csv.writer(csvfile1, delimiter='\t')
       for row in my_reader:
           row[0] = Clust_dict[row[0]]
           row[0] = row[0].replace('Cancer code, self-reported: ', '')
           row[0] = row[0].replace('Non-cancer illness code, self-reported: ', '')
           row[0] = row[0].replace('Treatment/medication code:', 'Treatment')
           row[0] = row[0].replace('Diagnoses - main ICD10: ', '')
           print(row)
           my_writer.writerow(row)
        
            