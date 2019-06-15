# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 22:11:50 2018

@author: anton
"""


import csv
row_dict = {}

with open('C:\\Downloads\\herit_Manhattan_09_no_0.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        row_dict[row[2]]= row[4] 
        
csvwfile = open('C:\\Downloads\\herit_Var_MAF_jaccard_all_SNP.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')

       
with open('C:\\Downloads\\Var_MAF.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if row[0] in row_dict.keys():
            cur_line = row
            cur_line.append(row_dict[row[0]])
        else:
            cur_line = row
            cur_line.append(0)
        my_writer.writerow(cur_line)