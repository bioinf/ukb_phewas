# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 21:32:53 2018

@author: anton
"""
from collections import defaultdict
import csv
import os
SNP_names_list = dict()
SNP_single_eqltl = defaultdict(dict)
SNP_single_eqltl_pairs = defaultdict(dict)
Multi_SNP_eqltl = defaultdict(dict)


path = 'C:\\Downloads\\single-tissue_eqtl'

for filename in os.listdir(path):
#    if 'egenes' in filename:
#        with open (path+"\\"+filename, 'r',newline='') as csvfile1:
#            my_reader1 = csv.reader(csvfile1, delimiter='\t')
#            next(my_reader1)
#            for row in my_reader1:
#                if row[18]!='.' and float(row[28]) <= 0.05:
#                    SNP_single_eqltl[row[18]][filename[0:len(filename)-14]]=row[28]
    if 'signifpairs' in filename:
         with open (path+"\\"+filename, 'r',newline='') as csvfile1:
            my_reader1 = csv.reader(csvfile1, delimiter='\t')
            next(my_reader1)
            for row in my_reader1:
                if row[0]!='.' and float(row[6]) <= 0.05:
                    SNP_single_eqltl_pairs[row[0]][filename[0:len(filename)-16]]=row[6]
          

with open('C:\\Downloads\\1249.assoc.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
#        print(row[1], "_".join(row[0].split(":"))+'_b37')
        SNP_names_list[row[1]] = "_".join(row[0].split(":"))+'_b37'
        



with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    with open('C:\\Downloads\\herit_eqtls_for_clusters_pairs.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for row in my_reader:
            SNP = SNP_names_list[row[2]]
            if len(SNP_single_eqltl_pairs[SNP])==0:
                my_writer.writerow(row[2:5]+['NA',"NA", '0'])
            elif len(SNP_single_eqltl_pairs[SNP])==1:
                my_writer.writerow(row[2:5]+list(SNP_single_eqltl_pairs[SNP].keys())+list(SNP_single_eqltl_pairs[SNP].values())+["1"])
            else: 
                my_writer.writerow(row[2:5]+["; ".join(list(SNP_single_eqltl_pairs[SNP].keys()))]+["; ".join(list(SNP_single_eqltl_pairs[SNP].values()))]+[len(SNP_single_eqltl_pairs[SNP])])


with open('C:\\Downloads\\herit_valuable_SNP_jaccard_for_func.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    with open('C:\\Downloads\\herit_eqtls_for_signif_pairs.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for row in my_reader:
            SNP = SNP_names_list[row[2]]
            if len(SNP_single_eqltl_pairs[SNP])==0:
                my_writer.writerow(row[2:5]+['NA',"NA", '0'])
            elif len(SNP_single_eqltl_pairs[SNP])==1:
                my_writer.writerow(row[2:5]+list(SNP_single_eqltl_pairs[SNP].keys())+list(SNP_single_eqltl_pairs[SNP].values())+["1"])
            else: 
                my_writer.writerow(row[2:5]+["; ".join(list(SNP_single_eqltl_pairs[SNP].keys()))]+["; ".join(list(SNP_single_eqltl_pairs[SNP].values()))]+[len(SNP_single_eqltl_pairs[SNP])])


clust_man_dict = defaultdict(list)
with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if row[2] not in clust_man_dict:   
            clust_man_dict[row[2]].append(row[3])
            clust_man_dict[row[2]].append(row[4])   
        else:
            clust_man_dict[row[2]][0]=clust_man_dict[row[2]][0]+'; '+row[3]
            clust_man_dict[row[2]][1]=str(int(clust_man_dict[row[2]][1])+1)
            

with open('C:\\Downloads\\1249.assoc.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    with open('C:\\Downloads\\herit_eqtls_for_all_pairs.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for row in my_reader:
            SNP = SNP_names_list[row[1]]
            if row[1] in clust_man_dict.keys():
#                print(row[1])
                row_list=[row[1]]+clust_man_dict[row[1]]
#                print(row_list)
                if len(SNP_single_eqltl_pairs[SNP])==0:
                    written=row_list+['NA',"NA", '0']
                elif len(SNP_single_eqltl_pairs[SNP])==1:
                    written=row_list+list(SNP_single_eqltl_pairs[SNP].keys())+list(SNP_single_eqltl_pairs[SNP].values())+["1"]
                else: 
                   written=row_list+["; ".join(list(SNP_single_eqltl_pairs[SNP].keys()))]+["; ".join(list(SNP_single_eqltl_pairs[SNP].values()))]+[len(SNP_single_eqltl_pairs[SNP])]
            else:
                row_list=[row[1],'0','0']
                if len(SNP_single_eqltl_pairs[SNP])==0:
                    written=row_list+['NA',"NA", '0']
                elif len(SNP_single_eqltl_pairs[SNP])==1:
                    written=row_list+list(SNP_single_eqltl_pairs[SNP].keys())+list(SNP_single_eqltl_pairs[SNP].values())+["1"]
                else: 
                    written=row_list+["; ".join(list(SNP_single_eqltl_pairs[SNP].keys()))]+["; ".join(list(SNP_single_eqltl_pairs[SNP].values()))]+[len(SNP_single_eqltl_pairs[SNP])]
#            print(written)
#            if len(written)!=6:
#               print(written)
            my_writer.writerow(written)


