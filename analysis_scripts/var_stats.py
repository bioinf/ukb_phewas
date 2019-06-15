# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 20:11:58 2018

@author: anton
"""
from collections import Counter, defaultdict
import csv
import os
import heapq


cnt = Counter()

with open('C:\\Downloads\\variant stats\\annotated_clustered_jaccard.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        cnt [row[0]]+=1
        

with open('C:\\Downloads\\variant stats\\all_variants_stats.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for var_type in cnt:
            row = [None]*2
            row[0], row[1] =  var_type, cnt[var_type]
            my_writer.writerow(row)
            
cnt2 = Counter()

with open('C:\\Downloads\\variant stats\\pleio_all_stats.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        cnt2 [row[0]]+=1
        

with open('C:\\Downloads\\variant stats\\pleio_variants_stats.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for var_type in cnt2:
            row = [None]*2
            row[0], row[1] =  var_type, cnt2[var_type]
            my_writer.writerow(row)
            

snp_all = defaultdict(dict)

with open('C:\\Downloads\\herit_merged_09.tsv','r', newline='') as csvfile1:
        my_reader = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader:
            for i in range(len( row[2].split("; "))):
                snp_all[row[1].split('; ')[0]][row[2].split("; ")[i]]=row[9].split('; ')[i]

clusters_dict = defaultdict(list)
with open('C:\\Downloads\\Cluster_diagram_jaccard.tsv','r', newline='') as csvfile1:
        my_reader = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader:
                clusters_dict[row[0]].append(row[1])

                
path = 'C:\\Downloads\\union_n2\\VCF'
for filename in os.listdir(path):
     with open ("C:\\Downloads\\for_DEPICT\\Unions_n2"+"\\"+filename + "_depict.tsv", 'w',newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter='\t')
        phen_list = clusters_dict[filename.split("_")[5][0:len(filename.split("_")[5])-4]] 
        with open (path+"\\"+filename, 'r',newline='') as csvfile1:
            my_reader1 = csv.reader(csvfile1, delimiter='\t')
            my_writer.writerow(['SNP','Chr', 'Pos', 'P', "phen_list"])
            for row in my_reader1:
                names = []
                pvals = []
                for phen in phen_list:
                    if phen in snp_all[row[2]]:
                        pvals.append(snp_all[row[2]][phen])
                        names.append(phen)
                my_writer.writerow([row[2], row[0], row[1], min(heapq.nlargest(2, pvals)), "; ".join(names)])


with open('C:\\Downloads\\herit_valuable_SNP_jaccard.tsv','r', newline='') as csvfile1:
        my_reader = csv.reader(csvfile1, delimiter='\t')
        with open('C:\\Downloads\\herit_signif_SNP_jaccard.vcf','r', newline='') as csvfile2:
            my_reader2 = csv.reader(csvfile2, delimiter='\t')
            with open('C:\\Downloads\\herit_signif_SNP_jaccard_depict.tsv','w', newline='') as csvfile3:
                my_writer1 = csv.writer(csvfile3, delimiter='\t')
                my_writer1.writerow(['SNP','Chr', 'Pos', 'P'])
                for row1,row2 in zip(my_reader, my_reader2):
                    row=[None]*4
                    row[0], row[1], row[2], row [3]= row1[0], row2[0], row2[1], row1[2]
                    my_writer1.writerow(row)


snp_coords_dict = defaultdict(list)
with open('C:\\Downloads\\herit_Clustered_mahattan_pval.tsv','r', newline='') as csvfile1:
        my_reader = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader:
            snp_coords_dict[row[2]].append(row[0])
            snp_coords_dict[row[2]].append(row[1])


for cluster in clusters_dict:
    for phen in clusters_dict[cluster]:
        with open ("C:\\Downloads\\for_DEPICT\\clumping"+"\\"+phen+ "_phen_to_clump_depict.tsv", 'w',newline='') as csvfile:
            my_writer = csv.writer(csvfile, delimiter='\t')
            my_writer.writerow(['SNP','Chr', 'Pos', 'P'])
            for snp in snp_all:
                if phen in snp_all[snp]:
                    if len(snp_coords_dict[snp])!=0:
                        my_writer.writerow([snp,snp_coords_dict[snp][0], snp_coords_dict[snp][1], snp_all[snp][phen]])

                        
                        

cliques_dict = defaultdict(list)
path = 'C:\\Downloads\\clique_gene_sets'
for filename in os.listdir(path):
    print(filename.split('.')[0][0:len(filename.split('.')[0])-10])
    with open(path+"\\"+filename, 'r', newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        genes = []
        for row in my_reader:
            if len(row)!=0: 
                genes.append(row[0])
        if len(genes)!=0 and len(genes)<=5:
            cliques_dict[filename.split('.')[0][0:len(filename.split('.')[0])-10]].append(len(genes))
            cliques_dict[filename.split('.')[0][0:len(filename.split('.')[0])-10]].append('; '.join(genes))
            cliques_dict[filename.split('.')[0][0:len(filename.split('.')[0])-10]].append(len(filename.split('.')[0][0:len(filename.split('.')[0])-10].split('_'))-1)
  


with open('C:\\Downloads\\variant stats\\small_cliques_stats.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for clique in cliques_dict:
            my_writer.writerow([clique, cliques_dict[clique][0], cliques_dict[clique][1]])

  

with open('C:\\Downloads\\variant stats\\all_cliques_stats.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for clique in cliques_dict:
            my_writer.writerow([clique, cliques_dict[clique][0], cliques_dict[clique][1], cliques_dict[clique][2]])


path = 'C:\\Downloads\\depict_gene_sets'
for filename in os.listdir(path):   
     with open ("C:\\Downloads\\depict_reenrich"+"\\"+filename + "_extracted.tsv", 'w',newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter='\t')
        with open (path+"\\"+filename, 'r',newline='') as csvfile1: 
            my_reader1 = csv.reader(csvfile1, delimiter='\t')
            next(my_reader1)
            for row in my_reader1:
                if row[5] != '-':
                    my_writer.writerow([str(row[5])])
        

clust_man_dict = defaultdict(list)
with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if row[2] not in clust_man_dict:   
            clust_man_dict[row[2]].append(row) 
            

with open('C:\\Downloads\\1249.assoc.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    with open('C:\\Downloads\\for_all_SNP_stats.tsv','w', newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')
        for row in my_reader:
            if row[1] in clust_man_dict.keys():
                written=clust_man_dict[row[1]][0]
            else:
                written = [None]*5
                written[0:2],written[2],written[3:5]=row[0].split(":")[0:2],row[1], ['0', '0']
            my_writer.writerow(written)











        