# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 09:38:40 2018

@author: anton
"""

import csv
import os
import collections
from collections import defaultdict
import copy


path = '/media/DISK1/phegwas/Downloads/Clumped_enrich/filtered_enrich'
pathways_dict = defaultdict(dict)
gene_sets_dict = defaultdict(list)
cluster_list = []
pathways_list = []

len(gene_sets_dict)
pathways_to_count = []
with open ('C:\\Downloads\h_all_v6_2_symbols.gmt', 'r',newline='') as csvfile2:
        my_reader3 = csv.reader(csvfile2, delimiter='\t')
        for row in my_reader3:
           pathways_to_count.append(row[0])
with open ('C:\\Downloads\c2.cp.v6.2.symbols.gmt', 'r',newline='') as csvfile3:
        my_reader4 = csv.reader(csvfile3, delimiter='\t')
        for row in my_reader4:
           pathways_to_count.append(row[0])


#
with open ('C:\\Downloads\\filtered_check_unions_6.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            pathways_list.append(row[0])

with open ('C:\\Downloads\\enrichment_stats_clusters.tsv', 'w',newline='') as csvfile1:
        my_writer1 = csv.writer(csvfile1, delimiter='\t')
        for pathway in pathways_dict:
            row=[None]*2
            row[0]=pathway
            row[1]=len(pathways_dict[pathway])
            my_writer1.writerow(row)

path1 = 'C:\\Downloads\\GSEA-Fisher-master\\data\\Unions_n2'
for filename in os.listdir(path1):
    cluster_list.append(filename[28:len(filename)-27])

colling_dict=defaultdict(dict)

for filename in os.listdir(path):
    with open(path+"/" +filename,newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')  
        next(my_reader)
        counter=0
        for row in my_reader:
            print(row[3].split(',')[3])
            name = filename.split('.')[0]
            counter+=1
            gene_sets_dict[str(row[0])].append(name)
            if str(row[0]) not in pathways_dict:
                pathways_dict[str(row[0])] = defaultdict(list)
            pathways_dict[str(row[0])][name].append(row[2])
            pathways_dict[str(row[0])][name].append(row[3])
#            pathways_dict[str(row[0])][filename.split('.')[0]].append(row[4])
            #if counter == 500:
                              #  break
 pathways_dict['NABA_MATRISOME']

count=0
for gene_set in gene_sets_dict:
        if len(gene_sets_dict[gene_set])>=2 and gene_set in pathways_to_count :
            print(gene_set, len(gene_sets_dict[gene_set]))
            count+=1

count=0
for gene_set in gene_sets_dict:
        if len(gene_sets_dict[gene_set])>=2 :
            print(gene_set, len(gene_sets_dict[gene_set]))
            count+=1




with open ('C:\\Downloads\\filtered_check_unions_6.tsv', 'w',newline='') as csvfile1:
    my_writer = csv.writer(csvfile1, delimiter='\t')
    for gene_set in gene_sets_dict:
            if len(gene_sets_dict[gene_set]) >= 6:
                row = [None]*2
                row[0] = gene_set
                row[1] = "; ".join(gene_sets_dict[gene_set])
                print(row)
                my_writer.writerow(row) 
##            if gene_set in pathways_list and len(gene_sets_dict[gene_set]) >= 2:
#                for cluster in gene_sets_dict[gene_set]:
#                    if cluster not in cluster_list:
#                        cluster_list.append(cluster)
#           

            
with open ('C:\\Downloads\\filtered_heatmap_clump_set_unions_6.tsv', 'w',newline='') as csvfile1:
    my_writer = csv.writer(csvfile1, delimiter='\t')
    first_row = copy.copy(cluster_list)
    first_row.insert(0,"category")
    my_writer.writerow(first_row)
    for pathway in pathways_dict:
        if pathway in pathways_list and len(pathways_dict[pathway]) >= 2:
            row = [None]*30
            row[0] = pathway
            for i in range(1,len(cluster_list)+1):
               if cluster_list[i-1] in pathways_dict[pathway]:
                   row[i]=pathways_dict[pathway][cluster_list[i-1]][0][0]
               else:
                   row[i]='NA'
            print(row)
            my_writer.writerow(row)


phen_dict = defaultdict(dict)
phen_path = 'C:\\Downloads\\Clumped_enrich\\Phenotypes_enrich'
path_phen_dict = defaultdict(dict)
for phen in  os.listdir(phen_path):
    name =  phen.split(".")[0].replace('_phen', '')
    if name not in phen_dict:
                phen_dict[name] = defaultdict(dict)
    with open(phen_path +"\\"+phen,newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')  
        next(my_reader)
        for row in my_reader:
            if row[0] not in path_phen_dict:
                path_phen_dict[row[0]] = defaultdict(dict)
            path_phen_dict[str(row[0])][name]=float(row[2])  
            phen_dict[name][str(row[0])]=float(row[2])   

with open ('C:\\Downloads\\enrichment_stats_total.tsv', 'w',newline='') as csvfile1:
        my_writer1 = csv.writer(csvfile1, delimiter='\t')
        row0=['pathway', 'n_phen', 'n_clust', "enrich_un_n3", "un_n3_erich"]
        my_writer1.writerow(row0)
        for pathway in path_phen_dict:
            if True: 
#            if pathway in pathways_to_count:   
               clusters=set()
               clusters_e_u=set()
               clusters_u_e=set()
               for phen in path_phen_dict[pathway].keys():
                    clusters.add(phen_to_clust_dict[phen])
                    if len(clust_dict[phen_to_clust_dict[phen]])>=3:
                        clusters_e_u.add(phen_to_clust_dict[phen])
               row=[None]*5
               row[0]=pathway
               row[1]=len(path_phen_dict[pathway])
               row[2]=len(clusters)
               row[3]=len(clusters_e_u)
               row[4]=len(pathways_dict[pathway])
               my_writer1.writerow(row)


with open ('C:\\Downloads\\enrichment_cluster_dist.tsv', 'w',newline='') as csvfile1:
        my_writer1 = csv.writer(csvfile1, delimiter='\t')
        row0=['pathway', 'n_phen', 'n_clust', "enrich_un_n3", "un_n3_erich"]
        my_writer1.writerow(cluster_list)
        for cluster in path_phen_dict:
            if True: 
#            if pathway in pathways_to_count:   
               clusters=set()
               clusters_e_u=set()
               clusters_u_e=set()
               for phen in path_phen_dict[pathway].keys():
                    clusters.add(phen_to_clust_dict[phen])
                    if len(clust_dict[phen_to_clust_dict[phen]])>=3:
                        clusters_e_u.add(phen_to_clust_dict[phen])
               row=[None]*5
               row[0]=pathway
               row[1]=len(path_phen_dict[pathway])
               row[2]=len(clusters)
               row[3]=len(clusters_e_u)
               row[4]=len(pathways_dict[pathway])
               my_writer1.writerow(row)


with open ('C:\\Downloads\\pleio_enrich_no_HLA_filtered_2.tsv', 'w',newline='') as csvfile1:
        my_writer1 = csv.writer(csvfile1, delimiter='\t')
        with open ('C:\\Downloads\\pleio_enrich_no_HLA_2.tsv', 'r',newline='') as csvfile1:
            my_reader1 = csv.reader(csvfile1, delimiter='\t')
            for row in my_reader1:
                if row[0] in  pathways_to_count or row[0]=="HSIAO_LIVER_SPECIFIC_GENES":
                    my_writer1.writerow(row)



clust_dict = defaultdict(list)
phen_to_clust_dict=dict()
with open ('C:\\Downloads\\Cluster_diagram_jaccard.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            phen_to_clust_dict[row[1]]=row[0]
            clust_dict[row[0]].append(row[1])        

unioned_n2_set = defaultdict(dict)
def second_smallest(numbers):
    m1, m2 = float('inf'), float('inf')
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2 

for cluster in cluster_list:
    if cluster  not in unioned_n2_set:
       unioned_n2_set[cluster] = defaultdict(list)
    cash=defaultdict(list)
    for phen in clust_dict[cluster]:
        for enrich in phen_dict[str(phen)]:
            cash[enrich].append(phen)
    for enrich in cash:
            pvals = []
            for phen in cash[enrich]:
                try:
                    pvals.append(float(phen_dict[str(phen)][enrich]))
                except(KeyError):
                    pass
            print(pvals)
            if len(pvals)>2:
                unioned_n2_set[cluster][enrich].append(second_smallest(pvals))
            else:
                unioned_n2_set[cluster][enrich].append(pvals[0])

pathways_dict = defaultdict(dict)
gene_sets_dict = defaultdict(list)
for cluster in unioned_n2_set: 
#    counter=0         
    sorted_by_value = sorted(unioned_n2_set[cluster].items(), key=lambda kv: kv[1])
    for pair in sorted_by_value:
        if pair[0] not in pathways_dict:
            pathways_dict[pair[0]] = defaultdict(list)
        pathways_dict[pair[0]][cluster].append(pair[1])   
        if pair[0] in pathways_to_count:
#            counter+=1
            gene_sets_dict[pair[0]].append(cluster)   
#        if counter == 20:
#            break
count=0
for gene_set in gene_sets_dict:
        if len(gene_sets_dict[gene_set]) >= 2:
            count+=1
for phen in clust_dict['153']:
    print(phen,phen_dict[str(phen)])
    
    
    
phen_dict['21001']['REACTOME_IMMUNE_SYSTEM']
phen_dict['21001']['REACTOME_IMMUNE_SYSTEM']
for cluster in clust_dict.keys():
    phen_list = list()
    if len(clust_dict[str(cluster)])==1:
        phen_pathways_dict = dict()
        phen_list.append(clust_dict[cluster][0])
        for enrich in phen_dict[clust_dict[str(cluster)][0]]:
            if enrich in pathways_to_count:
                phen_pathways_dict[enrich]=phen_dict[str(clust_dict[str(cluster)][0])][enrich]
        if len(phen_pathways_dict)>1:
            with open ('C:\\Downloads\\Phen_heatmaps\\1_{}heatmap.tsv'.format(cluster), 'w',newline='') as csvfile1:
                my_writer = csv.writer(csvfile1, delimiter='\t')
                first_row = copy.copy(phen_list)
                first_row.insert(0,"category")
                my_writer.writerow(first_row)
                counter=0
                for pathway in phen_pathways_dict:
                    row = [None]*2
                    row[0] = pathway
                    row[1] = phen_pathways_dict[pathway]
                    my_writer.writerow(row)
                    counter+=1
                    if counter>=40:
                        break                   
    elif len(clust_dict[cluster])==2:
        phen_pathways_dict = defaultdict(dict)
        phen_list.append(clust_dict[cluster][0])
        phen_list.append(clust_dict[cluster][1])
        for phen in clust_dict[str(cluster)]:
            for enrich in phen_dict[phen]:
                if enrich in pathways_to_count:
                    phen_pathways_dict[enrich][phen]=phen_dict[str(phen)][enrich]
        if len(phen_pathways_dict)>0:

            with open ('C:\\Downloads\\Phen_heatmaps\\2_{}heatmap.tsv'.format(cluster), 'w',newline='') as csvfile1:
                my_writer = csv.writer(csvfile1, delimiter='\t')
                first_row = copy.copy(phen_list)
                first_row.insert(0,"category")
                my_writer.writerow(first_row)
                counter=0
                for pathway in phen_pathways_dict:
                    row = [None]*3
                    row[0] = pathway
                    for i in range(1,3):
                        if phen_list[i-1] in phen_pathways_dict[pathway]:
                            row[i]=phen_pathways_dict[pathway][phen_list[i-1]]
                            if str(row[i])=='{}':
                                 row[i]='NA'
                        else:
                            row[i]='NA'
                    counter+=1
                    my_writer.writerow(row)
                    if counter>=40:
                        break  
    else:
        phen_pathways_dict = defaultdict(dict)
        phen_list = clust_dict[cluster]
        for phen in clust_dict[str(cluster)]:
            for enrich in phen_dict[phen]:
                if enrich in pathways_to_count:
                    phen_pathways_dict[enrich][phen]=phen_dict[str(phen)][enrich]
        print(phen_list)
        for enrich in phen_pathways_dict:
            if len(phen_pathways_dict[enrich])>=2:
                print(phen_pathways_dict[enrich])
        counter_dict=dict()
        for enrich in phen_pathways_dict:
            if len(phen_pathways_dict[enrich])>=2:
                pvals=[]
                for phen in phen_pathways_dict[enrich]:
                    pvals.append(phen_pathways_dict[enrich][phen])
                counter_dict[sum(pvals)/len(pvals)]=enrich
        for enr in collections.OrderedDict(sorted(counter_dict.items())):
            print(counter_dict[enr],enr)
                
        if len(phen_pathways_dict)>0:
            with open ('C:\\Downloads\\Phen_heatmaps\\n2\\n2_{}heatmap.tsv'.format(cluster), 'w',newline='') as csvfile1:
                my_writer = csv.writer(csvfile1, delimiter='\t')
                first_row = copy.copy(phen_list)
                first_row.insert(0,"category")
                my_writer.writerow(first_row)
                counter=0

                for pathway in collections.OrderedDict(sorted(counter_dict.items())):
                    if len(phen_pathways_dict[counter_dict[pathway]])>=2:
                        row = [None]*(len(phen_list)+1)
                        row[0] = counter_dict[pathway]
                        for i in range(1,len(phen_list)+1):
                            if phen_list[i-1] in phen_pathways_dict[counter_dict[pathway]]:
                                row[i]=phen_pathways_dict[counter_dict[pathway]][phen_list[i-1]]
                                if str(row[i])=='{}':
                                    row[i]='NA'
                            else:
                                row[i]='NA'
                        counter+=1
                        my_writer.writerow(row)
                        if counter>=50:
                            break      
filtered_heatmap_list[2]
loci_heatmap[2]
filtered_heatmap_list=list()

with open ('/media/DISK1/phegwas/Downloads/filtered_heatmap_clump_CP.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            print(row)
            filtered_heatmap_list.append(row)
clusters.index('33')


clusters = filtered_heatmap_list[0][1:]
pathways = [item[0] for item in filtered_heatmap_list[1:]]

loci_heatmap=[]
loci_heatmap.append(filtered_heatmap_list[0])

for pathway  in pathways:
    test_list=['NA']*len(filtered_heatmap_list[0])
    test_list[0]=pathway
    for cluster in pathways_dict[pathway].keys():
        if cluster!=1:
                #print(clusters.index(cluster)+1)
            test_list[clusters.index(cluster)+1]=pathways_dict[pathway][cluster][3].split(',')[3]
                #print(pathway,cluster,pathways_dict[pathway][cluster][3].split(',')[3])
    loci_heatmap.append(test_list)
 

with open ('/media/DISK1/phegwas/Downloads/filtered_heatmap_clump_CP_loci.tsv', 'w',newline='') as csvfile1:
        my_reader1 = csv.writer(csvfile1, delimiter='\t')
        for row in loci_heatmap:
            print(row)
            my_reader1.writerow(row)