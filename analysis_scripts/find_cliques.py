# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 22:53:35 2018

@author: anton
"""

from collections import defaultdict
import csv 
import networkx as nx
import statistics
import copy



graph_vertexes = []

#make dict for vcf

SNP_dict = defaultdict(list)

with open('C:\\Downloads\\for_vcf.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        SNP_dict[row[2]].append(row[0])
        SNP_dict[row[2]].append(row[1])


#create list of edges
with open('C:\\Downloads\\herit_graph_names_jaccard.tsv', 'r',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        graph_vertexes.append(row[0])

graph = defaultdict(list)

adjacency_matrix = []

#build adjacency matrix
with open("C:\\Downloads\\herit_clump_array_for_graph_jaccard.out") as csvfile:
    reader = csv.reader(csvfile, delimiter='\t') 
    for row in reader: 
        adjacency_matrix.append(list([x.replace(" ", "") for x in row]))

#convert to adjacency list
for i  in range(len(graph_vertexes)):
    for  u in range(len(adjacency_matrix[i])):
        if int(adjacency_matrix[i][u]) != 0 :
            graph[graph_vertexes[i]].append({graph_vertexes[u]: adjacency_matrix[i][u]})


#convert to nx graph
G=nx.Graph()
for node in graph:
    for neigh in graph[node]:      
        G.add_edge(node, list(neigh.keys())[0],weight=int(list(neigh.values())[0]))

G.number_of_edges()

nx.write_gml(G, 'C:\\Downloads\\clump_herit_graph')


#perform clique recursive reduction    
clique_list = list(nx.find_cliques(G))
clique_dict = defaultdict(list)
value_list = defaultdict(list)
new_clique_list = []


for clique in clique_list:
    name = "clique_"+"_".join(clique)
    value = 0
    nodes_set = set()
    print(name)
    for node1 in clique:
       for node2 in clique:
           for node3 in graph[node1]:
               if list(node3.keys())[0] == node2:
                   if frozenset([node1,list(node3.keys())[0]]) not in nodes_set:
                       nodes_set.add(frozenset([node1,list(node3.keys())[0]]))
                       value+=int(list(node3.values())[0])
#                       value_list[name][str("_".join([node1,list(node3.keys())[0]]))] = int(list(node3.values())[0])
                       value_list[name].append(int(list(node3.values())[0]))
    clique_dict[value]=clique

def Median_search(clique, graph):
    values = []
    nodes_set = set()
    for node1 in clique:
        for node2 in clique:
           for node3 in graph[node1]:
               if list(node3.keys())[0] == node2:
                   if frozenset([node1,list(node3.keys())[0]]) not in nodes_set:
                      nodes_set.add(frozenset([node1,list(node3.keys())[0]]))
                      values.append(int(list(node3.values())[0]))
    return(statistics.mean(values)/2)
    
def return_check(clique, graph):
    clique_med = Median_search(clique, graph)
    first_flag = 0
    for node1 in clique:
       inner_flag=0
       for node2 in clique:
           for node3 in graph[node1]:
               if list(node3.keys())[0] == node2:
                   if int(list(node3.values())[0]) >= clique_med:
#                       print(node1, list(node3.keys())[0], list(node3.values())[0], clique_med)
                       inner_flag=1
       if inner_flag == 0:
           first_flag = 1 
    if first_flag == 0:
        return('yes')
    else:
        return('no')

                
def Clique_search(clique, graph):
    if return_check(clique, graph) == 'yes':
        return(clique)
    else:
        cash_clique = copy.copy(clique)
        clique_med =  Median_search(cash_clique, graph)
        for node1 in cash_clique:
            flag = 0
            for node2 in cash_clique:
                for node3 in graph[node1]:
                    if list(node3.keys())[0] == node2:
#                        print("checked",node1, list(node3.keys())[0], list(node3.values())[0], clique_med)
                        if int(list(node3.values())[0]) >= clique_med:
                            flag=1
            if flag == 0:
                cash_clique.remove(node1)
        return(Clique_search(cash_clique, graph))
      
new_clique_list = set()
for clique in clique_list:
    print(clique)
    print(Clique_search(clique, graph))
    new_clique_list.add(frozenset(Clique_search(clique, graph)))

#finding all cluster sets 

sets_dict = {}
with open ('C:\\Downloads\\herit_valuable_SNP_jaccard_for_func.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
#            print("_".join([str(y) for y in sorted([int(x) for x in row[3].split('; ')])]))
            if "_".join([str(y) for y in sorted([int(x) for x in row[3].split('; ')])]) not in sets_dict:
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[3].split('; ')])])]=0
#counting cluster sets 
for snp_set in sets_dict:
    with open ('C:\\Downloads\\herit_valuable_SNP_jaccard_for_func.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            if len(set(row[3].split('; ')).intersection(set(snp_set.split('_')))) >= len(snp_set.split('_')):
                sets_dict[snp_set]+=1

#file for sets stats
with open ('C:\\Downloads\\herit_signif_SNP_sets.tsv', 'w',newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')    
        for snp_set in sets_dict:
            row = [None]*2
            row[0] = snp_set
            row[1] = sets_dict[snp_set]
            my_writer.writerow(row)


#whole graph analysis
with open ('C:\\Downloads\\clump_herit_graph_analysis.tsv', 'w',newline='') as csvfile1:
        my_writer = csv.writer(csvfile1, delimiter='\t')    
        for cluster in graph_vertexes:
            row = [None]*5
            row[0] = cluster
            row[1] = len(graph[str(cluster)])
            row[2] = sum([int(list(value.values())[0]) for value in graph[str(cluster)]])
            row[3] = row[2]/row[1]
            row[4] = "; ".join([str(list(value.keys())[0]) for value in graph[str(cluster)]])
            print(row)
            my_writer.writerow(row)


#VCFS and BED for cliques  
for clique in sets_dict:
    name = "clique_"+clique
    csvwfile = open('C:\\Downloads\\cliques\\Bed\\{}.bed'.format(name), 'w',newline='')
    my_writer = csv.writer(csvwfile, delimiter='\t')
    csvwfile2 = open('C:\\Downloads\\cliques\\VCF\\{}.vcf'.format(name), 'w',newline='')
    my_writer2 = csv.writer(csvwfile2, delimiter='\t')  
    with open ('C:\\Downloads\\herit_valuable_SNP_jaccard_for_func.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            if set(list(row[3].split('; '))).intersection(set(list(clique.split('_')))) == set(list(clique.split('_'))):
                row_to_write = [None]*4
                row_to_write[0] = "chr"+str(row[0])
                row_to_write[1] = int(row[1])-1000
                row_to_write[2] = int(row[1])+1000
                row_to_write[3] = len(row[3].split('; '))
                my_writer.writerow(row_to_write)               
                
                row_to_write = [None]*5
                row_to_write[0] = str(row[0])
                row_to_write[1] = int(row[1])
                row_to_write[2] = row[2]
                row_to_write[3] = SNP_dict[row[2]][0]
                row_to_write[4] = SNP_dict[row[2]][1]
                for i in range(4):
                    row_to_write.append('.')
                row_to_write[8] = len(row[3].split('; '))
                my_writer2.writerow(row_to_write)
                

