# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 21:22:55 2018

@author: anton
"""

import csv
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy as sp
from scipy import interpolate
from scipy.stats import hypergeom
from statistics import median
Graph_dict = defaultdict(list)



SNP_names_list = dict()
with open('C:\\Downloads\\1249.assoc.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        SNP_names_list[row[1]] = "_".join(row[0].split(":"))+'_b37'

refset=set()
with open('C:\\Downloads\\GSEA-Fisher-master\\backgroundset\\all_human_genes.txt',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if len(row)==1:
            refset.add(row[0]) 
with open('C:\\Downloads\\GSEA-Fisher-master\\genesets\\sets\\all_sets.gmt',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        elems = (list)
        for el in range(2, len(row)):
            refset.add(row[el])       
comp_names = list([gene for gene in refset if "-" in gene])   
coord_gene_dict = defaultdict(list)
chr_snp_dict = defaultdict(dict)
coord_chr_snp_dict = defaultdict(dict)
print('creating additional dicts')
with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        chr_snp_dict[row[0]][row[2]]=row[1]
        coord_chr_snp_dict[row[0]][row[1]]=row[2]
        coord_gene_dict[row[2]].append(row[0]+"_"+row[1])
        coord_gene_dict[row[2]].append(str(0))

print('adding genes to dict')
with open('C:\\Downloads\\all_annotated.vcf','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if "#" not in row[0]: 
            if row[7].split("|")[3] in refset:
                coord_gene_dict[row[2]].append(row[7].split("|")[3])
                coord_gene_dict[row[2]][1] = str(int(coord_gene_dict[row[2]][1])+1)
            elif row[7].split("|")[3] not in refset and "-" not in row[7].split("|")[3]:      
                coord_gene_dict[row[2]].append(row[7].split("|")[3])
                coord_gene_dict[row[2]][1] = " not_set"
            else:
                testing_str = row[7].split("|")[3]
                genes = list([gene for gene in comp_names if gene in testing_str])
                if len(genes) > 0:
                    for gene in genes:
                        coord_gene_dict[row[2]].append(gene)
                    coord_gene_dict[row[2]][1] = str(int(coord_gene_dict[row[2]][1])+len(genes)) 
                else:
                    for el in testing_str.split("-"):
                        if el in refset:
                            coord_gene_dict[row[2]].append(el)
                            coord_gene_dict[row[2]][1] = str(int(coord_gene_dict[row[2]][1])+1)

old_SNP=list()
with open('C:\\Downloads\\herit_for_graph_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        print(row)
        old_SNP.append(row[0])
        for clust in row[1].split('; '):
            Graph_dict[clust].append(row[0])
            
clumps = defaultdict(list) 
snp_list = list()
tested_snps = list()
with open('C:\\Downloads\\pleio-clump_clumped.clumped ',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    i=0
    for row in my_reader:
        if len(row)>0 : 
            print(row)
            if row[0].split()[11]!="NONE":
                #clumps["clump"+str(i)].append(row[0].split()[2])
                #clumps["clump"+str(i)].extend([snp[0: len(snp)-3] for snp in row[0].split()[11].split(',')])
                add_list = []
                add_list.append(row[0].split()[2])
                add_list.extend([snp[0: len(snp)-3] for snp in row[0].split()[11].split(',')])
                clumps[row[0].split()[2]].append(add_list)
                tested_snps.append(row[0].split()[2])
                tested_snps.extend([snp[0: len(snp)-3] for snp in row[0].split()[11].split(',')])
                i+=1
            if row[0].split()[11]=="NONE":
                snp_list.append(row[0].split()[2])  

    

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
for clump in clumps.keys():
    intervals=[]
    try:
        for snp in clumps[clump][0]:
            intervals.append(int(coord_gene_dict[snp][0].split('_')[1]))
        clumps[clump].insert(0,coord_gene_dict[snp][0].split('_')[0])
        clumps[clump].insert(2,max(intervals)-min(intervals))
        clumps[clump].insert(3,min(intervals))
        clumps[clump].insert(4,max(intervals))
    except:
        pass  

to_add = list(set(snp_list+list(set(old_SNP).difference(set(tested_snps+snp_list)))))
not_passed=list()
len(not_passed)
for snp in to_add:
        check_flag = 0
        for clump in clumps.keys():
            if clumps[clump][0]==coord_gene_dict[snp][0].split('_')[0]:
                a = [clumps[clump][2],clumps[clump][3]]
                if int(a[0]) <= int(coord_gene_dict[snp][0].split('_')[1]) <= int(a[1]):
                    clumps[clump].append(snp)
                    check_flag=1 
                    print(snp)         
            break
        if check_flag ==0:
            not_passed.append((snp,coord_gene_dict[snp][0].split('_')[0]))        

clump_med = median(list(clumps[clump][2] for index in clumps.keys() if len(clumps[clump])>=2))  #361285 
passing_dict = defaultdict(list) 
for snp in not_passed:
    passing_dict[snp[0]].append(coord_gene_dict[snp[0]][0].split('_')[0])
    passing_dict[snp[0]].append([snp[0]])
    passing_dict[snp[0]].append(clump_med)
    passing_dict[snp[0]].append(int(coord_gene_dict[snp[0]][0].split('_')[1])-clump_med)
    passing_dict[snp[0]].append(int(coord_gene_dict[snp[0]][0].split('_')[1])+clump_med)

def clump_reducing(clumpig_universe):
    tot_set=set()
    ret_universe=defaultdict(list) 
    to_reduce=defaultdict(list)
    n=0
    for index1 in clumpig_universe:
        n+=1
        ind_check=0
        for index2 in clumpig_universe:
            if index1!=index2:
                if clumpig_universe[index1][0]==clumpig_universe[index2][0]:
                    a = [int(clumpig_universe[index1][3]),int(clumpig_universe[index1][4])]
                    b = [int(clumpig_universe[index2][3]),int(clumpig_universe[index2][4])]
                    if int(getOverlap(a, b)) != 0:
                        if max([int(getOverlap(a, b))/int(clumpig_universe[index1][2]),int(getOverlap(a, b))/int(clumpig_universe[index2][2])])>0.9:
                            ind_check=1
        if n%10==0:
            print('first', n)
        if ind_check==0:
            tot_set.update(clumpig_universe[index1][1])
            ret_universe[index1]=clumpig_universe[index1]
        else:
            to_reduce[index1]=clumpig_universe[index1]
    i=0
    passed = set()
    for index1 in list(to_reduce.keys()):
        i+=1
        if index1 not in passed:
            snp_set=set()
            new_clump=list()
            new_clump.append(index1)
            snp_set.update(to_reduce[index1][1])
            intervals=[to_reduce[index1][3],to_reduce[index1][4]]
            for index2 in list(to_reduce.keys()):
                if index2 not in passed:
                    if index1!=index2:
                        if to_reduce[index1][0]==to_reduce[index2][0]:
                            a = [int(to_reduce[index1][3]),int(to_reduce[index1][4])]
                            b = [int(to_reduce[index2][3]),int(to_reduce[index2][4])]
                            if int(getOverlap(a, b)) !=0:
                                if max([int(getOverlap(a, b))/int(to_reduce[index1][2]),int(getOverlap(a, b))/int(to_reduce[index2][2])])>=0.9:
                                    intervals[0]=min([int(x) for x in [intervals[0],a[0],b[0]]])
                                    intervals[1]=max([int(x) for x in [intervals[1],a[1],b[1]]])
                                    snp_set.update(to_reduce[index1][1])
                                    snp_set.update(to_reduce[index2][1]) 
                                    passed.add(index2)
            
            ret_universe[index1]=([to_reduce[index1][0],list(snp_set),int(intervals[1])-int(intervals[0]),int(intervals[0]),int(intervals[1])])
        passed.add(index1)
        if i% 10==0:
            print('second',i)
    return(ret_universe)
red = clump_reducing(passing_dict) 
 
for new_clump in red:
    clumps[new_clump] = red[new_clump]
    #clumps[new_clump].append(red[new_clump][0])
    #clumps[new_clump].append(red[new_clump][2])
    #clumps[new_clump].append(red[new_clump][3])
    #clumps[new_clump].append(red[new_clump][4])
    #clumps[new_clump].extend(red[new_clump][1])
    #i+=1
clumps = clump_reducing(clumps)
masked =defaultdict(list)
i=0
for clump in clumps:
    masked['clump'+str(i)].append(clumps[clump][0])
    masked['clump'+str(i)].append(clumps[clump][2])
    masked['clump'+str(i)].append(clumps[clump][3])
    masked['clump'+str(i)].append(clumps[clump][4])
    masked['clump'+str(i)].extend(clumps[clump][1])
    i+=1
clumps =masked

with open('C:\\Downloads\\clumps_and_genes.tsv','w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t')
    for clump in clumps:
        genes=set()
        for snp in clumps[clump][4:]:
            if " not_set" not in coord_gene_dict[snp][1]:
                    indexes = int(coord_gene_dict[snp][1])
                    for i in range(indexes):
                        genes.add(coord_gene_dict[snp][i+2])
        row=[clump, "; ".join(genes)]
        my_writer.writerow(row)



snp_to_clumps=dict()
for clump in clumps:
    for snp in clumps[clump][4:]:
        snp_to_clumps[snp]=clump
       
clump_graph_dict = defaultdict(list)
for key in Graph_dict.keys():
        clumpset=set()
        for snp in Graph_dict[key]:
            if snp in snp_to_clumps:
                clumpset.add(snp_to_clumps[snp])
        clump_graph_dict[key]=list(clumpset)

with open ('C:\\Downloads\\herit_clump_graph_names_jaccard.tsv', 'w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t')              
    row = [None]
    for key in Graph_dict.keys():
        row[0] = str(key)
        print(row)
        my_writer.writerow(row)

array_for_graph = np.zeros((len(clump_graph_dict), len(clump_graph_dict)))
n,m = 0,0
for clust1 in clump_graph_dict:
    for clust2 in clump_graph_dict:
        if n>len(clump_graph_dict)-1:
            n = 0
        array_for_graph[m][n] = int(len(set(clump_graph_dict[clust1]).intersection((clump_graph_dict[clust2]))))
        n += 1
    m += 1    

np.savetxt('C:\\Downloads\\herit_clump_array_for_graph_jaccard.out',array_for_graph,fmt='%10.0f', delimiter='\t') 


genesetdb= defaultdict(list)
print('creating gene set db')
with open('C:\\Downloads\\GSEA-Fisher-master\\genesets\\sets\\all_sets.gmt',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
#        print(row[0])
        elems = (list)
        for el in range(2, len(row)):
            genesetdb[row[0]].append(row[el])
            
clumpsetbd = defaultdict(list)
print('reading clump set db')
with open('C:\\Downloads\\clumpsetdb_merged_no_HLA.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        clumpsetbd[row[0]].extend(row[1].split(","))  
        
def estimate(pv, m=None, verbose=False, lowmem=False, pi0=None):
    """
    Estimates q-values from p-values
    Args
    =====
    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1
    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = sp.zeros_like(qv)
        qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv

print("reading universe from file")
clumpig_universe = defaultdict(list) 
with open('C:\\Downloads\\clumpverse_merged.tsv', 'r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        clumpig_universe[row[0]].append(row[1])
        clumpig_universe[row[0]].append(row[2].split(','))
        clumpig_universe[row[0]].append(row[3])
        clumpig_universe[row[0]].append(row[4])
        clumpig_universe[row[0]].append(row[5])


filtered_set_db=defaultdict(list)     
for cl_set in clumpsetbd:
    n1=len(clumpsetbd[cl_set])
    n2=0
    n3=0
    passed_genes=set()
    for index in clumpsetbd[cl_set]:
        index_set = set()
        for snp in clumpig_universe[index][1]:
            if " not_set" not in coord_gene_dict[snp][1]:
                indexes = int(coord_gene_dict[snp][1])
                for i in range(indexes):
                    if coord_gene_dict[snp][i+2] in genesetdb[cl_set]:
                            index_set.add(coord_gene_dict[snp][i+2])
        if len(index_set)!=0 and len(passed_genes.union(index_set))>len(passed_genes):
            filtered_set_db[cl_set].append(index)
            passed_genes = passed_genes.union(index_set)   

def GSEA(refset, quiery, clumpsetbd, genesetdb, coord_gene_dict,gen='no'):
    pval_arrays_new=list()
    N=len(refset)
    n=len(set(quiery.keys()))
    genes = defaultdict(list)
    passed_genes=set()
    nums =list()
    n_adj=0
    for index in quiery.keys():
        index_set = set()
        for snp in quiery[index]:
            if " not_set" not in str(coord_gene_dict[snp][1]):
                indexes = int(coord_gene_dict[snp][1])
                for i in range(indexes):
                    index_set.add(coord_gene_dict[snp][i+2])
                    genes[index].append(coord_gene_dict[snp][i+2])  
        if len(index_set)!=0 and len(passed_genes.union(index_set))>len(passed_genes):
            n_adj+=1
            passed_genes = passed_genes.union(index_set)
    print(n,n_adj)
    n = n_adj                           
#    for gene_set in clumpsetbd:
#        K=len(set(clumpsetbd[gene_set]))
#        count = set()
#        for index in  quiery.keys(): 
#            for snp in quiery[index]:
#                if " not_set" not in coord_gene_dict[snp][1]:
#                    indexes = int(coord_gene_dict[snp][1])
#                    for i in range(indexes):
#                        if coord_gene_dict[snp][i+2] in genesetdb[gene_set]:
#                            count.add(index)
#        k=len(count)
#        nums.append((N,n,K,k))
#        pval_arrays_new.append(hypergeom.sf(k-1, N-k,n, K))
        
    for gene_set in clumpsetbd:
        K=len(set(clumpsetbd[gene_set]))
        count = set()
        igenes = set()
        for index in  genes.keys():
            if len(set(genesetdb[gene_set]).intersection(set(genes[index])))!=0:
                count.add(index)
                igenes.update(set(genesetdb[gene_set]).intersection(set(genes[index])))
        k=len(count)
        if gen =='yes':
            k=min(len(set(genesetdb[gene_set]).intersection(igenes)),n_adj,K,k)
#        n=n-(len(count)-k)
        nums.append((N,n,K,k,";".join(list(set(genesetdb[gene_set]).intersection(igenes)))))
#        print(N,n,K,k)
        pval_arrays_new.append(hypergeom.sf(k-1, N-k,n, K))
    qval_arrays = list(estimate(np.asarray(pval_arrays_new)))
    results=defaultdict(list)
    i=0
    for geneset in clumpsetbd:
        if qval_arrays[i]<=0.05:
            results[geneset].append(pval_arrays_new[i])
            results[geneset].append(qval_arrays[i])  
            results[geneset].append(nums[i]) 
        i+=1        
    if len(results) != 0:
        df = pd.DataFrame.from_dict(results, orient='index')
        df.reset_index(inplace=True)
        df = df.sort_values(by=[1])
        return(df)
    else:
        
        return(list())
        
def HLA_GSEA(refset, quiery, genesetdb):
    pval_arrays_new=list()
    N=len(refset)
    n=len(set(quiery))
    nums =list()
    Hla_k=0
    Hist_k=0
    for index in quiery.keys():
            if int(SNP_names_list[quiery[index][4]].split("_")[0])==6 and int(SNP_names_list[quiery[index][4]].split("_")[1]) <= 33775446 and int(SNP_names_list[quiery[index][4]].split("_")[1])>=28866528:
                Hla_k+=1
            if int(SNP_names_list[quiery[index][4]].split("_")[0])==6 and int(SNP_names_list[quiery[index][4]].split("_")[1]) <= 28000446 and int(SNP_names_list[quiery[index][4]].split("_")[1])>=25000528:
                Hist_k+=1
    for gene_set in genesetdb:
        if gene_set=='HLA_locus':
            k=min(len(set(genesetdb[gene_set])), Hla_k)
        else:
            k=min(len(set(genesetdb[gene_set])), Hist_k)
        K=len(set(genesetdb[gene_set]))
        nums.append((N,n,K,k))
        pval_arrays_new.append(hypergeom.sf(k-1, N-k,n, K))
    print(nums)
    print(pval_arrays_new)
    try:
        qval_arrays = list(estimate(np.asarray(pval_arrays_new)))
    except:
        return(list())
    results=defaultdict(list)
    i=0
    for geneset in genesetdb:
        if qval_arrays[i]<=0.05:
            results[geneset].append(pval_arrays_new[i])
            results[geneset].append(qval_arrays[i])
            results[geneset].append(nums[i])  
        i+=1        
    if len(results) != 0:
        df = pd.DataFrame.from_dict(results, orient='index')
        df.reset_index(inplace=True)
        df = df.sort_values(by=[1])
#        print(qval_arrays)
#        print(df)
        return(df)
    else:
#        print(qval_arrays)
        return(list())

snp_refset = set()
with open('C:\\Downloads\\herit_Var_MAF_jaccard_all_SNP.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        snp_refset.add(row[0])


for_enrich=defaultdict(list)
len(for_enrich)
for clump in clumps:
    flag = 0
    for snp in clumps[clump][4:]:
        if int(clumps[clump][0])==6 and  (int(SNP_names_list[snp].split("_")[1]) <= 33775446 and int(SNP_names_list[snp].split("_")[1])>=28866528 or int(SNP_names_list[snp].split("_")[1]) <= 28700446 and int(SNP_names_list[snp].split("_")[1]) >= 25000528):
            flag=1
    if flag != 1:
        for_enrich[clump]=clumps[clump][4:]
 
    


df1 = GSEA(clumpig_universe,for_enrich,filtered_set_db, genesetdb, coord_gene_dict,'yes')
df1.to_csv('C:\\Downloads\\pleio_enrich_no_HLA_2.tsv',sep='\t', index=False)
 
df2 = HLA_GSEA(snp_refset, clumps, hlasetdb)
df2.to_csv('C:\\Downloads\\pleio_enrich_HIST.tsv',sep='\t', index=False)
         

clique_SNP_dict=list()
with open('C:\\Downloads\\herit_for_graph_jaccard.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if row[0] in snp_to_clumps:
            row[0]=snp_to_clumps[row[0]]
        print(row)
        clique_SNP_dict.append(row)

with open('C:\\Downloads\\herit_for_graph_clumped.tsv','w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t')
    for row in clique_SNP_dict:
        my_writer.writerow(row)


sets_dict = defaultdict(list)
with open ('C:\\Downloads\\herit_for_graph_clumped.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            if "_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])]) not in sets_dict:
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])])].append(list())
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])])].append(0)
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])])].append(list())
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])])].append(0)
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])])].append(0)
                sets_dict["_".join([str(y) for y in sorted([int(x) for x in row[1].split('; ')])])].append(list())


#counting cluster sets 
for snp_set in sets_dict:
    with open ('C:\\Downloads\\herit_for_graph_clumped.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            if len(set(row[1].split('; ')).intersection(set(snp_set.split('_')))) >= len(snp_set.split('_')):
                sets_dict[snp_set][0].append(row[0])

for snp_set in sets_dict:
    with open ('C:\\Downloads\\herit_for_graph_jaccard.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            if len(set(row[1].split('; ')).intersection(set(snp_set.split('_')))) >= len(snp_set.split('_')):
                sets_dict[snp_set][2].append(row[0])

hlasetdb = defaultdict(list)

for index in clumpig_universe.keys():
   bp =  int(SNP_names_list[index].split("_")[1])
   chrom = int(SNP_names_list[index].split("_")[0])
   if chrom==6 and  bp <= 33775446 and bp>=28866528:
       hlasetdb['HLA_locus'].append(index)
   if chrom==6  and bp <= 28000446 and bp >= 25000528:
      hlasetdb['HIST1_complex'].append(index)

for snp_set in sets_dict:
    name = snp_set
    quiery = defaultdict(list)
    clump_set = list(set(sets_dict[snp_set][0]))
    for clump in clump_set:
        snp_list = list()
        for snp in sets_dict[snp_set][2]:
            if snp_to_clumps[snp] == clump:
                snp_list.append(snp)
        quiery[clump]=snp_list
    for_enrich = defaultdict(list)
    for clump in quiery:
        flag = 0
        for snp in clumps[clump][4:]:
            if int(SNP_names_list[snp].split("_")[0])==6 and  (int(SNP_names_list[snp].split("_")[1]) <= 34775446 and int(SNP_names_list[snp].split("_")[1])>=28700446 or int(SNP_names_list[snp].split("_")[1]) <= 28700445 and int(SNP_names_list[snp].split("_")[1]) >= 25000528):
                flag=1
        if flag != 1:
            for_enrich[clump]=quiery[clump]
    print(name,quiery)
    df1 = HLA_GSEA(snp_refset, quiery, hlasetdb)
    if len(df1)!=0:
        df1.to_csv('C:\\Downloads\\Cliques_enrich\\HIST_HLA\\{}_enrich.tsv'.format(name),sep='\t', index=False)


snp_phens_dict = defaultdict(list)
with open ('C:\\Downloads\\herit_Manhattan_09_no_0.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            snp_phens_dict[row[2]].extend(row[3].split('; '))


phen_codes_dict = dict()
with open ('C:\\Downloads\\Codes_ann.csv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter=';')
        for row in my_reader1:
            phen_codes_dict[row[0]]= row[1]


cluster_phen_dict =dict()
with open ('C:\\Downloads\Cluster_diagram_jaccard.tsv', 'r',newline='') as csvfile1:
        my_reader1 = csv.reader(csvfile1, delimiter='\t')
        for row in my_reader1:
            cluster_phen_dict[row[1]]= row[0]  

              
for snp_set in sets_dict:
    sets_dict[snp_set][0]=list(set(sets_dict[snp_set][0]))
    sets_dict[snp_set][1]=len(list(set(sets_dict[snp_set][0])))
    sets_dict[snp_set][3]=len(sets_dict[snp_set][2])
    index_set=set()
    for snp in sets_dict[snp_set][2]:
        if " not_set" not in str(coord_gene_dict[snp][1]):
                indexes = int(coord_gene_dict[snp][1])
                for i in range(indexes):
                    index_set.add(coord_gene_dict[snp][i+2])
    phen_set=set()
    for snp in sets_dict[snp_set][2]:
        phen_set.update(snp_phens_dict[snp])
    names=set()
    for phen in phen_set:
        names.add('('+cluster_phen_dict[phen]+')' +phen_codes_dict[phen])
    names =list(names)
    names.sort    
    sets_dict[snp_set][2]=list(index_set)
    sets_dict[snp_set][5]=names
    sets_dict[snp_set][4]=len(index_set)  
       
   


with open('C:\\Downloads\\clumped_cliques_analysis_phens.tsv','w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter='\t')
    zer_row=["Clique", "Clump_list", "Num_of_clumps", "Gene_list","Num_of_genes","Num_of_SNPs","Num_of_nodes"]
    my_writer.writerow(zer_row)
    for clique in sets_dict:
        row=[None]*8
        row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7] = clique,"; ".join(sets_dict[clique][0]),sets_dict[clique][1],"; ".join(sets_dict[clique][2]),sets_dict[clique][4],sets_dict[clique][3],len(clique.split("_")),"; ".join(sets_dict[clique][5])
        print(row)
        my_writer.writerow(row)
        

