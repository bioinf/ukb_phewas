# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 18:47:35 2018

@author: anton
"""
from collections import defaultdict
import csv
import pandas as pd
import scipy as sp
from scipy import interpolate
from scipy.stats import hypergeom
import numpy as np
import os
import copy


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


genesetdb = defaultdict(list)
refset = set()

with open('C:\\Downloads\\herit_Var_MAF_jaccard_all_SNP.tsv',newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        if row[0] not in refset:
            refset.add(row[0])

SNP_names_list = dict()
with open('C:\\Downloads\\1249.assoc.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    next(my_reader)
    for row in my_reader:
        SNP_names_list[row[1]] = "_".join(row[0].split(":"))+'_b37'


chr_snp_dict = defaultdict(dict)
with open('C:\\Downloads\\herit_Clustered_mahattan_jaccard.tsv','r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        chr_snp_dict[row[0]][row[2]]=row[1]

clumpig_universe = defaultdict(list) 
with open('C:\\Downloads\\clumpverse_merged.tsv', 'r', newline='') as csvfile:
    my_reader = csv.reader(csvfile, delimiter='\t')
    for row in my_reader:
        clumpig_universe[row[0]].append(row[1])
        clumpig_universe[row[0]].append(row[2].split(','))
        clumpig_universe[row[0]].append(row[3])
        clumpig_universe[row[0]].append(row[4])
        clumpig_universe[row[0]].append(row[5])
#def get_quiery(path):
#    quiery = defaultdict(list)
#    with open (path,newline='') as csvfile:
#        my_reader = csv.reader(csvfile, delimiter='\t')
#        next(my_reader)
#        for row in my_reader:
#            try:
#                row = list([el.replace('(1)','') for el in row[0].split()])
#                quiery[row[2]].append(row[2])
#                for snp in row[11].split(','):
#                    if snp != 'NONE':
#                        quiery[row[2]].append(snp)
#            except:
#                 pass
#    return(quiery)

def get_check(path, SNP_names_list,clumpig_universe):
    quiery = defaultdict(list)
    not_passed = list()
    qverse =  defaultdict(list)
    with open (path,newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        next(my_reader)
        for row in my_reader:
            try:
                row = list([el.replace('(1)','') for el in row[0].split()])
                qverse[row[2]].append(row[0])
                qverse[row[2]].append(list())
                qverse[row[2]][1].append(row[2])
                quiery[row[2]].append(row[2])
                for snp in row[11].split(','):
                    if snp != 'NONE':
                        quiery[row[2]].append(snp)  
                        qverse[row[2]][1].append(snp)
            except:
                pass                               
    for index in qverse.keys():
        intervals=[]
        try:
            for snp in qverse[index][1]:
                intervals.append(int(chr_snp_dict[qverse[index][0]][snp]))
            qverse[index].append(max(intervals)-min(intervals))
            qverse[index].append(min(intervals))
            qverse[index].append(max(intervals))
        except:
            pass   
    init_count=0 
    for index in quiery.keys():
        for snp in quiery[index]:
            init_count+=1         
#    for snp in old_clusters[path.split('\\')[3].split('.')[0]]:
#        if int(snp[0])==6 and  (int(SNP_names_list[snp[1]].split("_")[1]) <= 33775446 and int(SNP_names_list[snp[1]].split("_")[1])>=28866528 or int(SNP_names_list[snp[1]].split("_")[1]) <= 28700446 and int(SNP_names_list[snp[1]].split("_")[1]) >= 25000528):
#            pass
#        else:
#            flag=0
#            for index in  quiery.keys(): 
#                if snp[1] in quiery[index]: 
#                    flag = 1
#            if flag == 0:
#                try:
#                    check_flag = 0
#                    for clump in qverse.keys():
#                        if qverse[clump][0]==snp[0]:
#                            a = [clumpig_universe[clump][3],clumpig_universe[clump][4]]
#                        if int(a[0]) <= int(chr_snp_dict[snp[0]][snp[1]]) <= int(a[1]):
#                            quiery[clump].append(snp[1])
#                            all_not.append(snp)
#                            check_flag=1                           
#                            break
#                    if check_flag ==0:
#                        all_not.append(snp)
#                        not_passed.append((snp[0],snp[1]))  
#                except:
#                    pass 
    for snp in not_passed:
         for clump in clumpig_universe.keys():
            if len(clumpig_universe[clump])!=0:
                if clumpig_universe[clump][0]==snp[0]:
                    a = [clumpig_universe[clump][3],clumpig_universe[clump][4]]
                    if int(a[0]) <= int(chr_snp_dict[snp[0]][snp[1]]) <= int(a[1]):
                        quiery[clump].append(snp[1])
                        break
    return(qverse)

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))



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
                        if max([int(getOverlap(a, b))/int(clumpig_universe[index1][2]),int(getOverlap(a, b))/int(clumpig_universe[index2][2])])>0.6:
                            ind_check=1
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
                                if max([int(getOverlap(a, b))/int(to_reduce[index1][2]),int(getOverlap(a, b))/int(to_reduce[index2][2])])>=0.6:
                                    intervals[0]=min([int(x) for x in [intervals[0],a[0],b[0]]])
                                    intervals[1]=max([int(x) for x in [intervals[1],a[1],b[1]]])
                                    snp_set.update(to_reduce[index1][1])
                                    snp_set.update(to_reduce[index2][1]) 
                                    passed.add(index2)
            
            ret_universe[index1]=([to_reduce[index1][0],list(snp_set),int(intervals[1])-int(intervals[0]),int(intervals[0]),int(intervals[1])])
        passed.add(index1)
    return(ret_universe)

def GSEA(refset, quiery, genesetdb):
    pval_arrays_new=list()
    N=len(refset)
    n=len(set(quiery))
    nums =list()
    Hla_k=0
    Hist_k=0
    for index in quiery.keys():
            if int(SNP_names_list[index].split("_")[0])==6 and int(SNP_names_list[index].split("_")[1]) <= 33775446 and int(SNP_names_list[index].split("_")[1])>=28866528:
                Hla_k+=1
            if int(SNP_names_list[index].split("_")[0])==6 and int(SNP_names_list[index].split("_")[1]) <= 28000446 and int(SNP_names_list[index].split("_")[1])>=25000528:
                Hist_k+=1
    for gene_set in genesetdb:
        if gene_set=='HLA_locus':
            k=min(len(set(genesetdb[gene_set])), Hla_k)
        else:
            k=min(len(set(genesetdb[gene_set])), Hist_k)
        K=len(set(genesetdb[gene_set]))
        nums.append((N,n,K,k))
        pval_arrays_new.append(hypergeom.sf(k-1, N-k,n, K))
    qval_arrays = list(estimate(np.asarray(pval_arrays_new)))
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
        
#bound_HLA = [28866528, 33775446]
#bound_HIST1 = [26000528, 28000446]

genesetdb = defaultdict(list)

for index in clumpig_universe.keys():
   bp =  int(SNP_names_list[index].split("_")[1])
   chrom = int(SNP_names_list[index].split("_")[0])
   if chrom==6 and  bp <= 33775446 and bp>=28866528:
       genesetdb['HLA_locus'].append(index)
   if chrom==6  and bp <= 28000446 and bp >= 25000528:
       genesetdb['HIST1_complex'].append(index)
            
path = 'C:\\Downloads\\clumped'
print('GSEA')
for filename in os.listdir(path):
    if 'clumped.clumped' in filename and "_phen" not in filename:
#        print(len(get_check(path+"\\"+filename, SNP_names_list,clumpig_universe)))
#        print(len(clump_reducing(clump_reducing(get_check(path+"\\"+filename, SNP_names_list,clumpig_universe)))))
        out_name = 'C:\\Downloads\\Clumped_enrich\\HLA_2\\'+filename+'.tsv'
        df1 = GSEA(refset, clump_reducing(clump_reducing(get_check(path+"\\"+filename, SNP_names_list,clumpig_universe))), genesetdb)
        if len(df1)!=0:
            df1.to_csv(out_name,sep='\t', index=False)
        
cluster_list = []
path1 = 'C:\\Downloads\\GSEA-Fisher-master\\data\\Unions_n2'
for filename in os.listdir(path1):
    cluster_list.append(filename[28:len(filename)-27])    

path = 'C:\\Downloads\\Clumped_enrich\\HLA_2'
pathways_dict = defaultdict(dict)
print('heatmap')
for filename in os.listdir(path):
    with open(path+"\\"+filename,newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')  
        next(my_reader)
        for row in my_reader:
            if str(row[0]) not in pathways_dict:
                pathways_dict[str(row[0])] = defaultdict(list)
            pathways_dict[str(row[0])][filename.split('.')[0]].append(row[1])
            pathways_dict[str(row[0])][filename.split('.')[0]].append(row[2])
#            pathways_dict[str(row[0])][filename[28:len(filename)-8]].append(row[2])    

print('heatmap')

with open ('C:\\Downloads\\clumped_HLA_heatmap_filtered.tsv', 'w',newline='') as csvfile1:
    my_writer = csv.writer(csvfile1, delimiter='\t')
    first_row = copy.copy(cluster_list)
    first_row.insert(0,"category")
    my_writer.writerow(first_row)
    for pathway in pathways_dict:
        row = [None]*30
        row[0] = pathway
        for i in range(1,len(cluster_list)+1):
            if cluster_list[i-1] in pathways_dict[pathway]:
                row[i]=pathways_dict[pathway][cluster_list[i-1]][0]
            else:
                row[i]='NA'
        my_writer.writerow(row)
