# -*- coding: utf-8 -*-
"""
Created on Thu May  2 23:22:15 2019

@author: anton
"""

from collections import defaultdict
from scipy import interpolate
from scipy.stats import hypergeom
from statistics import median
import csv
import pandas as pd
import scipy as sp
import numpy as np
import os
import copy
import sys
import argparse
import subprocess

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

def make_vcf(SNP_dict,name):
    with open(os.getcwd() + "/"+name.split('.')[0] +".vcf" ,"w",newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter='\t')
        for key in SNP_dict:
            row_to_write = []
            for i in range(8):
                row_to_write.append('.')
            row_to_write[0] = SNP_dict[key][0]
            row_to_write[1] = SNP_dict[key][1]
            row_to_write[2] = key
            row_to_write[3] = SNP_dict[key][2]
            row_to_write[4] = SNP_dict[key][3]
            my_writer.writerow(row_to_write)

def get_check(path, HLA, SNP_names_list,clumpig_universe):
    quiery = defaultdict(list)
    not_passed = list()
    qverse =  defaultdict(list)
    all_not = list()
    if HLA == "yes":
        with open (path,newline='') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            next(my_reader)
            for row in my_reader:
                try:
                    row = list([el.replace('(1)','') for el in row[0].split()])
                    quiery[row[2]].append(row[2])
                    for snp in row[11].split(','):
                        if snp != 'NONE':
                            quiery[row[2]].append(snp)
                except:
                    pass
    if HLA == "no":
        with open (path,newline='') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            next(my_reader)
            for row in my_reader:
                try:
                    row = list([el.replace('(1)','') for el in row[0].split()])
                    flag = 0
                    if int(SNP_names_list[row[2]].split("_")[0])==6 and  (int(SNP_names_list[row[2]].split("_")[1]) <= 33775446 and int(SNP_names_list[row[2]].split("_")[1])>=28866528 or int(SNP_names_list[row[2]].split("_")[1]) <= 28700446 and int(SNP_names_list[row[2]].split("_")[1]) >= 25000528):
                        flag=1
                    if flag != 1:
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
            
    for snp in not_passed:
         for clump in clumpig_universe.keys():
            if len(clumpig_universe[clump])!=0:
                if clumpig_universe[clump][0]==snp[0]:
                    a = [clumpig_universe[clump][3],clumpig_universe[clump][4]]
                    if int(a[0]) <= int(chr_snp_dict[snp[0]][snp[1]]) <= int(a[1]):
                        quiery[clump].append(snp[1])
                        break
    count=0 
    for index in quiery.keys():
        for snp in quiery[index]:
            count+=1  
    return(quiery)


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
            if snp in coord_gene_dict. keys() and" not_set" not in coord_gene_dict[snp][1]:
                indexes = int(coord_gene_dict[snp][1])
                for i in range(indexes):
                    index_set.add(coord_gene_dict[snp][i+2])
                    genes[index].append(coord_gene_dict[snp][i+2])  
        if len(index_set)!=0 and len(passed_genes.union(index_set))>len(passed_genes):
            n_adj+=1
            passed_genes = passed_genes.union(index_set)
    n = n_adj                             
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
        pval_arrays_new.append(hypergeom.sf(k-1, N-k,n, K))
    qval_arrays = list(estimate(np.asarray(pval_arrays_new)))
    results=defaultdict(list)
    i=0
    for geneset in clumpsetbd:
        if qval_arrays[i]<=0.1:
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

def create_clumpverse(flag, chr_snp_dict, coord_chr_snp_dict,all_snps,cl_dir='no'):
    clumpig_universe = defaultdict(list) 
    clumped_snps = set() 
    NONE_clumps = defaultdict(list) 
    if cl_dir=='no':
        cl_dir=os.getcwd()
    for filename in os.listdir(cl_dir):
        if flag in filename:
            if cl_dir=='no':
                name=filename
            else:
                name=cl_dir+"/"+filename
            with open (name,newline='') as csvfile:
                my_reader = csv.reader(csvfile, delimiter='\t')
                next(my_reader)
                for row in my_reader:
                    try:
                        row = list([el.replace('(1)','') for el in row[0].split()])
                        clumped_snps.add(row[2])
                        for snp in row[11].split(','):
                                clumped_snps.add(snp)
                        if row[2] not in clumpig_universe and "NONE" not in row[11].split(','):
                            clumpig_universe[row[2]].append(row[0])
                            clumpig_universe[row[2]].append(list())
                            clumpig_universe[row[2]][1].append(row[2])
                            for snp in row[11].split(','):
                                clumpig_universe[row[2]][1].append(snp)
                        elif row[2] in clumpig_universe and "NONE" not in row[11].split(','):
                            for snp in row[11].split(','):
                                if snp not in clumpig_universe[row[2]][1]:
                                    clumpig_universe[row[2]][1].append(snp)  
                        elif row[2] not in clumpig_universe and "NONE" in row[11].split(',') and row[2] not in NONE_clumps:
                            NONE_clumps[row[2]].append(row[0])
                            NONE_clumps[row[2]].append(list())
                            NONE_clumps[row[2]][1].append(row[2])
                            for snp in row[11].split(','):
                                NONE_clumps[row[2]][1].append(snp)          
                    except:
                        pass
    for index in clumpig_universe.keys():
        intervals=[]
        try:
            for snp in clumpig_universe[index][1]:
                intervals.append(int(chr_snp_dict[clumpig_universe[index][0]][snp]))
            clumpig_universe[index].append(max(intervals)-min(intervals))
            clumpig_universe[index].append(min(intervals))
            clumpig_universe[index].append(max(intervals))
        except:
            pass
    clump_med = median(list([clumpig_universe[index][2] for index in clumpig_universe.keys() if len(clumpig_universe[index][1])>=2 and len(clumpig_universe[index])>3]))  
    for index in NONE_clumps:
        try:
            intervals = []
            test_clump = copy.deepcopy(NONE_clumps[index])
            test_clump[1].remove('NONE')
            intervals.append(int(chr_snp_dict[test_clump[0]][test_clump[1][0]]))
            for coord in coord_chr_snp_dict[test_clump[0]].keys():
                if int(chr_snp_dict[test_clump[0]][test_clump[1][0]]) - clump_med/2 <= int(coord) <= int(chr_snp_dict[test_clump[0]][test_clump[1][0]]) + clump_med/2 and coord != chr_snp_dict[test_clump[0]][test_clump[1][0]]:
                    intervals.append(int(coord))
                    test_clump[1].append(coord_chr_snp_dict[test_clump[0]][coord])
            test_clump.append(max(intervals)-min(intervals))
            test_clump.append(min(intervals))
            test_clump.append(max(intervals))
            check_flag = 0
        
            for clump in clumpig_universe.keys():
                if clumpig_universe[clump][0]==test_clump[0]:
                    a = [clumpig_universe[clump][3],clumpig_universe[clump][4]]
##                b = [test_clump[3], test_clump[4]]
##                if getOverlap(a, b)!=0:
##                    check_flag=1
                    if a[0] <= int(chr_snp_dict[test_clump[0]][index]) <= a[1]:
                         check_flag=1
            if check_flag ==0:
                clumpig_universe[index]=test_clump
        except:
            pass
    unclumped = all_snps.difference(clumped_snps)
    for index in unclumped:
        try:
            intervals = []
            test_clump = [coord_gene_dict[index][0].split("_")[0], [index]]
            intervals.append(int(chr_snp_dict[test_clump[0]][test_clump[1][0]]))
            for coord in coord_chr_snp_dict[test_clump[0]].keys():
                if int(chr_snp_dict[test_clump[0]][test_clump[1][0]]) - clump_med/2 <= int(coord) <= int(chr_snp_dict[test_clump[0]][test_clump[1][0]]) + clump_med/2 and coord != chr_snp_dict[test_clump[0]][test_clump[1][0]]:
                    intervals.append(int(coord))
                    test_clump[1].append(coord_chr_snp_dict[test_clump[0]][coord])
            test_clump.append(max(intervals)-min(intervals))
            test_clump.append(min(intervals))
            test_clump.append(max(intervals))
            check_flag = 0
            for clump in clumpig_universe.keys():
                if clumpig_universe[clump][0]==test_clump[0]:
                    a = [clumpig_universe[clump][3],clumpig_universe[clump][4]]
##                b = [test_clump[3], test_clump[4]]
##                if getOverlap(a, b)!=0:
##                    check_flag=1
                    if a[0] <= int(chr_snp_dict[test_clump[0]][index]) <= a[1]:
                         check_flag=1
            if check_flag ==0:
                clumpig_universe[index]=test_clump
        except:
            pass
    med_verse = defaultdict(list)
    for index in clumpig_universe:
        if len(clumpig_universe[index])>=4:
            med_verse[index]=clumpig_universe[index]
    new_verse = clump_reducing(med_verse) 
    non_hla = defaultdict(list)
    for clump in new_verse: 
        l1=len(new_verse[clump][1])
        l2=0
        for snp in new_verse[clump][1]:
            try:
                if int(SNP_names_list[snp].split("_")[0])==6 and  (int(SNP_names_list[snp].split("_")[1]) <= 33775446 and int(SNP_names_list[snp].split("_")[1])>=28866528 or int(SNP_names_list[snp].split("_")[1]) <= 28700446 and int(SNP_names_list[snp].split("_")[1]) >= 25000528):
                    l2+=1
            except:
                pass
        if l2!=0:
            pass
        else:
            non_hla[clump]=new_verse[clump]
    return(non_hla)

def make_clumps(input_dict,file_name, plink_dir, bfile):
    print('Calculating independant loci with PLINK')
    with open(os.getcwd()+ "/"+file_name.split('.')[0]+'_for_plink.tsv','w', newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter='\t')
        init_row = ["SNP","Chr","Pos","P"]
        my_writer.writerow(init_row)
        for snp in input_dict:
            row = [snp] + input_dict[snp][0:2]+[input_dict[snp][4]]
            my_writer.writerow(row)
        snp_cmd = subprocess.call('{0}/plink --bfile {1} --clump {2} --clump-field P --clump-p1 5e-08 --clump-r2 0.1 --clump-snp-field SNP --clump-kb 500 --out {3} --silent 2> {4}'.format(plink_dir,plink_dir +'/'+ bfile,os.getcwd()+ "/"+file_name.split('.')[0]+'_for_plink.tsv',os.getcwd()+ "/"+file_name.split('.')[0]+'.clumped', 'PLINK_clumping.log'),shell=True)

def annotate_vcf(file_name, snpeff, genome):
    print('Annotating file with snpEff')
    snp_cmd = subprocess.call('java -Xmx4g -jar {0} -v {1} {2} > {3} 2> {4}'.format(snpeff+'/snpEff.jar',genome,os.getcwd() + "/"+file_name,os.getcwd()+"/"+file_name.split('.')[0]+"_annotated.vcf", 'snpeff.log'),shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LSEA')
    parser.add_argument('-af', help='Input file in tsv-format', metavar='File',
                        type=str, required=True)
    parser.add_argument('-vf', help='Annotated vcf-file if it is already prepared ', metavar='File',
                        type=str, required=False)
    parser.add_argument('-pl', help='PlINK  result of clumping (.clumped file) if it is already prepared', metavar='File',
                        type=str, required=False)
    #parser.add_argument('-cf', help='Do not use this argument, Full Manhattan cluster file', metavar='File',
    #                    type=str, required=False)
    #parser.add_argument('-fa', help='Do not use this argument, Full association file', metavar='File',
    #                    type=str, required=False)
    #parser.add_argument('-va', help='Do not use this argument, Full VCF annotation file', metavar='File',
    #                    type=str, required=False)
    parser.add_argument('--precompiled', '-p',action='store_true',help='Use precompiled loci')
    parser.add_argument('-sn', help='Path to SnpEff directory', metavar='File',
                        type=str, required=True)
    parser.add_argument('-g', help='Flag for specifying genome for SnpEff annotation', metavar='Str',
                        type=str, required=False, default='GRCh37.75')
    parser.add_argument('-pld', help='Path to PLINK directory', metavar='Dir',
                        type=str, required=True)
    parser.add_argument('-bf', help='Bfile for PLINK', metavar='Str',
                        type=str, required=False, default='1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01')
    parser.add_argument('-cld', help='directory with clumped files for creating universe', metavar='Dir',
                        type=str, required=False, default='')
    parser.add_argument('--create_universe', '-cu',action='store_true',help='Create universe for GWAS data')
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    #af,p,vf,pl,cf,fa,va,sn,g,pld,bf,cu,cld =  args.af,args.precompiled, args.vf, args.pl,args.cf,args.fa,args.va,args.sn,args.g,args.pld,args.bf,args.create_universe,args.cld
    af,p,vf,pl,sn,g,pld,bf,cu,cld =  args.af,args.precompiled, args.vf, args.pl,args.sn,args.g,args.pld,args.bf,args.create_universe,args.cld
    all_snps = set()   
    clumpsetbd = defaultdict(list)
    genesetdb = defaultdict(list)
    coord_gene_dict = defaultdict(list)
    chr_snp_dict = defaultdict(dict)
    coord_chr_snp_dict = defaultdict(dict)        
    refset=set()
    clumpig_universe = defaultdict(list) 
    SNP_names_list = dict()
    input_dict=defaultdict(list)

    with open(os.getcwd()+'/precompiled_data/all_human_genes.txt',newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        for row in my_reader:
            if len(row)==1:
                refset.add(row[0]) 
    with open(os.getcwd()+'/precompiled_data/all_sets.gmt',newline='') as csvfile:
        my_reader = csv.reader(csvfile, delimiter='\t')
        next(my_reader)
        for row in my_reader:
            elems = (list)
            for el in range(2, len(row)):
                refset.add(row[el])  
    comp_names = list([gene for gene in refset if "-" in gene])   

    with open(os.getcwd()+'/precompiled_data/all_sets.gmt',newline='') as csvfile:
       my_reader = csv.reader(csvfile, delimiter='\t')
       next(my_reader)
       for row in my_reader:
            elems = (list)
            for el in range(2, len(row)):
                genesetdb[row[0]].append(row[el])

    if p:
        with open(os.getcwd()+'/precompiled_data/clumpverse_no_hla.tsv', 'r', newline='') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            for row in my_reader:
                clumpig_universe[row[0]].append(row[1])
                clumpig_universe[row[0]].append(row[2].split(','))
                clumpig_universe[row[0]].append(row[3])
                clumpig_universe[row[0]].append(row[4])
                clumpig_universe[row[0]].append(row[5])

        with open(os.getcwd()+'/precompiled_data/clumpsetdb_filtered.tsv',newline='') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            for row in my_reader:
                clumpsetbd[row[0]].extend(row[1].split(","))    

#    if not cf and not fa and not va:
#        with open(af,'r', newline='') as csvfile:
#            my_reader = csv.reader(csvfile, delimiter='\t')
#            for row in my_reader:
#                all_snps.add(row[2])
#                chr_snp_dict[row[0]][row[2]]=row[1]
#                coord_chr_snp_dict[row[0]][row[1]]=row[2]
#                coord_gene_dict[row[2]].append(row[0]+"_"+row[1])
#                coord_gene_dict[row[2]].append(str(0))
#                if not vf:
#                    input_dict[row[2]]=row[0:2]+row[3:]

    if vf:
        with open(vf,'r', newline='') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            for row in my_reader:
                if "#" not in row[0]: 
                    SNP_names_list[row[2]] = "_".join([row[0]]+[row[1]]+[row[3]]+[row[4]])+'_b37'

        with open(vf,'r', newline='') as csvfile:
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
    if not vf and not pl and p:
        make_vcf(input_dict,af)
        annotate_vcf(af.split('.')[0]+".vcf", sn, g)
        make_clumps(input_dict,af,pld,bf)
        with open(af.split('.')[0]+"_annotated.vcf",'r', newline='') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            for row in my_reader:
                if "#" not in row[0]:
                    SNP_names_list[row[2]] = "_".join([row[0]]+[row[1]]+[row[3]]+[row[4]])+'_b37'

        with open(af.split('.')[0]+"_annotated.vcf",'r', newline='') as csvfile:
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

        df1 = GSEA(clumpig_universe,get_check(os.getcwd()+"/"+af.split('.')[0]+'.clumped.clumped', 'no', SNP_names_list,clumpig_universe),clumpsetbd, genesetdb, coord_gene_dict,'yes')
        if len(df1)!=0:
            df1.to_csv(af.split('.')[0]+'.enrich.tsv',sep='\t', index=False)
            rename_cmd =subprocess.call("var=\"gene_set\tp-value\tq-value\tenrich_description\";  sed -i \"1s/.*/$var/\" {}".format(af.split('.')[0]+'.enrich.tsv'),shell=True)
    if cu:
        if cld:
            clumpig_universe = create_clumpverse('clumped.clumped', chr_snp_dict, coord_chr_snp_dict,all_snps,cl_dir=cld)
        else:
            clumpig_universe = create_clumpverse('clumped.clumped', chr_snp_dict, coord_chr_snp_dict,all_snps)
        #print(clumpig_universe)
    if pl:
        df1 = GSEA(clumpig_universe,get_check(pl, 'no', SNP_names_list,clumpig_universe),clumpsetbd, genesetdb, coord_gene_dict,'yes')
        if len(df1)!=0:
            df1.to_csv(pl.split('.')[0]+'.enrich.tsv',sep='\t', index=False)
            rename_cmd =subprocess.call("var=\"gene_set\tp-value\tq-value\tenrich_description\";  sed -i \"1s/.*/$var/\" {}".format(af.split('.')[0]+'.enrich.tsv'),shell=True)

#"var=\"gene_set\tp-value\tq-value\tenrich_description\";echo $var;  sed -i '1s/.*/$var/' {}"
 
