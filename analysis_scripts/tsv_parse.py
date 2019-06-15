# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 00:26:14 2018

@author: anton
"""

import csv

csvwfile = open('C:\\Downloads\\herit_merged_09.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')
    
csvfile = open ('C:\\Downloads\\herit_sorted_09.tsv')
my_reader = csv.reader(csvfile, delimiter='\t')
row = next(my_reader)
while(True):
    try:
        cur_line = next(my_reader)
        if row[0] != cur_line[0]:
            my_writer.writerow(row)
        else:
            cur_line[10] = int(row[10])+1
            for i in range (1,10):
                cur_line[i] = row[i] + '; ' + cur_line[i]
        row = cur_line
    except:
        my_writer.writerow(cur_line)
        break
    