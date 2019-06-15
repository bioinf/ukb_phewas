# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 14:40:34 2018

@author: anton
"""

import csv

csvwfile = open('C:\\Downloads\\herit_merged_data_for_heatmap.tsv', 'w',newline='')
my_writer = csv.writer(csvwfile, delimiter='\t')
    
csvfile = open ('C:\\Downloads\\herit_data_for_heatmap.tsv')
my_reader = csv.reader(csvfile, delimiter='\t')
row = next(my_reader)
while(True):
    try:
        cur_line = next(my_reader)
        if row[1] != cur_line[1]:
            my_writer.writerow(row)
        else:
            cur_line[0] = row[0] + '; ' + cur_line[0]
        row = cur_line
    except:
        my_writer.writerow(cur_line)
        break
