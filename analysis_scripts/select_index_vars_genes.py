#!/usr/bin/env python3

import sys
import re

clfile, vcf = sys.argv[1:]

index_list = []
with open(clfile, 'r') as cf:
    next(cf)
    for line in cf:
        try:
            index_list.append(line.split()[2])
        except:
            pass

index_coords = []
with open(vcf, 'r') as vf:
    for line in vf:
        if line.startswith('#'):
            continue
        content = line.strip().split('\t')
        if content[0] == '6' and int(content[1]) > 25000528 and int(content[1]) < 33775446:
            continue
        if content[2] in index_list:
            print(content[0] + '\t' + str(int(content[1]) - 5000) + '\t' + str(int(content[1]) + 5000))
