#!/usr/bin/env python
import csv

file = open("bin/data5997x2000/factors_csc.txt")

reader = csv.reader(file, delimiter=',')
csc = []
means = []
variances = []


for row in reader:
    csc.append(row[0] + ' ' + row[1] + ' ' + row[2])
    # if(row[3] != ''):
    #     means.append(row[3])
    #     variances.append(row[4])
    
# mean_file = open("means.txt", mode="w")
# for i in means:
#     mean_file.write("%s\n" % i)

# var_file = open("variances.txt", mode="w")
# for i in variances:
#     var_file.write("%s\n" % i)
    
factors_file = open("factors_csc.txt", mode="w")

for i in csc:
    factors_file.write("%s\n" % i)