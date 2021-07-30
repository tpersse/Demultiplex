#! /usr/bin/env python

import argparse
import gzip
import matplotlib.pyplot as plt

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="name of file to be opened", type=str)
    parser.add_argument("-o", help="output file name and destination", type=str)
    return parser.parse_args()

args = getArgs()

f=args.f
o=args.o

def init_list(lst, cell_num):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    lst = [0]*cell_num
    return lst

my_list= []
#my_list = init_list(my_list)

def convert_phred(ch,num):
    """Converts a single character into a phred score"""
    return ord(ch) - num

def populate_list(f):
    """This function populates the lists created via init_list 
    with the quality score values of each nucleotide at each position
    in the file"""
    new_list1 = []
    i = 0
    #new_list = init_list(new_list1)
    with gzip.open(f, "rt") as fastq:
        for line in fastq:
            line = line.strip('\n')
            i+=1 
            if i%4==0:
                if i == 4:
                    qual_score_len = len(line)
                    new_list = init_list(new_list1, qual_score_len)
                for y, ch in enumerate(line):
                    new_list[y]+=(convert_phred(ch,33))
        return(new_list, i)

my_list, num_lines = populate_list(f)
#print(my_list)
t_list = []
for t, phred_sum in enumerate(my_list):
    t_list.append(t)
    my_list[t] = (phred_sum / (num_lines/4))
    #print(t, '\t', my_list[t])

# print(t_list)
# print(my_list)

plt.bar(t_list, my_list) #bins=t_list[-1]
plt.xlabel("Position in Sequence")
plt.ylabel("Average Quality Score")
plt.title("Average Quality Score at Each Position in Sequence")
plt.savefig(o)
