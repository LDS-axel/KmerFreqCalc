# -*- coding:utf-8 -*-
import sys
import re
import os
import time
import multiprocessing

#Create directory at present workplace if this directory does not exist
def crtdir(work_station, filename):
    if not os.path.exists(work_station + '/' + filename):
        cmd = 'mkdir ' + work_station + '/' + filename
        os.system(cmd)

#Obtain Inverted complement Kmer of a DNA Kmer
def reverse_complement(s):
    basecomplement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)[::-1]

#Inverted complement Kmer count as default
def Kmercount(freq, Kmer):
    if Kmer in freq:
        freq[Kmer] += 1
    elif reverse_complement(Kmer) in freq:
        freq[reverse_complement(Kmer)] += 1
    else:
        freq[Kmer] = 1
    return freq

#The total number of kmers
def length(freq):
    Kmer_num = 0
    for i in freq:
        Kmer_num += freq[i]
    return Kmer_num

#Summarize the number of kmers of subsequences
def Kmer_statistics(seq_dict):
    freq = {}
    for i in seq_dict:
        my_seq = seq_dict[i].upper()
        for t in range(len(my_seq)-K+1):
            Kmer = my_seq[t:t+K]
            deviation = K-Kmer.count('A')-Kmer.count('T')\
            -Kmer.count('C')-Kmer.count('G')
            if deviation == 0:
                freq = Kmercount(freq, Kmer)
    return freq

#General statistics of Kmer number of all subsequences
def integrate(sub_freq, integrated_freq):
    for Kmer in sub_freq:
        forward, backward = str2int(Kmer), str2int(reverse_complement(Kmer))
        if forward in integrated_freq:
            if forward < backward:
            	integrated_freq[forward] += sub_freq[Kmer]
            else:
                integrated_freq[backward] += sub_freq[Kmer]
        elif backward in integrated_freq:
            if forward < backward:
                integrated_freq[forward] += sub_freq[Kmer]
            else:
                integrated_freq[backward] += sub_freq[Kmer]
        else:
            if forward < backward:
                integrated_freq[forward] = sub_freq[Kmer]
            else:
                integrated_freq[backward] = sub_freq[Kmer]
    return integrated_freq

#General statistics of Kmer number of all subsequences
def summary(sub_freq):
    freq = {}
    for item in sub_freq:
        freq = integrate(item, freq)
    return freq

def str2int(Kmer):
    nucleotide = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    Kmer_list = list(Kmer)
    location = len(Kmer)-1
    Kmer_value = 0
    for i in Kmer_list:
        Kmer_value += nucleotide[i]*4**location
        location -= 1
    return Kmer_value


if __name__ == '__main__':
    help = "python3 Kmer_frequence_vector.py -n -fna.path -speciesname -CPU.NUM\n\
-n\t-num\tthe length of Kmer\tn>1\n\
-fna.path\t-str\tthe path of dna sequence file as an input, presented suffix with fna,fa,ffn,or fasta\n\
-speciesname\t-str\tthe name of species\n\
-CPU.NUM\t-num\tthe number of CPU to be used\n"

    # Input Parameter:
    if len(sys.argv) != 5:
        sys.exit(help)
    K, filement, speciesname, coreNum = \
    int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
    startTime = time.time()

    print('Read the sequences into a dictionary from the sequence filement!!!')
    seq_dict = {}
    for row in open(filement, 'r'):
        if row[0] == '>':
            key = row[1:-1]
            seq_dict[key] = []
        else:
            seq_dict[key].append(row.strip())
    for key, value in seq_dict.items():
        seq_dict[key] = ''.join(value)

    queue = [{} for i in range(coreNum)]
    i = 0
    for seqid in seq_dict:
        queue[i][seqid] = seq_dict[seqid]
        i += 1
        if i == coreNum:
            i = 0
    p = multiprocessing.Pool(coreNum)
    freq_lst = p.map(Kmer_statistics, queue)
    p.close()
    p.join()
    
    print("Make a summary statistic for all kmers' number in all subsequences")
    grass_freq = summary(freq_lst)
    
    index_lst = []
    for i in grass_freq:
        index_lst.append(i)
    index_lst.sort()
    inner_product = 0
    for i in grass_freq:
        inner_product += grass_freq[i]**2
    
	# Formatted output
    kmerNum = len(index_lst)
    output = str(K)+'\n'+str(kmerNum)+'\n'+str(inner_product)+'\n'
    for i in index_lst:
        output += str(i)+' '+str(grass_freq[i])+'\n'
    crtdir('./', 'NUM_K_'+str(K))
    open('NUM_K_'+str(K)+'/'+speciesname+'.cv.txt', 'w').write(output)
    print(time.time()-startTime)
