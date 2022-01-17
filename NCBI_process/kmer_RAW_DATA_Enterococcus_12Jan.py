from itertools import product
from Bio import SeqIO
import sys
import pandas as pd
import os
import glob
import re
import time

def generate_kmer(m,n):
    kmer = []
    for k in range(m,n+1):
        kmer.extend([''.join(i) for i in product('ATGC',repeat=k)])
    return kmer

def calculate_freq(s,kmer):
    freq = []
    for i in kmer:
        freq.append(s.count(i))
    return freq


organisms = ["Enterococcus"]   #["Acinetobacter_baumannii"] #["Escherichia_coli"] #["Klebsiella_pneumoniae"] #"Escherichia_coli",
#organisms = ["Acinetobacter_baumannii", "Escherichia_coli", "Klebsiella_pneumoniae", "Pseudomonas_aeruginosa","Staphylococcus_aureus","Enterococcus",
#                         "Enterococcus_faecalis","Enterococcus_faecium"]

seq_dir = "Sequence_files"
if os.path.exists(seq_dir) == False:
    os.mkdir(seq_dir)
test_dir = "TEST_Kmer"
if os.path.exists(test_dir) == False:
    os.mkdir(test_dir)
for organism in organisms:
    start_time = time.time()
    sequence = []
    filenames = []

    seq_path = os.path.join(seq_dir,organism)
    if os.path.exists(seq_path) == False:
        os.mkdir(seq_path)
    
    paths = os.listdir(organism)
    count = 1
    for p in paths:
        seq = ">" + str(count) +"\n"
        #seq2 = ">" +str(count) +"\n"
        path = organism + "/*.fastq"
        data = glob.glob(path)
        for d in data:
            f=open(d,"r+")
            line = f.read(-1)
            f.close()
            line = line.replace(" ","")
            seq_temp = re.split('(.*)length=\d\d\d\n',line)[1:]
            #special_char = re.compile('[@_!#$%^&*()<>?/\|}{~:+]')
            #other_char = re.compile('[^ATCGN]')

            for i in range(1,len(seq_temp),2):
                if seq_temp[i-1].startswith("@"):
                    seq += seq_temp[i].replace("\n","")

            
            temp_file = d.split("\\")[-1]
            seq_file = os.path.join(seq_path,temp_file)
            f1 = open(seq_file, "w+")
            f1.write(seq)
            f1.close()

            '''    
            for i in seq_temp:
                i = i.replace("\n","")
                if(special_char.search(i) == None and other_char.search(i) == None):
                    seq2 += i
            '''
        sequence.append(seq)
        filenames.append(d.split("\\")[-1])

    kmer = generate_kmer(6,7) #length 6 to 7
    df = pd.DataFrame([calculate_freq(i,kmer) for i in sequence])
    df.columns = kmer
    df.insert(0,"Organism",organism, True)
    df.insert(1,"Filename",filenames,True)
    file_path = os.path.join(test_dir,organism+"_kmer.csv")
    df.to_csv(file_path, index = False)
    print("Time taken to create one feature set of ", organism, " is :", time.time() - start_time)

