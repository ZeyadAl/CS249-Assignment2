import os

tool = 'sweep_dbg'
reads = ['no_error_ont_hq_50x.fastq', 'no_error_reads_hiseq_5k.fastq', 'ont_hq_50x.fastq', 'reads_hiseq_5k.fastq']

dataset = "../synthetic_dataset/"
dataset = "../toy_dataset/"

num = 51
for num in range(25, 71, 5):
    for read in reads:
        os.system(f"./{tool}.py -i ../../../synthetic_dataset/reads/{read} -k {num} -o sweep/k_{read}")
