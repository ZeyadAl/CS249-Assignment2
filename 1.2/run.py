import os
ont_reads = ['no_error_ont_hq_50x.fastq', 'ont_hq_50x.fastq']
hi_reads = ['no_error_reads_hiseq_5k.fastq', 'reads_hiseq_5k.fastq']

for read in hi_reads:
        os.system(f"/usr/bin/time python3 olc.py ../../2/synthetic_dataset/reads/{read} -o contigs_{read} -n 10")

for read in ont_reads:
        os.system(f"/usr/bin/time python3 olc.py --ont ../../2/synthetic_dataset/reads/{read} -o contigs_{read} -n 100")
