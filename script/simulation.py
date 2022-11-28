import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

def polyT(length):
    return 'T'*length
    # do I need to randomize the length?

def reverse_complement(seq):
    '''
    Args: <str>
        queried seq
    Returns: <str>
        reverse_complement seq
    '''
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    letters = \
        [comp[base] if base in comp.keys() else base for base in seq]
    return ''.join(letters)[::-1]

def random_seq(length):
    return ''.join(np.random.choice(['A','G','C', 'T'],length))



    

def artificial_template(BC_list, repeats = None,
                        adaptor='ACTAAAGGCCATTACGGCCTACACGACGCTCTTCCGATCT',
                        TSO='TGTACTCTGCGTTGATACCACTGCTT',
                        output_fn='template.fa'):
    '''
    Generate a barcode 
        BC_list: list of unique BC
        repeats: list of count each BC is simulated
        adaptor/TSO: adapter/TSO sequence to attach before/after the barcodes, 
                            (use 10X 3' solution V3 adaptor by default)
        output_fn: output filename (append to the end of the file if exist)
    output:
        FASTA:
            seqname: barcode sequence
            seq: adapter + BC + 12nt random seq + 15nt polyT + 500 random sequence + TSO
    ''' 
    if repeats is None:
        repeats = [1]*len(BC_list)
   
        
    # v3 adaptor
    with open(output_fn, 'a') as f:
        random_frag = get_random_frag(
                    '/home/ubuntu/vol_data/project/SC_analysis/data/Genome_n_Anno/gencode.v31.transcripts.fa',200)
        read_index = 0
        for BC,count in tqdm(zip(BC_list, repeats)):
            for i in tqdm(range(count), leave=False):

                # seq name
                read_index += 1
                f.write(">" + BC + '_' +  str(read_index) + '\n')

                seq = adaptor+BC+random_seq(12) + polyT(15) + next(random_frag) +TSO
                # seq (randomly reversed) 
                if np.random.choice([False,True]):
                    f.write(seq + '\n')
                else:
                    f.write(reverse_complement(seq) + '\n')


def get_random_frag(ref, frag_len):                    
    # There should be one and only one record, the entire genome:
    records = list(SeqIO.parse(ref, "fasta"))
    num_records = len(records)
    limits = [len(r.seq) for r in records]
    while True:
        i = randint(0,num_records-1)
        record = records[i]
        limit = limits[i]
        if frag_len > limit:
            frag = record.seq
            if 'N' in frag:
                continue
            yield reverse_complement(frag)
        else:
            start = randint(0, limit - frag_len)
            end = start + frag_len
            frag = record.seq[start:end]
            if 'N' in frag:
                continue
            yield reverse_complement(frag)
                  
                    
                    
def main():
    # read BC list
    BC_df = pd.read_csv('../data/SR_bc.csv')
    artificial_template(BC_df.BC, BC_df.counts)
    
    
if __name__=='__main__':
    main()