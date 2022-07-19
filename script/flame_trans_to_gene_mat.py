import sys
import argparse
import pandas as pd
import textwrap

def parse_arg():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent(
        '''
        This script can be used to convert the transcipt count matrics from Flame to a gene count
        matrix with similar format (cell barcode for columns and gene for rows)
        '''),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    
    parser.add_argument('in_file', type=str,
                        help='Filename of the transcript count CSV file from FLAME')
    parser.add_argument('--out_file', type=str, default='gene_count.csv',
                        help='Output filename.')
    args = parser.parse_args()
    return args

def main(args):
    d=pd.read_csv(args.in_file)
    d_gene = d.groupby(['gene_id']).sum()
    d_gene.to_csv(args.out_file)

if __name__ == '__main__':
    args = parse_arg()
    main(args)