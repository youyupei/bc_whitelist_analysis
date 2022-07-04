import collections
import scipy.sparse as sp_sparse
import tables
import pysam

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


def get_BC_BQ_in_read(read_seq, bc_start, bc_end):
    '''
    Parameters
    ----------
    read_seq : STR
        Read sequence.
    bc_start : int
        1-based coordinate from sicelore bam.
    bc_end : int
        1-based coordinate from sicelore bam.

    Returns
    -------
    barcode within read
    '''
    if not bc_start or not bc_end:
        return None
    strand = '+'
    if bc_start > bc_end:
        bc_start, bc_end = bc_end, bc_start
        strand = '-'
        return reverse_complement(read_seq[bc_start-1: bc_end]), strand
    return read_seq[bc_start-1: bc_end], strand


def get_BC_in_read_from_adp_pos(read_seq, adp_start, adp_end, flank_size=16):
    '''
    Parameters
    ----------
    read_seq : STR
        Read sequence.
    bc_start : int
        1-based coordinate from sicelore bam.
    bc_end : int
        1-based coordinate from sicelore bam.

    Returns
    -------
    barcode within read
    '''
    if not adp_start or not adp_end:
        return None, None
    strand = '-'
    if adp_start > adp_end:
        strand = '+'
    
    if strand == '-':
    # 1-base to 0-based
        bc_start = adp_end #+ 1 - 1
        bc_end = bc_start+flank_size
        return read_seq[bc_start: bc_end], strand
    if strand == '+':
        bc_end = adp_end-1
        bc_start = bc_end-flank_size
        return reverse_complement(read_seq[bc_start: bc_end]), strand

    else:
        return None, None


def get_BQ_in_read_from_adp_pos(Q_seq, adp_start, adp_end, flank_size=16):
    '''
    Parameters
    ----------
    read_seq : STR
        Read sequence.
    bc_start : int
        1-based coordinate from sicelore bam.
    bc_end : int
        1-based coordinate from sicelore bam.

    Returns
    -------
    barcode within read
    '''
    if not adp_start or not adp_end:
        return None
    strand = '-'
    if adp_start > adp_end:
        strand = '+'
    
    if strand == '-':
    # 1-base to 0-based
        bc_start = adp_end #+ 1 - 1
        bc_end = bc_start+flank_size
        return Q_seq[bc_start: bc_end]
    if strand == '+':
        bc_end = adp_end-1
        bc_start = bc_end-flank_size
        return Q_seq[bc_start: bc_end][::-1]

    else:
        return None
    
    
CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
 
def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
         
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()
         
        return CountMatrix(feature_ref, barcodes, matrix)

def get_3_prime_genome_loc(read, strand):
    '''
    Get the 0-base location of the where the 3' end of the mapped
    
    Input:
        read:
            <class:pysam.AlignedSegment> 
        strand:
            '+' or '-':
                strand of the reads determined by poly A/T
    output:
        <int>
    '''
    if not strand or strand not in '+-':
           return None
    
    if strand == '+':
        # read with poly-A tail at the end of the reads
        if read.is_reverse:
            return read.reference_start
        elif read.is_forward:
            return read.reference_end-1
        else:
            return None
    if strand == '-':
        # read with poly-A tail at the end of the reads
        if read.is_reverse:
            return read.reference_end-1
        elif read.is_forward:
            return read.reference_start
        else:
            return None
        
        
def write_bam_by_read_ids(inbam, outbam, read_ids, chrID = None):
    import pysam
    samfile = pysam.AlignmentFile(inbam, "rb")
    outreads = pysam.AlignmentFile(outbam, "wb", template=samfile)

    for read in tqdm(samfile.fetch(chrID)):
        if read.qname in read_ids:
            outreads.write(read)

    outreads.close()
    samfile.close()
    