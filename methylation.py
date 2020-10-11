import pandas
from Bio import SeqIO


def which_chromosomes(bed_file):
    '''returns the number(or letter) for each chromosome in the bed file'''
    dataframe = pandas.read_csv(bed_file, delimiter='\t', header=None)
    dataframe.columns = [
        'chr',
        'start_pos',
        'end_pos',
        'name',
        'score',
        'strand'
    ]
    chromosomes = dataframe.chr.unique()
    crhomosome_numbers = [chromosome[3:] for chromosome in chromosomes
                                    if chromosome != 'chrM' and chromosome != 'rRNA']

    return crhomosome_numbers


def get_extended_sequence(chr_n, extend, bed_file):
    '''
    extend sequences using reference genome.

    :param chr_n: int
    chromosome number or letter

    :param extend: int
    +- base pairs to extend from both sides of the sequence

    :param bed_file: path
    A BED format file containing the mapped sequences.

    :return:
    '''

    chromosome = f'chr{chr_n}'

    # reference genome
    ref_filename = f'reference_genome/{chromosome}.fasta'
    reference = SeqIO.read(ref_filename, "fasta") # a SeqRecord object

    # read bed file into a dataframe
    mapped = pandas.read_csv(bed_file, delimiter='\t', header=None)
    mapped.columns = [
        'chr',
        'start_pos',
        'end_pos',
        'name',
        'score',
        'strand'
    ]

    # sequences mapped to the specific chromosome
    chromosome_mapped = mapped.loc[mapped.chr == chromosome]

    # iterate through the mapped sequences and extend each one by 'extend'
    extended_seq_records = []
    for row in chromosome_mapped.index:
            # get the positions on the reference genome
            start_position = chromosome_mapped.loc[row, 'start_pos']
            end_position = chromosome_mapped.loc[row, 'end_pos']
            # extend the sequence
            extended_sequence = \
                reference[(start_position - extend):(end_position + extend)]
            # append to list
            extended_seq_records.append(extended_sequence)

    return extended_seq_records

def count_sites_and_base_pairs(records):
    '''

    :param records: list
    extended sequences. a list of SeqRecoreds
    :return:
    '''

    # initialize dictionary to store all the results
    counts = {}

    # initialize the count for recognition sites and base pairs
    total_site_count = 0
    total_base_pair_count = 0

    # count
    for record in records:
        sequence = record.seq

        site_count = sequence.count('TCGA')\
                             + sequence.count('AGCT') # count recognition sequence and complementary sequence
        base_pair_count = len(sequence)

        total_site_count += site_count
        total_base_pair_count += base_pair_count
        # total_sequences_count += 1

    counts['sites'] = total_site_count
    counts['base_pairs'] = total_base_pair_count

    return counts

