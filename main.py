import pandas as pd
from tabulate import tabulate

from methylation import *

EXTEND_BY = 200
INPUT_FILES = ['template.bed', 'notemplate.bed']
CHROMOSOME_NUMBERS = which_chromosomes('template.bed')

def write_result(data):
    '''helper function to write the result.'''
    df = pd.DataFrame(data, index=['recognition_sites', 'base_pairs', 'ratio'])
    df.to_excel('output.xlsx')

# count total sites and total base pairs in each 'bed' file
total_counts  = {}
for bed_file in INPUT_FILES:

    # initialize count
    total_sites = 0
    total_base_pairs = 0

    # count in each chromosome
    for number in CHROMOSOME_NUMBERS:
        records = get_extended_sequence(number, EXTEND_BY, bed_file)
        counts = count_sites_and_base_pairs(records)

        n_sites = counts['sites']
        n_base_pairs = counts['base_pairs']

        total_sites += n_sites
        total_base_pairs += n_base_pairs

    # compute ratio
    ratio = total_sites / total_base_pairs * 100

    # write to dict
    file_name = bed_file.split('.')[0]
    total_counts[file_name] = (total_sites, total_base_pairs, ratio)


# write result to an excel file
write_result(total_counts)


## write the result to output file as a table
# df = pd.DataFrame(total_counts, index=['recognition_sites', 'base_pairs', 'ratio'])
# with open('output', 'w') as file:
#     file.write(tabulate(df.round(3), headers='keys', tablefmt='psql'))
