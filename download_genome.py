import os

import pandas
from Bio import Entrez

# list of tuples with chromosome name and accession_number.'''
input_file = 'GH38_chromosomes.csv'
dataframe = pandas.read_csv(input_file, delimiter='\t')
dataframe =  dataframe[['UCSC-style-name', 'RefSeq-Accn']]
ACCESSION_NUMBERS = list(dataframe.itertuples(index=False, name=None))
ACCESSION_NUMBERS.pop(0) # chr1 allready downloaded

Entrez.email = 'braudelan@gmail.com'


def download_reference_genome():

    for chr_name, accession_numebr in ACCESSION_NUMBERS:

        filename = f'{chr_name}.fasta'
        
        net_handle = Entrez.efetch(
                db="nucleotide", id=accession_numebr, rettype="fasta", retmode="text")
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print("Saved")


if __name__ == '__main__':
    download_reference_genome()
