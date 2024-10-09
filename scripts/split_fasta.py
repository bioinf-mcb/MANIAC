
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

INPUT_PATH = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]
FRAGMENT_SIZE = snakemake.params.FRAGMENT_SIZE 

def split_into_fragments(seq, fragment_length = 1020):
    '''A generator to divide a sequence into fragments of n units.'''
    while seq:
        yield seq[:fragment_length]
        seq = seq[fragment_length:]

with open(OUTPUT_PATH, "w") as output_handle:
    for record in SeqIO.parse(INPUT_PATH, "fasta"):
        for i, fragment in enumerate(split_into_fragments(record.seq, FRAGMENT_SIZE)):
            fragment_record = SeqRecord(fragment, id=record.id + "_FRAGMENT_" + str(i), description="")
            SeqIO.write(fragment_record, output_handle, "fasta")