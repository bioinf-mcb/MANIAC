
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

output_records = []

for record in SeqIO.parse(INPUT_PATH, "fasta"):
    for i, fragment in enumerate(split_into_fragments(record.seq, FRAGMENT_SIZE)):
        fragment_record = SeqRecord(fragment, id = record.id + "_FRAGMENT_" + str(i), description = "")
        output_records.append(fragment_record)

with open(OUTPUT_PATH, "w") as output_handle:
    SeqIO.write(output_records, output_handle, "fasta")
