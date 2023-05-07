
HEADER_LINE_START = "Seq1,Seq2"

def read_pair_list(list_path):
    query_genomes = []
    reference_genomes = []

    with open(list_path) as f:
        for line in f:
            line = line.strip()
            
            if not line or line[0] == "#" or line.startswith(HEADER_LINE_START):
                continue
            
            genome_pair = line.split(",")
            query_genomes.append(genome_pair[0].strip())
            reference_genomes.append(genome_pair[1].strip())
    
    return query_genomes, reference_genomes


def format_mmseqs_params(params):
    return " ".join([f"-{key} {value}" if len(key) == 1 else f"--{key} {value}" 
                     for key, value in params.items()])


def display_settings(INPUT_FILE, OUTPUT_DIR, TMP_DIR, FRAGMENT_SIZE, CDS_BASED, MEMORY_EFFICIENT, MMSEQS_THREADS, MMSEQS_PARAMS, EVALUE, IDENTITY, COVERAGE):    
    """ Print settings to console """
    
    print('PATHS:')
    print(f'Input file: {INPUT_FILE}')
    print(f'Output directory: {OUTPUT_DIR}')
    print(f'Temporary directory: {TMP_DIR}\n')

    print('PARAMETERS:')
    print(f'Fragment size: {FRAGMENT_SIZE}')
    print(f'CDS based: {CDS_BASED}')
    print(f'Memory efficient mode: {MEMORY_EFFICIENT}\n')

    print(f'MMSEQS threads: {MMSEQS_THREADS}')
    print(f'MMSEQS params: {MMSEQS_PARAMS}\n')

    print(f'Maximum e-value: {EVALUE}')
    print(f'Minimum identity: {IDENTITY}')
    print(f'Minimum query & target coverage: {COVERAGE}')

