
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


def input_checkpoint(INPUT_FILE, SEPARATOR, CDS_BASED):

    print('Processing input... ', end='')

    # verify separator for CDS BASED
    if CDS_BASED:
        print('Veryfiying separator... ', end='')
        with open(INPUT_FILE) as f:
            first_line = f.readline()
        
        if SEPARATOR in first_line: pass
        elif '_CDS' in first_line: 
            print(f'\nWARNING: CHANGING SEPARATOR! "{SEPARATOR}" not found, but "_CDS" found!')
            SEPARATOR = '_CDS'
        elif '_cds' in first_line:
            print(f'\nWARNING: CHANGING SEPARATOR! "{SEPARATOR}" not found, but "_cds" found!')
            SEPARATOR = '_cds'
        else:
            print(f'\nFAILED! "{SEPARATOR}" nor "_CDS" nor "_cds" found in first header!')
            exit()


    # print('File does not exists! Abort!')
    # print('File is empty! Abort!')
    # print('File is not a fasta file!')
    # print('Fasta file is corrupted! [less then 10 characters | less then 2 proteins]')
    # print('You are trying to run protein comparisons, but your file do not contain proteins!')
    # print('You are trying to run nucleotide comparisons, but your file do not contain nucleotides!')
    # print('Extracting phage identifiers... ') 
    # print('Parsing protein IDs... ') # [PHAGEID_"CDS"_NUMBER]


    # print('Phage IDs are not unique!')
    # print('Protein IDs are not unique!')
    

    print('Success!')

    return SEPARATOR

def get_params(config, cores, modes_available = ['FRAGMENTS_NT', 'CDS_NT', 'CDS_AA']):
    """ extract information about parameters from config file """


    # global params
    MODE = nested_get(config, ['MODE'])
    MMSEQS_THREADS = cores                  # mmseqs threads (equivalent to CPU cores)
    FAST = config.get("FAST", False)
    DELETE_INTERMEDIATE_FILES = config.get("DELETE_INTERMEDIATE_FILES", True)
    MEMORY_GB = config.get("MEMORY_GB", 16)

    # checkpoint for mode
    modes_fstring = ' | '.join(modes_available)
    if MODE in modes_available: print(f'MODE: "{MODE}"\t[{modes_fstring}]\nFAST: {FAST}\t[True | False]')
    else: 
        print(f'\nFATAL ERROR! Incorrect mode "{MODE}"\t(FAST: {FAST}.\t Available modes: [{modes_fstring}]\nFAST: [True | False]')
        exit()


    # mmseqs params per mode 
    if MODE == 'CDS_AA':

        # CDS
        CDS_BASED = True
        FRAGMENT_SIZE = None
        SEPARATOR = config.get("SEPARATOR", "_CDS")

        # homologous proteins
        HOMOLOGS_IDENTITY = nested_get(config, ['HOMOLOGS', 'IDENTITY'], default=0.3)
        HOMOLOGS_COVERAGE = nested_get(config, ['HOMOLOGS', 'COVERAGE'], default=0.7)

        # conservative proteins
        CONSERVED_IDENTITY = nested_get(config, ['CONSERVED', 'IDENTITY'], default=0.8)
        CONSERVED_COVERAGE = nested_get(config, ['CONSERVED', 'COVERAGE'], default=0.5)

        # mmseqs params
        EVALUE = nested_get(config, ['MMSEQS_PARAMS', 'EVALUE'], default='1e-15')
        SEARCH_TYPE = 1
        SENSITIVITY = nested_get(config, ['MMSEQS_PARAMS', 'SENSITIVITY'], default=7.5)
        ZDROP = nested_get(config, ['MMSEQS_PARAMS', 'ZDROP'], default=40)
        MAX_SEQS = nested_get(config, ['MMSEQS_PARAMS', 'MAX_SEQS'], default=10000)
        MAX_SEQ_LEN = nested_get(config, ['MMSEQS_PARAMS', 'MAX_SEQ_LEN'], default=65000)
            
        MMSEQS_PARAMS = f"--search-type {SEARCH_TYPE} -a --max-seqs {MAX_SEQS} --max-seq-len {MAX_SEQ_LEN} -s {SENSITIVITY} --mask 0 -e {EVALUE} --zdrop {ZDROP} -c {HOMOLOGS_COVERAGE} --cov-mode 2"


    elif MODE == 'CDS_NT':

        # CDS
        CDS_BASED = True
        FRAGMENT_SIZE = None
        SEPARATOR = config.get("SEPARATOR", "_CDS")

        # homologous proteins
        HOMOLOGS_IDENTITY = nested_get(config, ['HOMOLOGS', 'IDENTITY'], default=0.3)
        HOMOLOGS_COVERAGE = nested_get(config, ['HOMOLOGS', 'COVERAGE'], default=0.7)

        # conservative proteins
        CONSERVED_IDENTITY = nested_get(config, ['CONSERVED', 'IDENTITY'], default=0.8)
        CONSERVED_COVERAGE = nested_get(config, ['CONSERVED', 'COVERAGE'], default=0.5)

        # mmseqs params
        EVALUE = nested_get(config, ['MMSEQS_PARAMS', 'EVALUE'], default='1e-15')
        SEARCH_TYPE = 3
        SENSITIVITY = nested_get(config, ['MMSEQS_PARAMS', 'SENSITIVITY'], default=7.5)
        ZDROP = nested_get(config, ['MMSEQS_PARAMS', 'ZDROP'], default=40)
        MAX_SEQS = nested_get(config, ['MMSEQS_PARAMS', 'MAX_SEQS'], default=10000)
        MAX_SEQ_LEN = nested_get(config, ['MMSEQS_PARAMS', 'MAX_SEQ_LEN'], default=65000)
        KMER = nested_get(config, ['MMSEQS_PARAMS', 'KMER'], default=11)
        SEED_SUB_MATRIX = nested_get(config, ['MMSEQS_PARAMS', 'SEED_SUB_MATRIX'], default='scoring/blastn-scoring.out')
        SUB_MATRIX = nested_get(config, ['MMSEQS_PARAMS', 'SUB_MATRIX'], default='scoring/blastn-scoring.out')

        if FAST:
            KMER = 15

        MMSEQS_PARAMS = f'--search-type {SEARCH_TYPE} -a --max-seqs {MAX_SEQS} --max-seq-len {MAX_SEQ_LEN} -s {SENSITIVITY} --mask 0 -e {EVALUE} -k {KMER} --zdrop {ZDROP} -c {HOMOLOGS_COVERAGE} --cov-mode 2 --seed-sub-mat "{SEED_SUB_MATRIX}" --sub-mat "{SUB_MATRIX}"'



    elif MODE == 'FRAGMENTS_NT':

        # CDS
        CDS_BASED = False
        SEPARATOR = config.get("SEPARATOR", "_FRAGMENT")
        FRAGMENT_SIZE = config.get("FRAGMENT_SIZE", 500)

        # filter significant proteins
        HOMOLOGS_IDENTITY = nested_get(config, ['IDENTITY'], default=0.3)
        HOMOLOGS_COVERAGE = nested_get(config, ['COVERAGE'], default=0.7)

        CONSERVED_IDENTITY = None
        CONSERVED_COVERAGE = None

        # mmseqs params
        EVALUE = nested_get(config, ['MMSEQS_PARAMS', 'EVALUE'], default='1e-15')
        SEARCH_TYPE = 3
        SENSITIVITY = nested_get(config, ['MMSEQS_PARAMS', 'SENSITIVITY'], default=7.5)
        ZDROP = nested_get(config, ['MMSEQS_PARAMS', 'ZDROP'], default=40)
        MAX_SEQS = nested_get(config, ['MMSEQS_PARAMS', 'MAX_SEQS'], default=10000)
        MAX_SEQ_LEN = nested_get(config, ['MMSEQS_PARAMS', 'MAX_SEQ_LEN'], default=65000)
        KMER = nested_get(config, ['MMSEQS_PARAMS', 'KMER'], default=11)
        SEED_SUB_MATRIX = nested_get(config, ['MMSEQS_PARAMS', 'SEED_SUB_MATRIX'], default='scoring/blastn-scoring.out')
        SUB_MATRIX = nested_get(config, ['MMSEQS_PARAMS', 'SUB_MATRIX'], default='scoring/blastn-scoring.out')

        if FAST:
            KMER = 15
            FRAGMENT_SIZE = 1020

        MMSEQS_PARAMS = f'--search-type {SEARCH_TYPE} -a --max-seqs {MAX_SEQS} --max-seq-len {MAX_SEQ_LEN} -s {SENSITIVITY} --mask 0 -e {EVALUE} -k {KMER} --zdrop {ZDROP} -c {HOMOLOGS_COVERAGE} --cov-mode 2 --seed-sub-mat "{SEED_SUB_MATRIX}" --sub-mat "{SUB_MATRIX}"'

    else: 
        print('FATAL MODE ERROR!')
        exit()

    return MODE, FAST, CDS_BASED, FRAGMENT_SIZE, SEPARATOR, DELETE_INTERMEDIATE_FILES, HOMOLOGS_IDENTITY, HOMOLOGS_COVERAGE, CONSERVED_IDENTITY, CONSERVED_COVERAGE, MMSEQS_PARAMS, MEMORY_GB, MMSEQS_THREADS


def display_settings(MODE, INPUT_FILE, OUTPUT_DIR, LOG_DIR, INTERMEDIATE_FILES_DIR, FRAGMENT_SIZE, CDS_BASED, SEPARATOR, MMSEQS_THREADS, MEMORY_GB, MMSEQS_PARAMS, HOMOLOGS_IDENTITY, HOMOLOGS_COVERAGE, CONSERVED_IDENTITY, CONSERVED_COVERAGE, FAST, DELETE_INTERMEDIATE_FILES):    
    """ print settings to console """

    print('\nPATHS:')
    print(f'Input file: {INPUT_FILE}')
    print(f'Output directory: {OUTPUT_DIR}')
    print(f'Log directory: {LOG_DIR}')
    print(f'Intermediate directory: {INTERMEDIATE_FILES_DIR}\n')

    if MODE == 'CDS_AA' or MODE == 'CDS_NT':
        print(f'PARAMETERS:')
        print(f'CDS based: {CDS_BASED}')
        print(f'Separator: {SEPARATOR}')
        print(f'Delete intermediate files and fragment/CDS alignments: {DELETE_INTERMEDIATE_FILES})\n')

        print('Homologs proteins definition:')
        print(f'Minimum identity: {HOMOLOGS_IDENTITY}')
        print(f'Minimum query & target coverage: {HOMOLOGS_COVERAGE}\n')

        print('Conserved proteins definition:')
        print(f'Minimum identity: {CONSERVED_IDENTITY}')
        print(f'Minimum query & target coverage: {CONSERVED_COVERAGE}\n')
    
    elif MODE == 'FRAGMENTS_NT':
        print(f'PARAMETERS ({MODE}):')
        print(f'Fragment size: {FRAGMENT_SIZE}')
        print(f'CDS based: {CDS_BASED}')

        print('DNA significant hits definition:')
        print(f'Minimum identity: {HOMOLOGS_IDENTITY}')
        print(f'Minimum query & target coverage: {HOMOLOGS_COVERAGE}\n')

    else: 
        print('FATAL MODE ERROR!')
        exit()

    print(f'Available memory: {MEMORY_GB}G')
    print(f'MMSEQS CPU cores: {MMSEQS_THREADS}')
    print(f'MMSEQS params: {MMSEQS_PARAMS}\n')
    


def nested_get(d, keys, default=None):
    """
    Recursively get a value from a nested dictionary using a list of keys.
    
    Parameters:
    - d (dict): The dictionary to search.
    - keys (list): A list of keys representing the path to the value.
    - default: A default value if the key path doesn't exist.
    
    Returns:
    - The value if found, else the default.
    """
    for key in keys:
        if isinstance(d, dict):
            d = d.get(key, default)
        else:
            return default
    return d

