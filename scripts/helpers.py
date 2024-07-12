
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

def get_params(config, cores):
    """ extract information about parameters from config file """

    # mode
    if config["CDS_BASED"]:
        try: 
            SEED_SUB_MATRIX = config['MMSEQS_PARAMS']['SEED_SUB_MATRIX']
            MODE = 'CDS_BASED_NUCLEOTIDES'
        except KeyError: MODE = 'CDS_BASED_AMINOACIDS'

    elif not config["CDS_BASED"]:
        MODE = 'FRAGMENT_BASED_NUCLEOTIDES'
    else:
        print(f'ERROR! CDS_BASED param error.')
        exit()

    # params 
    if MODE == 'CDS_BASED_AMINOACIDS':
        
        # homologous proteins
        HOMOLOGS_IDENTITY = config['HOMOLOGS']['IDENTITY']
        HOMOLOGS_COVERAGE = config['HOMOLOGS']['COVERAGE']

        # conservative proteins
        CONSERVED_IDENTITY = config['CONSERVED']['IDENTITY']
        CONSERVED_COVERAGE = config['CONSERVED']['COVERAGE']

        # mmseqs threads (equivalent to CPU cores)
        MMSEQS_THREADS = cores 

        # mmseqs params
        EVALUE = config['MMSEQS_PARAMS']['EVALUE']
        SEARCH_TYPE = config['MMSEQS_PARAMS']['SEARCH_TYPE']
        SENSITIVITY = config['MMSEQS_PARAMS']['SENSITIVITY']
        ZDROP = config['MMSEQS_PARAMS']['ZDROP']
        MAX_SEQS = config['MMSEQS_PARAMS']['MAX_SEQS']
        MAX_SEQ_LEN = config['MMSEQS_PARAMS']['MAX_SEQ_LEN']

        MMSEQS_PARAMS = f"--search-type {SEARCH_TYPE} -a --max-seqs {MAX_SEQS} --max-seq-len {MAX_SEQ_LEN} -s {SENSITIVITY} --mask 0 -e {EVALUE} --zdrop {ZDROP}"


    elif MODE == 'CDS_BASED_NUCLEOTIDES':
        # homologous proteins
        HOMOLOGS_IDENTITY = config['HOMOLOGS']['IDENTITY']
        HOMOLOGS_COVERAGE = config['HOMOLOGS']['COVERAGE']

        # conservative proteins
        CONSERVED_IDENTITY = config['CONSERVED']['IDENTITY']
        CONSERVED_COVERAGE = config['CONSERVED']['COVERAGE']

        # mmseqs threads (equivalent to CPU cores)
        MMSEQS_THREADS = cores 

        # mmseqs params
        EVALUE = config['MMSEQS_PARAMS']['EVALUE']
        SEARCH_TYPE = config['MMSEQS_PARAMS']['SEARCH_TYPE']
        SENSITIVITY = config['MMSEQS_PARAMS']['SENSITIVITY']
        ZDROP = config['MMSEQS_PARAMS']['ZDROP']
        MAX_SEQS = config['MMSEQS_PARAMS']['MAX_SEQS']
        MAX_SEQ_LEN = config['MMSEQS_PARAMS']['MAX_SEQ_LEN']
        KMER = config['MMSEQS_PARAMS']['KMER']
        SEED_SUB_MATRIX = config['MMSEQS_PARAMS']['SEED_SUB_MATRIX']
        SUB_MATRIX = config['MMSEQS_PARAMS']['SUB_MATRIX']

        MMSEQS_PARAMS = f'--search-type {SEARCH_TYPE} -a --max-seqs {MAX_SEQS} --max-seq-len {MAX_SEQ_LEN} -s {SENSITIVITY} --mask 0 -e {EVALUE} -k {KMER} --zdrop {ZDROP} --seed-sub-mat "{SEED_SUB_MATRIX}" --sub-mat "{SUB_MATRIX}"'



    elif MODE == 'FRAGMENT_BASED_NUCLEOTIDES':

        # filter significant proteins
        HOMOLOGS_IDENTITY = config['IDENTITY']
        HOMOLOGS_COVERAGE = config['COVERAGE']

        CONSERVED_IDENTITY = None
        CONSERVED_COVERAGE = None

        # mmseqs threads (equivalent to CPU cores)
        MMSEQS_THREADS = cores 

        # mmseqs params
        EVALUE = config['MMSEQS_PARAMS']['EVALUE']
        SEARCH_TYPE = config['MMSEQS_PARAMS']['SEARCH_TYPE']
        SENSITIVITY = config['MMSEQS_PARAMS']['SENSITIVITY']
        ZDROP = config['MMSEQS_PARAMS']['ZDROP']
        MAX_SEQS = config['MMSEQS_PARAMS']['MAX_SEQS']
        MAX_SEQ_LEN = config['MMSEQS_PARAMS']['MAX_SEQ_LEN']
        KMER = config['MMSEQS_PARAMS']['KMER']
        SEED_SUB_MATRIX = config['MMSEQS_PARAMS']['SEED_SUB_MATRIX']
        SUB_MATRIX = config['MMSEQS_PARAMS']['SUB_MATRIX']

        MMSEQS_PARAMS = f'--search-type {SEARCH_TYPE} -a --max-seqs {MAX_SEQS} --max-seq-len {MAX_SEQ_LEN} -s {SENSITIVITY} --mask 0 -e {EVALUE} -k {KMER} --zdrop {ZDROP} --seed-sub-mat "{SEED_SUB_MATRIX}" --sub-mat "{SUB_MATRIX}"'


    else: 
        print(f'CDS_BASED [True | False] cannot be {CDS_BASED} of {type(CDS_BASED)}! Abort!')
        exit()

    return MODE, HOMOLOGS_IDENTITY, HOMOLOGS_COVERAGE, CONSERVED_IDENTITY, CONSERVED_COVERAGE, MMSEQS_PARAMS, MMSEQS_THREADS


def display_settings(MODE, INPUT_FILE, OUTPUT_DIR, LOG_DIR, INTERMEDIATE_FILES_DIR, FRAGMENT_SIZE, CDS_BASED, MEMORY_EFFICIENT, SEPARATOR, MMSEQS_THREADS, MMSEQS_PARAMS, HOMOLOGS_IDENTITY, HOMOLOGS_COVERAGE, CONSERVED_IDENTITY, CONSERVED_COVERAGE):    
    """ print settings to console """

    print('\nPATHS:')
    print(f'Input file: {INPUT_FILE}')
    print(f'Output directory: {OUTPUT_DIR}')
    print(f'Log directory: {LOG_DIR}')
    print(f'Intermediate directory: {INTERMEDIATE_FILES_DIR}\n')


    if MODE == 'CDS_BASED_AMINOACIDS' or MODE == 'CDS_BASED_NUCLEOTIDES':
        print(f'PARAMETERS ({MODE}):')
        print(f'CDS based: {CDS_BASED}')
        print(f'Memory efficient mode: {MEMORY_EFFICIENT}')
        print(f'Separator: {SEPARATOR}\n')

        print('Homologs proteins definition:')
        print(f'Minimum identity: {HOMOLOGS_IDENTITY}')
        print(f'Minimum query & target coverage: {HOMOLOGS_COVERAGE}\n')

        print('Conserved proteins definition:')
        print(f'Minimum identity: {CONSERVED_IDENTITY}')
        print(f'Minimum query & target coverage: {CONSERVED_COVERAGE}\n')
    
    elif MODE == 'FRAGMENT_BASED_NUCLEOTIDES':
        print(f'PARAMETERS ({MODE}):')
        print(f'Fragment size: {FRAGMENT_SIZE}')
        print(f'CDS based: {CDS_BASED}')
        print(f'Memory efficient mode: {MEMORY_EFFICIENT}\n')

        print('DNA significant hits definition:')
        print(f'Minimum identity: {HOMOLOGS_IDENTITY}')
        print(f'Minimum query & target coverage: {HOMOLOGS_COVERAGE}\n')

    else: 
        print(f'CDS_BASED [True | False] cannot be {CDS_BASED} of {type(CDS_BASED)}! Abort!')
        exit()
    
    
    print(f'MMSEQS CPU cores: {MMSEQS_THREADS}')
    print(f'MMSEQS params: {MMSEQS_PARAMS}\n')
    

