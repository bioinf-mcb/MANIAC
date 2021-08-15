
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
