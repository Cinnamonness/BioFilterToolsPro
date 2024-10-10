def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = '') -> None:
    """
    Convert a multiline FASTA file to a single-line FASTA file.
    Args:
    :param input_fasta: str:
    The input FASTA file path containing
    the names of sequences and sequences to be modified
    to a single line.
    :param output_fasta: str = '':
    The output FASTA file path containing the names
    of sequences and sequences modified to a single line.
    If not provided, a default output file will be created
    in the same directory as the input file.
    return: None: The function writes the converted
    single-line FASTA file.
    """
    if output_fasta == '':
        directory = input_fasta.rsplit('/', 1)[0]
        output_fasta_name = 'output_fasta.fasta'
        output_fasta = f"{directory}/{output_fasta_name}"

    with open(input_fasta, 'r') as f, open(output_fasta, 'w') as o:
        sequence = []
        sequence_name = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sequence_name is not None:
                    output = f"{sequence_name}\n{''.join(sequence)}\n"
                    o.write(output)
                sequence_name = line
                sequence = []
            else:
                sequence.append(line)
        if sequence_name is not None:
            o.write(f"{sequence_name}\n{''.join(sequence)}\n")


def parse_blast_output(input_file: str,
                       output_file: str = '') -> None:
    """
    Parse a BLAST output file and write the sorted
    descriptions to a specified output file.
    Args:
    :param input_file: str:
    The path to the input BLAST output file
    that will be parsed and sorted.
    :param output_file: str:
    The path for the output parsed BLAST file.
    If not provided, a default output file will be
    created in the same directory as the input file.
    :return: None: The function writes the sorted
    descriptions of proteins by QUERY to the output file.
    """
    if output_file == '':
        directory = input_file.rsplit('/', 1)[0]
        output_file_name = 'parse_blast.txt'
        output_file = f"{directory}/{output_file_name}"

    with open(input_file, 'r') as f:
        lines = f.readlines()

    descriptions = []
    in_significant_alignments_block = False
    first_description_found = False

    for line in lines:
        if 'Query #' in line:
            first_description_found = False
            continue
        if 'Sequences producing significant alignments:' in line:
            in_significant_alignments_block = True
            continue
        if in_significant_alignments_block:
            if any(header in line for header in ['Scientific', 'Description',
                                                 'Common', 'Taxid', 'Query', 'Accession']):
                continue
            if line.strip() == '':
                in_significant_alignments_block = False
                continue
            if not first_description_found:
                if '[' in line:
                    description = line[:line.find('[')].strip()
                elif '...' in line:
                    description = line[:line.find('...')].strip()
                else:
                    description = line.strip()
                descriptions.append(description)
                first_description_found = True

    sorted_descriptions = sorted(descriptions)

    with open(output_file, 'w') as o:
        for description in sorted_descriptions:
            o.write(f"{description}\n")


def select_genes_from_gbk_to_fasta(input_gbk: str,
                                   *genes: str,
                                   n_before: int = 1,
                                   n_after: int = 1,
                                   output_fasta: str = '') -> None:
    """
    Select genes of interest from a genome annotation file (.gbk file)
    and write their names and translation sequences into a FASTA file.
    Args:
    :param input_gbk: str:
    The path to the input genome annotation file from which
    genes will be selected.
    :param *genes: str:
    Genes of interest; neighbors will be selected around these genes.
    :param n_before: int = 1:
    Number of genes to select before each gene of interest (>= 1).
    Default is 1.
    :param n_after: int = 1:
    Number of genes to select after each gene of interest (>= 1).
    Default is 1.
    :param output_fasta: str = '':
    The path for the output FASTA file.
    If not provided, defaults to the same directory as the input file.
    :return: None: The function writes the selected genes and their
    translation sequences to the output FASTA file.
    """
    with open(input_gbk, 'r') as f:
        lines = f.readlines()

    genes_on_gbk = [s.split('"')[1] for s in lines if '/gene=' in s]
    selected_genes = []

    for gene in genes:
        gene_index = genes_on_gbk.index(gene)
        for up in range(1, n_before + 1):
            if gene_index - up >= 0:
                selected_genes.append(genes_on_gbk[gene_index - up])
        for down in range(1, n_after + 1):
            if gene_index + down < len(genes_on_gbk):
                selected_genes.append(genes_on_gbk[gene_index + down])

    translation_seqs = []
    for selected_gene in selected_genes:
        indexes = [i for i, line in enumerate(lines) if selected_gene in line]
        for index in indexes:
            start_pos = index
            while '/translation=' not in lines[start_pos]:
                start_pos += 1
            translation_seq = lines[start_pos].split('"')[1].replace('\n', '')
            end_pos = start_pos + 1
            while '"' not in lines[end_pos]:
                translation_seq += lines[end_pos].replace('\n', '').replace(' ', '')
                end_pos += 1
            translation_seq += lines[end_pos].split('"')[0].replace(' ', '')
            translation_seqs.append(translation_seq)

    if output_fasta == '':
        directory = input_gbk.rsplit('/', 1)[0]
        output_fasta_name = 'gbk.fasta'
        output_fasta = f"{directory}/{output_fasta_name}"

    with open(output_fasta, 'w') as o:
        for idx in range(len(selected_genes)):
            o.write(f">{selected_genes[idx]}\n{translation_seqs[idx]}\n")
