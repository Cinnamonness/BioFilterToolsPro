from modules import module_rna_dna_tools as tools
from modules import module_filter_fastq as module_filter
import os
from typing import Union

"""
The utility for working with DNA and RNA and of
filtering FASTQ sequences.
This utility provides two main functions:
1. Operations on the transmitted DNA or RNA sequences:
transcription, reversing, complementing, reversing
complementing, identifying the molecule type.
2. Filtering FASTQ sequences in FASTQ file
based on GC-content, length and quality score.

Modules:
1. module_rna_dna_tools: contains the operations on the
transmitted DNA or RNA sequences.
2. module_filter_fastq: contains functions to calculate
sequence properties: GC-content, length
and quality score.
"""


def run_dna_rna_tools(*arguments):
    """
    This utility performs one of the operations on
    the transmitted DNA or RNA sequences.
    Supported operations:
    - 'transcribe' - determine transcribed sequence.
    - 'reverse' - reverse the sequence.
    - 'complement' - determine complement sequence.
    - 'reverse_complement' - determine reverse complement sequence.
    Addition operation:
    - 'molecule_type' - determine the type of molecule
    Args:
    :param arguments: DNA/RNA sequences and the operation to be performed.
    last argument: str: determine the operation to be performed.
    other arguments: str: sequences.
    :return: modified sequences depending on the selected operation.
     If one sequence was passed, one line is returned.
     If several, a list of modified sequences is returned.
    """
    *sequences, operation = arguments
    res = []
    for seq in sequences:
        if tools.molecule_type(seq) in ['DNA', 'RNA']:
            if operation == 'transcribe':
                res.append(tools.transcribe(seq))
            elif operation == 'reverse':
                res.append(tools.reverse(seq))
            elif operation == 'complement':
                res.append(tools.complement(seq))
            elif operation == 'reverse_complement':
                res.append(tools.reverse_complement(seq))
    return res[0] if len(res) == 1 else res


def filter_fastq(input_fastq: str,
                 output_fastq: str,
                 gc_bounds: Union[tuple[float, float], float] = (0, 100),
                 length_bounds: Union[tuple[int, int], int] = (0, 2 ** 32),
                 quality_threshold: float = 0) -> None:
    """
    This utility filter a file of FASTQ sequences
    based on GC content, sequence length, and
    quality threshold.
    Criteria for filtering:
    1. GC-content: sequences must fall within the
    specified GC range.
    2. Length: sequences must fall within the specified
    length range.
    3. Quality threshold: sequences must fall within
    the specified quality.
    Args:
    :param input_fastq: str:
    The input FASTQ file path containing
    the sequences to be filtered.
    :param output_fastq: str:
    The output FASTQ file path where filtered
    sequences will be saved.
    :param gc_bounds: tuple[float, float] | float = (0, 100):
    A tuple containing the lower and upper bound of the
    GC content percentage.
    If a single value (float) is provided, it will be used
    as the upper bound and 0 as the lower bound.
    Default value is (0, 100).
    :param length_bounds: tuple[int, int] | int = (0, 2**32):
    A tuple containing the minimum and maximum length
    of the sequence.
    If a single value (int) is provided, it will
    be used as the maximum length and 0 as
    the minimum length. Default value is (0, 2**32).
    :param quality_threshold: float = 0:
    Contain the minimum average quality score of the sequence.
    Default value is 0, which means no quality score is considered.
    :return: None: The function writes the filtered sequences
    directly to the specified output FASTQ file.
    Raises:
    FileExistsError: If the output file already exists,
    a message will be printed, and the function will exit.
    """
    if isinstance(gc_bounds, float):
        gc_bounds = (0.0, gc_bounds)
    elif isinstance(gc_bounds, tuple):
        if len(gc_bounds) == 1:
            gc_bounds = (0.0, gc_bounds[0])
    else:
        return
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    elif isinstance(length_bounds, tuple):
        if len(length_bounds) != 2:
            return
    else:
        return
    if os.path.exists(output_fastq):
        print(f"Error: file '{output_fastq}' already exists.")
        return
    for name, sequence, comment, quality in module_filter.read_fastq(input_fastq):
        gc_content = module_filter.calculate_gc_content(sequence)
        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        sequence_length = module_filter.length_calculation(sequence)
        if not (length_bounds[0] <= sequence_length <= length_bounds[1]):
            continue
        average_quality = module_filter.quality_calculation(quality)
        if average_quality < quality_threshold:
            continue
        module_filter.write_fastq(output_fastq, name, sequence, comment, quality)
