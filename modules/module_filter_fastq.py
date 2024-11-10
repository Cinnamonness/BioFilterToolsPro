import os
from typing import Generator, Tuple


def calculate_gc_content(sequence: str) -> float:
    """
    Function to calculate GC-content of sequence.
    :param sequence: str.
    :return: float: Percentage of GC sequence.
    """
    gc_content = (sequence.count('G') + sequence.count('C'))
    gc_content_percent = gc_content / len(sequence) * 100
    return gc_content_percent


def length_calculation(sequence: str) -> int:
    """
    Function to calculate length of sequence.
    :param sequence: str.
    :return: int: Length of sequence.
    """
    return len(sequence)


def quality_calculation(quality: str) -> float:
    """
    Function to calculate quality of sequence.
    :param quality: str.
    :return: float: Percentage of quality sequence.
    """
    return sum(ord(quality[char]) - 33 for char in range(len(quality))) / len(quality)


def read_fastq(input_fastq: str) \
        -> Generator[Tuple[str, str, str, str], None, None]:
    """
    Read a FASTQ file and yield the name,
    sequence, comment, and quality.
    :param input_fastq: str.
    Path to the FASTQ file to be read.
    :return: A generator that yields tuples containing:
     - name: str - The name of the sequence.
     - seq: str - The sequence.
     - comment: str - The comment.
     - quality: str - The quality of the sequence.
    """
    with open(input_fastq, 'r') as f:
        while True:
            name = f.readline().rstrip()
            if not name:
                break
            seq = f.readline().rstrip()
            comment = f.readline().rstrip()
            quality = f.readline().rstrip()
            yield name, seq, comment, quality


def write_fastq(output_fastq: str,
                sequence_name: str,
                sequence: str,
                comment: str,
                quality: str) -> None:
    """
    Write a sequence name, sequence, comment,
    and quality sequence to a FASTQ file.
    :param output_fastq: str.
    Path to the output FASTQ file.
    :param sequence_name: str.
    The name of the sequence (starts with '@').
    :param sequence: str.
    The sequence to be written.
    :param comment: str.
    The comment to be written.
    :param quality: str.
    The quality sequence to be written.
    :return: None.
    The function writes the sequence
    information to the FASTQ file.
    """
    output_dir = os.path.dirname(output_fastq)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_fastq, 'a') as f:
        f.write(f"{sequence_name}\n{sequence}\n{comment}\n{quality}\n")
