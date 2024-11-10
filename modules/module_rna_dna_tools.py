"""
This module contains procedures for working with RNA and DNA sequences.
"""


def molecule_type(sequence: str):
    """
    Function to determine the type of molecule.
    :param sequence: str.
    :return: 'DNA' or 'RNA'.
    """
    seq = set(sequence)
    if seq.issubset({'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}):
        return 'DNA'
    elif seq.issubset({'A', 'C', 'G', 'U', 'a', 'c', 'g', 'u'}):
        return 'RNA'
    else:
        return 'This is not a DNA or RNA sequence'


def transcribe(sequence: str) -> str:
    """
    Function to determine transcribed sequence.
    :param sequence: str.
    :return: str.
    """
    return sequence.replace('t', 'u').replace('T', 'U')


def reverse(sequence: str) -> str:
    """
    Function to reverse the sequence.
    :param sequence: str.
    :return: str.
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Function to determine complement sequence.
    :param sequence: str.
    :return: str.
    """
    dna_map = str.maketrans('ATGCatgc', 'TACGtacg')
    rna_map = str.maketrans('AUGCaugc', 'UACGuacg')
    type_of_molecule = molecule_type(sequence)
    if type_of_molecule == 'DNA':
        return sequence.translate(dna_map)
    elif type_of_molecule == 'RNA':
        return sequence.translate(rna_map)


def reverse_complement(sequence: str) -> str:
    """
    Function to determine reverse complement sequence.
    :param sequence: str.
    :return: str.
    """
    return reverse(complement(sequence))
