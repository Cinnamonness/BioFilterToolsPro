from abc import ABC, abstractmethod
from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os


class BiologicalSequence(ABC):
    """
    Abstract class for biological sequences
    that defines the following interface:

    - Support for the len function;
    - Ability to access elements by index
    and slice the sequence;
    - Pretty-printing the sequence and the
    ability to convert it to a string;
    - Ability to check the validity of
    the sequence alphabet.
    """
    @abstractmethod
    def __len__(self) -> int:
        """
        Returns the length of the
        biological sequence.

        :returns: int: length of sequence.
        """
        pass

    @abstractmethod
    def __getitem__(self, index: int) -> str:
        """
        Supports indexing and slicing
        the sequence.

        :param index: int or slice: index
        or slice of the sequence.
        :returns: str: sequence at the index
        or slice.
        """
        pass

    @abstractmethod
    def __str__(self) -> str:
        """
        Returns a string of sequence.

        :returns: str: string of sequence.
        """
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """
        Returns a pretty format of sequence.

        :returns: str: pretty string of sequence.
        """
        pass

    @abstractmethod
    def is_valid_alphabet(self) -> bool:
        """
        Check that sequence contains valid characters.

        :returns: bool: True if all characters are valid,
        False if not.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    This class implements the BiologicalSequence
    interface.

    Attributes:
    sequence: str: the biological sequence.
    """
    valid_alphabet: set
    complement_map: tuple

    def __init__(self, sequence: str):
        """
        Initialize a NucleicAcidSequence.

        :param sequence: str: input string.
        :raises ValueError: if the sequence
        contains invalid characters.
        """
        self.sequence = sequence.upper()
        if not self.is_valid_alphabet():
            raise ValueError("Invalid alphabet for biological sequence.")

    def __len__(self) -> int:
        """
        Returns the length of the sequence.

        :returns: int: length of sequence.
        """
        return len(self.sequence)

    def __getitem__(self, i) -> str:
        """
        Supports indexing and slicing
        the sequence.

        :param index: int or slice: index
        or slice of the sequence.
        :returns: str: sequence at the index
        or slice.
        """
        if isinstance(i, slice):
            return self.__class__(self.sequence[i])
        else:
            return self.sequence[i]

    def __str__(self) -> str:
        """
        Returns a string of sequence.

        :returns: str: string of sequence.
        """
        return f"Sequence: {self.sequence}"

    def __repr__(self) -> str:
        """
        Returns a pretty format of sequence.

        :returns: str: pretty string of sequence.
        """
        return f"{self.__class__.__name__}('{self.sequence}')"

    def is_valid_alphabet(self) -> bool:
        """
        Check that sequence contains valid characters.

        :returns: bool: True if all characters are valid,
        False if not.
        """
        return all(char in self.valid_alphabet for char in self.sequence)

    def complement(self) -> "NucleicAcidSequence":
        """
        Returns the complement of the sequence.

        :returns: NucleicAcidSequence: complement sequence.
        :raises NotImplementedError: if complement_map
        is not defined.
        """
        if not hasattr(self, 'complement_map'):
            raise NotImplementedError("Complement map must be defined in subclass.")
        complement_sequence = self.sequence.translate(str.maketrans(*self.complement_map))
        return self.__class__(complement_sequence)

    def reverse(self) -> "NucleicAcidSequence":
        """
        Returns the reversed sequence.

        :returns: NucleicAcidSequence: reversed sequence.
        """
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """
        Returns the reverse complement
        of the sequence.

        :returns: NucleicAcidSequence: reversed
        and complement sequence.
        """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    DNA sequence class that extends
    NucleicAcidSequence.

    Attributes:
    valid_alphabet: set: valid characters
    for DNA sequence ('A', 'T', 'G', 'C').
    complement_map: tuple: mapping for
    complement of DNA sequence.
    """
    valid_alphabet = {'A', 'T', 'G', 'C'}
    complement_map = ('ATGC', 'TACG')

    def transcribe(self) -> str:
        """
        Transcribes a DNA sequence to
        RNA.

        :returns: str: transcribed RNA sequence.
        """
        transcribed_sequence = self.sequence.replace('T', 'U')
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
    """
    RNA sequence class that extends
    NucleicAcidSequence.

    Attributes:
    valid_alphabet: set: valid characters
    for RNA sequence ('A', 'U', 'G', 'C').
    complement_map: tuple: mapping for
    complement of RNA sequence.
    """
    valid_alphabet = {'A', 'U', 'G', 'C'}
    complement_map = ('AUGC', 'UACG')


class AminoAcidSequence(BiologicalSequence):
    """
    This class implements the
    BiologicalSequence interface.
    """
    valid_alphabet = {'A', 'R', 'N', 'D', 'C',
                      'Q', 'E', 'G', 'H', 'I',
                      'L', 'K', 'M', 'F', 'P',
                      'S', 'T', 'W', 'Y', 'V'}

    def __init__(self, sequence: str):
        """
        """
        self.sequence = sequence.upper()
        if not self.is_valid_alphabet():
            raise ValueError("Invalid alphabet for biological sequence.")

    def __len__(self) -> int:
        """
        Returns the length of the sequence.

        :returns: int: length of sequence.
        """
        return len(self.sequence)

    def __getitem__(self, i) -> str:
        """
        Supports indexing and slicing
        the sequence.

        :param index: int or slice: index
        or slice of the sequence.
        :returns: str: sequence at the index
        or slice.
        """
        if isinstance(i, slice):
            return self.__class__(self.sequence[i])
        else:
            return self.sequence[i]

    def __str__(self) -> str:
        """
        Returns a string of sequence.

        :returns: str: string of sequence.
        """
        return f"Sequence: {self.sequence}"

    def __repr__(self) -> str:
        """
        Returns a pretty format of sequence.

        :returns: str: pretty string of sequence.
        """
        return f"{self.__class__.__name__}('{self.sequence}')"

    def is_valid_alphabet(self) -> bool:
        """
        Check that sequence contains valid characters.

        :returns: bool: True if all characters are valid,
        False if not.
        """
        return all(char in self.valid_alphabet for char in self.sequence)

    def classify_aminoacids(self) -> str:
        """
        Classifies amino acids into categories: non-polar,
        polar uncharged, polar negatively charged and
        polar positively charged.

        :returns: str: Count and percentage of each category
        of amino acids.
        """
        classification_result = {
            'non-polar': 0,
            'polar uncharged': 0,
            'polar negatively charged': 0,
            'polar positively charged': 0
        }

        classification_map = {
            'non-polar': ['G', 'A', 'L', 'I', 'V', 'P', 'M', 'W', 'F'],
            'polar uncharged': ['N', 'Q', 'S', 'T', 'Y', 'C'],
            'polar negatively charged': ['D', 'E'],
            'polar positively charged': ['R', 'H', 'K']
        }

        for aa in self.sequence:
            for category, amino_acids in classification_map.items():
                if aa in amino_acids:
                    classification_result[category] += 1
                    break

        total_count = len(self.sequence)
        result_str = "\n".join(
            [
                f"{category}: Count = {count}, Percentage = {count / total_count * 100:.2f}%"
                for category, count in classification_result.items()
            ]
        )

        return result_str


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

    :raises ValueError: If the provided `gc_bounds` or `length_bounds` are invalid.
    :raises FileExistsError: If the output file already exists.
    """
    if isinstance(gc_bounds, float):
        gc_bounds = (0.0, gc_bounds)
    elif isinstance(gc_bounds, tuple) and len(gc_bounds) == 1:
        gc_bounds = (0.0, gc_bounds[0])
    elif not isinstance(gc_bounds, tuple) or len(gc_bounds) != 2:
        raise ValueError("Invalid gc_bounds. It must be a tuple of two floats or a single float.")

    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    elif isinstance(length_bounds, tuple) and len(length_bounds) != 2:
        raise ValueError("Invalid length_bounds. It must be a tuple of two integers.")

    if os.path.exists(output_fastq):
        raise FileExistsError(f"Error: file '{output_fastq}' already exists.")

    if not output_fastq:
        filtered_dir = os.path.join(os.getcwd(), 'filtered')
        if not os.path.exists(filtered_dir):
            os.makedirs(filtered_dir)
        output_fastq = os.path.join(filtered_dir, 'filtered_sequences.fastq')

    with open(output_fastq, "w") as output:
        for record in SeqIO.parse(input_fastq, "fastq"):
            gc_content_percent = gc_fraction(record.seq) * 100
            if not (gc_bounds[0] <= gc_content_percent <= gc_bounds[1]):
                continue
            sequence_length = len(record.seq)
            if not (length_bounds[0] <= sequence_length <= length_bounds[1]):
                continue
            phred_quality = record.letter_annotations["phred_quality"]
            if (sum(phred_quality) / len(phred_quality)) < quality_threshold:
                continue

            SeqIO.write(record, output, "fastq")
