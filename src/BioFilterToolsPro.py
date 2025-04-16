from abc import ABC, abstractmethod
from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import argparse
import logging


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
        if not hasattr(self, "complement_map"):
            raise NotImplementedError("Complement map must be defined in subclass.")
        complement_sequence = self.sequence.translate(
            str.maketrans(*self.complement_map)
        )
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

    valid_alphabet = {"A", "T", "G", "C"}
    complement_map = ("ATGC", "TACG")

    def transcribe(self) -> str:
        """
        Transcribes a DNA sequence to
        RNA.

        :returns: str: transcribed RNA sequence.
        """
        transcribed_sequence = self.sequence.replace("T", "U")
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

    valid_alphabet = {"A", "U", "G", "C"}
    complement_map = ("AUGC", "UACG")


class AminoAcidSequence(BiologicalSequence):
    """
    This class implements the
    BiologicalSequence interface.
    """

    valid_alphabet = {
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    }

    def __init__(self, sequence: str):
        """ """
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
            "non-polar": 0,
            "polar uncharged": 0,
            "polar negatively charged": 0,
            "polar positively charged": 0,
        }

        classification_map = {
            "non-polar": ["G", "A", "L", "I", "V", "P", "M", "W", "F"],
            "polar uncharged": ["N", "Q", "S", "T", "Y", "C"],
            "polar negatively charged": ["D", "E"],
            "polar positively charged": ["R", "H", "K"],
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


logging.basicConfig(
    filename="logs.log",
    filemode="w",
    format="{levelname} | {asctime} --> {module} {funcName} --> {message}",
    datefmt="%Y-%m-%d %H:%M:%S",
    style="{",
    level=logging.DEBUG,
    force=True,
    encoding="utf-8",
    errors="ignore",
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[tuple[float, float], float] = (0, 100),
    length_bounds: Union[tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters a FASTQ file based on GC content, sequence length,
    and quality threshold.
    Logs all significant operations and errors to 'logs.log' file.

    Filtering Criteria:
    1. GC-content: Sequences must fall within specified GC percentage range
    2. Length: Sequences must be within specified length bounds
    3. Quality: Sequences must meet minimum average quality score

    Logging Behavior:
    - INFO level: Successful operations (parameter conversions,
      directory creation, completion stats)
    - ERROR level: Invalid parameters, file existence conflicts,
      processing errors
    - DEBUG level: Detailed sequence processing (if enabled)

    Args:
        input_fastq: Path to input FASTQ file
        output_fastq: Path for output FASTQ file (will raise error if exists)
        gc_bounds: GC content bounds as:
                  - float (upper bound, lower=0)
                  - tuple[lower, upper]
                  Default: (0, 100) - no filtering
        length_bounds: Length bounds as:
                      - int (upper bound, lower=0)
                      - tuple[min, max]
                      Default: (0, 2^32) - no filtering
        quality_threshold: Minimum average Phred quality score
                          Default: 0 - no quality filtering

    Returns:
        None: Writes filtered sequences to output file

    Raises:
        ValueError: For invalid gc_bounds or length_bounds format
        FileExistsError: If output file already exists

    Examples:
        Successful conversion logged:
        INFO: GC bounds converted to tuple (0, 80)
        INFO: Length bounds converted to tuple (0, 500)
        INFO: Filtering completed. Total: 100, passed: 75 (75.0%)

        Error cases logged:
        ERROR: Invalid gc_bounds. Must be tuple or float
        ERROR: Output file 'output.fastq' already exists
        ERROR: Error during filtering: [detailed traceback]
    """
    try:
        if isinstance(gc_bounds, float):
            gc_bounds = (0.0, gc_bounds)
            logging.info(f"GC bounds converted to tuple (0, {gc_bounds[1]})")
        elif isinstance(gc_bounds, tuple) and len(gc_bounds) == 1:
            gc_bounds = (0.0, gc_bounds[0])
            logging.info(f"GC bounds converted to tuple (0, {gc_bounds[1]})")
        elif not isinstance(gc_bounds, tuple) or len(gc_bounds) != 2:
            error_message = (
                "Invalid gc_bounds. It must be a tuple of two floats or a single float."
            )
            logging.error(error_message)
            raise ValueError(error_message)

        if isinstance(length_bounds, int):
            length_bounds = (0, length_bounds)
            logging.info(f"Length bounds converted to tuple (0, {length_bounds[1]})")
        elif isinstance(length_bounds, tuple) and len(length_bounds) != 2:
            error_message = "Invalid length_bounds. It must be a tuple of two integers."
            logging.error(error_message)
            raise ValueError(error_message)

        if os.path.exists(output_fastq):
            error_message = f"Output file '{output_fastq}' already exists."
            logging.error(error_message)
            raise FileExistsError(error_message)

        if not output_fastq:
            filtered_dir = os.path.join(os.getcwd(), "filtered")
            if not os.path.exists(filtered_dir):
                os.makedirs(filtered_dir)
                logging.info(f"Created directory: {filtered_dir}")
            output_fastq = os.path.join(filtered_dir, "filtered_sequences.fastq")

        total_sequences = 0
        passed_sequences = 0

        with open(output_fastq, "w") as output:
            for record in SeqIO.parse(input_fastq, "fastq"):
                total_sequences += 1
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
                passed_sequences += 1

        logging.info(
            f"Filtering completed. Total sequences: {total_sequences}, "
            f"passed sequences: {passed_sequences} ({passed_sequences/total_sequences:.1%})"
        )

    except Exception as e:
        logging.error(f"Error during filtering: {str(e)}", exc_info=True)
        raise


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments for FASTQ filtering.

    Parsed command line arguments with atributes:
        - input: str - input FASTQ file path.
        - output: str - output FASTQ file path.
        - length: List[int] - length bounds.
        - gc_content: List[float] - GC-content bounds.
        - quality: float - quality threshold.
    """
    parser = argparse.ArgumentParser(
        description="Filster FASTQ based on lenght, GC content and quality"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="This is the input FASTQ file path"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="This is the output FASTQ file path"
    )
    parser.add_argument(
        "-len",
        "--length",
        type=int,
        nargs="*",
        default=(0, 2**32),
        help="Length bounds. Either one value (upper) or two values (lower upper)",
    )
    parser.add_argument(
        "-gc",
        "--gc_content",
        type=float,
        nargs="*",
        default=(0, 100),
        help="GC-content bounds. Either one value (upper) or two values (lower upper)",
    )
    parser.add_argument(
        "-q",
        "--quality",
        type=float,
        default=0,
        help="Minimum average quality threshold",
    )
    return parser.parse_args()


def process_bounds(
    args_value: list, default_bounds: Union[tuple[float, float], tuple[int, int]]
) -> Union[tuple[float, float], tuple[int, int]]:
    """
    Process bounds arguments (gc_bounds and
    length_bounds).

    Args:
        - args_value: List of bound values from
        command line.
        - default_bounds: Default bounds to return
        if args_value is invalid.

    Returns:
        Tuple of (lower, upper) bounds processed
        from input arguments.
    """
    if len(args_value) == 1:
        result = (0, args_value[0])
        logging.info(f"Bounds converted to tuple {result}")
        return result
    elif len(args_value) == 2:
        return tuple(args_value)
    else:
        return default_bounds


def main():
    args = parse_args()
    gc_bounds = process_bounds(args.gc_content, (0, 100))
    length_bounds = process_bounds(args.length, (0, 2**32))

    filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality,
    )


if __name__ == "__main__":
    main()
