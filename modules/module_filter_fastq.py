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
