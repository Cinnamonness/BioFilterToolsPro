import os
import sys
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from BioFilterToolsPro import filter_fastq


@pytest.fixture
def test_fastq_file(tmp_path):
    """
    Create a temporary FASTQ file with three synthetic reads
    for testing filtering functions.
    """
    input_path = tmp_path / "test_input.fastq"
    records = [
        SeqRecord(
            Seq("ATGC" * 10),
            id="read1",
            letter_annotations={"phred_quality": [20] * 40},
        ),
        SeqRecord(
            Seq("AAAA" * 10),
            id="read2",
            letter_annotations={"phred_quality": [40] * 40},
        ),
        SeqRecord(
            Seq("GCGC" * 10),
            id="read3",
            letter_annotations={"phred_quality": [10] * 40},
        ),
    ]
    with open(input_path, "w") as f:
        SeqIO.write(records, f, "fastq")
    return input_path


class TestBasicFunctionality:
    """
    Tests for basic filtering functionality
    such as length, quality and GC content.
    """

    def test_length_filtering(self, test_fastq_file, tmp_path):
        """
        Test length filtering.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(str(test_fastq_file), str(output), length_bounds=(30, 50))
        records = list(SeqIO.parse(output, "fastq"))
        assert len(records) == 3

    def test_filter_quality(self, test_fastq_file, tmp_path):
        """
        Test quality filtering.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(str(test_fastq_file), str(output), quality_threshold=30)
        records = list(SeqIO.parse(output, "fastq"))
        assert len(records) == 1
        assert records[0].id == "read2"

    def test_gc_filtering(self, test_fastq_file, tmp_path):
        """
        Test GC content filtering.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(str(test_fastq_file), str(output), gc_bounds=(30, 70))
        records = list(SeqIO.parse(output, "fastq"))
        assert len(records) == 1
        assert records[0].id == "read1"

    def test_file_creation(self, test_fastq_file, tmp_path):
        """
        Test that the output file is created and is not empty.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(str(test_fastq_file), str(output))
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestError:
    """
    Tests for proper error handling
    when output file already exists.
    """

    def test_file_exists_error(self, test_fastq_file, tmp_path):
        """
        Test that FileExistsError is raised if
        output file already exists.
        """
        output = tmp_path / "output.fastq"
        output.touch()
        with pytest.raises(FileExistsError):
            filter_fastq(str(test_fastq_file), str(output))


class TestParameters:
    """
    Tests for advanced filter parameter behaviors such as
    single GC and single length
    """

    def test_single_gc_content(self, test_fastq_file, tmp_path):
        """
        Test that a single GC bound value
        (interpreted as lower bound) works.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(str(test_fastq_file), str(output), gc_bounds=50.0)
        records = list(SeqIO.parse(output, "fastq"))
        assert len(records) == 2
        assert records[0].id == "read1"

    def test_single_length(self, test_fastq_file, tmp_path):
        """
        Test that a single length bound value
        (interpreted as lower bound) works.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(str(test_fastq_file), str(output), length_bounds=50)
        records = list(SeqIO.parse(output, "fastq"))
        assert len(records) == 3

    def test_combine_parameters(self, test_fastq_file, tmp_path):
        """
        Test that a combine filter gc_content, length bound,
        quality works.
        """
        output = tmp_path / "output.fastq"
        filter_fastq(
            str(test_fastq_file),
            str(output),
            length_bounds=50,
            gc_bounds=50.0,
            quality_threshold=30,
        )
        records = list(SeqIO.parse(output, "fastq"))
        assert len(records) == 1
