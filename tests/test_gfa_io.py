import io
import os
import pathlib
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock

import pytest
from Bio import SeqIO

from gfa_io import (
    GFA,
    BaseHelper,
    Containment,
    Link,
    MissingFieldsError,
    OptionalField,
    Path,
    Segment,
    find_isolated_segments,
    update_segments,
)

GFA_IN_1 = "test_data/MT.gfa"
GFA_IN_2 = "test_data/sample.gfa"
FASTA_IN = "test_data/MT.fasta"


def test_gfa_read() -> None:
    gfa = GFA(GFA_IN_1)
    assert len(gfa.segments) == 8
    assert len(gfa.links) == 11
    assert len(gfa.paths) == 0
    assert len(gfa.containments) == 0


def test_gfa_read_no_paths() -> None:
    gfa = GFA(GFA_IN_2)
    assert len(gfa.segments) == 6
    assert len(gfa.links) == 4
    assert len(gfa.paths) == 0
    assert len(gfa.containments) == 0


def test_to_networkx() -> None:
    gfa = GFA(GFA_IN_1)
    graph = gfa.to_networkx()
    assert graph.number_of_nodes() == len(gfa.segments)
    assert graph.number_of_edges() == len(gfa.links)


def test_to_fasta() -> None:
    gfa = GFA(GFA_IN_1)
    with TemporaryDirectory() as tmpdir:
        fasta_out = os.path.join(tmpdir, "test_out.fasta")
        with open(fasta_out, "w") as f:
            gfa.to_fasta(f)
        assert os.path.exists(fasta_out)
        seqs = [s for s in SeqIO.parse(fasta_out, "fasta")]
        assert len(list(seqs)) == len(gfa.segments)
        for seq in seqs:
            assert str(seq.seq) == gfa.segments[seq.id].seq


def test_find_isolated_segments() -> None:
    isolated = find_isolated_segments(GFA_IN_1)
    assert isinstance(isolated, list)


def test_update_segments() -> None:
    with TemporaryDirectory() as tmpdir:
        gfa_out = os.path.join(tmpdir, "test_out.gfa")
        with open(gfa_out, "w") as f:
            update_segments(GFA_IN_1, FASTA_IN, f)
        assert os.path.exists(gfa_out)
        

# Test for a successful update
def test_update_segments_success(tmp_path: pathlib.Path) -> None:
    """Tests a successful segment update."""
    gfa_content = "H\tVN:Z:1.0\nS\t1\tAAAAA\tSN:Z:fake\tSO:i:0\tSR:i:0\n"
    fasta_content = ">1\nCCCCC\n"
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(gfa_content)
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)
    with open(tmp_path / "out.gfa", 'w+t') as output_stream:
        update_segments(str(gfa_file), str(fasta_file), output_stream)
        output_stream.seek(0)
        result = [s.strip() for s in output_stream.readlines()]
        # The sequence and length should be updated
        assert "S\t1\tCCCCC\tSN:Z:fake\tSO:i:0\tSR:i:0" in result


# Test case for a corrupt GFA file
def test_update_segments_corrupt_gfa(tmp_path: pathlib.Path) -> None:
    """Tests with a corrupt GFA file."""
    gfa_file = tmp_path / "corrupt.gfa"
    # This segment line is missing fields
    gfa_file.write_text("H\tVN:Z:1.0\nS\t1\n")
    fasta_file = tmp_path / "dummy.fasta"
    fasta_file.touch()
    with pytest.raises(MissingFieldsError), open(tmp_path / "out.gfa", 'w+t') as output_stream:
        update_segments(str(gfa_file), str(fasta_file), output_stream)


# Test case for a non-existent GFA file
def test_update_segments_gfa_does_not_exist(tmp_path: pathlib.Path) -> None:
    """Tests with a non-existent GFA file."""
    with pytest.raises(FileNotFoundError), open(tmp_path / "out.gfa", 'w+t') as output_stream:
        update_segments("non_existent.gfa", "dummy.fasta", output_stream)


# Test case for a GFA file with no segments
def test_update_segments_no_segments_in_gfa(tmp_path: pathlib.Path) -> None:
    """Tests with a GFA file that contains no segments."""
    gfa_file = tmp_path / "no_segments.gfa"
    gfa_file.write_text("H\tVN:Z:1.0\n")
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(">1\nACGT\n")
    with (pytest.raises(KeyError, match="Fasta sequence 1 was not found as a segment"),
          open(tmp_path / "out.gfa", 'w+t') as output_stream):
        update_segments(str(gfa_file), str(fasta_file), output_stream)


# Test case for a corrupt FASTA file
def test_update_segments_corrupt_fasta(tmp_path: pathlib.Path) -> None:
    """Tests with a corrupt FASTA file."""
    gfa_content = "H\tVN:Z:1.0\nS\t1\tAAAAA\tSN:Z:fake\tSO:i:0\tSR:i:0\n"
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(gfa_content)
    fasta_file = tmp_path / "corrupt.fasta"
    # Incomplete FASTA record
    fasta_file.write_text(">1\nACGT\n>")
    # SeqIO.parse raises a ValueError for corrupt FASTA files
    with pytest.raises(KeyError), open(tmp_path / "out.gfa", 'w+t') as output_stream:
        update_segments(str(gfa_file), str(fasta_file), output_stream)


# Test case for a non-existent FASTA file
def test_update_segments_fasta_does_not_exist(tmp_path: pathlib.Path) -> None:
    """Tests with a non-existent FASTA file."""
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text("H\tVN:Z:1.0\nS\t1\t*\n")
    with pytest.raises(FileNotFoundError), open(tmp_path / "out.gfa", 'w+t') as output_stream:
        update_segments(str(gfa_file), "non_existent.fasta", output_stream)


# Test case for a FASTA file with no records
def test_update_segments_no_records_in_fasta(tmp_path: pathlib.Path) -> None:
    """Tests with a FASTA file that contains no sequences."""
    gfa_content = "H\tVN:Z:1.0\nS\t1\tAAAAA\tSN:Z:fake\tSO:i:0\tSR:i:0\n"
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(gfa_content)
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.touch()

    with open(tmp_path / "out.gfa", 'w+t') as output_stream:
        update_segments(str(gfa_file), str(fasta_file), output_stream)
        output_stream.seek(0)
        # The output should be the same as the input GFA
        assert output_stream.read().strip() == gfa_content.strip()


# Test case when a segment from the FASTA is not in the GFA
def test_update_segments_segment_not_found(tmp_path: pathlib.Path) -> None:
    """Tests when a sequence in the FASTA file is not found in the GFA."""
    gfa_content = "H\tVN:Z:1.0\nS\t1\t*\n"
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(gfa_content)
    fasta_file = tmp_path / "other.fasta"
    fasta_file.write_text(">2\nGGGGG\n")
    with (pytest.raises(KeyError, match="Fasta sequence 2 was not found as a segment"),
        open(tmp_path / "out.gfa", 'w+t') as output_stream):
        update_segments(str(gfa_file), str(fasta_file), output_stream)


# Test case for writing failure
def test_update_segments_write_fails(tmp_path: pathlib.Path) -> None:
    """Tests the behavior when writing the output fails."""
    gfa_content = "H\tVN:Z:1.0\nS\t1\t*\n"
    fasta_content = ">1\nCCCCC\n"
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(gfa_content)
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    mock_stream = MagicMock(spec=io.TextIOWrapper)
    mock_stream.write.side_effect = IOError("Cannot write to stream")

    with pytest.raises(IOError, match="Cannot write to stream"):
        update_segments(str(gfa_file), str(fasta_file), mock_stream)
        

# A concrete implementation of BaseHelper for testing purposes
class ConcreteHelper(BaseHelper):
    """A minimal concrete implementation of the abstract BaseHelper class."""

    def __init__(self, min_fields: int) -> None:
        super().__init__(min_fields)

    def _id(self) -> str:
        """A dummy implementation for the abstract _id method."""
        return "test_id"


def test_read_optionals_duplicate_tags() -> None:
    """
    Tests that read_optionals raises a KeyError when duplicate tags are present.

    This test simulates a record with two optional fields that share the same
    tag ('LN'), which should trigger a name collision.
    """
    # We set min_fields to 2. Any fields beyond the second are considered optional.
    helper = ConcreteHelper(min_fields=2)
    fields = [
        "S",
        "segment1",
        "LN:i:100",  # First optional field with tag 'LN'
        "RC:i:500",
        "KC:i:12345",
        "LN:Z:another_val",  # A second field with the same 'LN' tag
    ]

    # We expect a KeyError because the tag 'LN' is duplicated.
    with pytest.raises(KeyError) as excinfo:
        helper.read_optionals(fields)

    # Assert that the correct exception was raised with a descriptive message.
    assert "optional fields with duplicate tag LN" in str(excinfo.value)


@pytest.mark.parametrize(
    "field_str",
    [
        "1a:i:123",  # Starts with a number
        "a:i:123",  # Single character tag
        "abc:i:123",  # More than two characters
        "a!:i:123",  # Contains a special character
    ],
)
def test_optional_field_unrecognised_tag_raises_error(field_str: str) -> None:
    """
    Tests that initializing OptionalField with an invalid tag raises an IOError.

    The GFA specification requires tags to be two characters long, starting
    with a letter, followed by a letter or a digit.
    """
    with pytest.raises(IOError, match="invalid field tag"):
        OptionalField(field_str)


def test_optional_field_unknown_type_raises_error() -> None:
    """
    Tests that an unknown data type identifier raises an IOError.
    """
    with pytest.raises(ValueError, match="unknown field type"):
        # 'X' is not a valid GFA optional field type
        OptionalField("RC:X:some_value")


@pytest.mark.parametrize(
    "field_str, expected_value, expected_type",
    [
        # Type 'A': Printable character
        ("T1:A:c", "c", str),
        # Type 'Z': Printable string
        ("T2:Z:a printable string", "a printable string", str),
        # Type 'J': JSON format string
        ('T3:J:{"key": "value", "num": 1}', '{"key": "value", "num": 1}', str),
        # Type 'i': Integer
        ("T4:i:12345", 12345, int),
        ("T5:i:-54321", -54321, int),
        # Type 'f': Float
        ("T6:f:3.14159", 3.14159, float),
        # Type 'H': Hex-encoded byte array
        ("T7:H:414243", bytearray(b"ABC"), bytearray),
        # Type 'B': Numeric array (integer)
        ("T8:B:c,1,-2,3,-4", [1, -2, 3, -4], list),
        ("T9:B:f,1.1,2.2,-3.3", [1.1, 2.2, -3.3], list),
    ],
)
def test_optional_field_handles_known_types_correctly(field_str: str,
                                                      expected_value,
                                                      expected_type) -> None:
    """
    Tests that OptionalField correctly parses and returns values for all
    supported data types.
    """
    opt_field = OptionalField(field_str)
    val = opt_field.get()

    assert val == expected_value
    assert isinstance(val, expected_type)


def test_optional_field_invalid_numeric_array_type() -> None:
    """
    Tests that an unknown numeric array subtype in a 'B' field raises an IOError.
    """
    # Create a field with an invalid numeric array type 'x'
    opt_field = OptionalField("NV:B:x,1,2,3")

    with pytest.raises(ValueError, match="unknown number array type"):
        opt_field.get()


def test_segment_with_sequence_and_optionals() -> None:
    """
    Tests parsing a standard segment line with a sequence and optional fields.
    Sequence data should be loaded.
    """
    # GFA line format: S <Name> <Sequence> [Optional Fields...]
    line = "S\tseg1\tGATTACA\tLN:i:7\tRC:i:100"
    segment = Segment(line, skip_seq=False)

    assert segment.name == "seg1"
    assert segment.seq == "GATTACA"
    assert segment.length == 7
    assert len(segment.optionals) == 2
    assert "LN" in segment.optionals
    assert segment.optionals["LN"].get() == 7
    assert "RC" in segment.optionals
    assert segment.optionals["RC"].get() == 100


def test_segment_with_sequence_and_skip_seq_true() -> None:
    """
    Tests parsing a segment with a sequence when skip_seq is True.
    The sequence should not be loaded into memory, but its length should be recorded.
    """
    line = "S\tseg2\tACGT"
    segment = Segment(line, skip_seq=True)

    assert segment.name == "seg2"
    assert segment.seq is None  # Sequence should not be stored
    assert segment.length == 4  # Length should be calculated
    assert not segment.optionals


def test_segment_with_dummy_sequence() -> None:
    """
    Tests parsing a segment where the sequence is a '*' (dummy).
    The sequence should be None and the length should be 0.
    """
    line = "S\tseg3\t*\tKC:i:500"
    # The value of skip_seq should not matter when the sequence is a dummy '*'
    segment = Segment(line, skip_seq=False)

    assert segment.name == "seg3"
    assert segment.seq is None
    assert segment.length == 0
    assert "KC" in segment.optionals
    assert segment.optionals["KC"].get() == 500


def test_segment_with_dummy_sequence_and_skip_seq_true() -> None:
    """
    Tests parsing a segment with a dummy sequence when skip_seq is True.
    The result should be the same as when skip_seq is False.
    """
    line = "S\tseg4\t*"
    segment = Segment(line, skip_seq=True)

    assert segment.name == "seg4"
    assert segment.seq is None
    assert segment.length == 0
    assert not segment.optionals


def test_segment_with_insufficient_fields_raises_error() -> None:
    """
    Tests that a segment line with fewer than the required 3 fields
    raises a MissingFieldsError.
    """
    # This line is missing the sequence field
    line = "S\tseg5"
    with pytest.raises(MissingFieldsError):
        Segment(line, skip_seq=False)


def test_link_parsing_with_overlap_and_optionals() -> None:
    """
    Tests the parsing of a standard link line that includes a CIGAR
    string for overlap and additional optional fields.
    """
    # GFA line format: L <From> <FromOr> <To> <ToOr> <Overlap> [Tags...]
    line = "L\tseg1\t+\tseg2\t-\t8M\tRC:i:50\tNM:i:1"
    link = Link(line)

    assert link.src == "seg1"
    assert link.src_orient == "+"
    assert link.dest == "seg2"
    assert link.dest_orient == "-"
    assert link.overlap == "8M"
    assert len(link.optionals) == 2
    assert "RC" in link.optionals
    assert link.optionals["RC"].get() == 50
    assert "NM" in link.optionals
    assert link.optionals["NM"].get() == 1


def test_link_parsing_with_dummy_overlap() -> None:
    """
    Tests parsing a link line where the overlap is a '*' (dummy value).
    The overlap attribute should be None, and optional fields should be parsed.
    """
    line = "L\tseg3\t-\tseg4\t+\t*\tKC:i:1024"
    link = Link(line)

    assert link.src == "seg3"
    assert link.src_orient == "-"
    assert link.dest == "seg4"
    assert link.dest_orient == "+"
    assert link.overlap is None  # Overlap should be None for '*'
    assert len(link.optionals) == 1
    assert "KC" in link.optionals
    assert link.optionals["KC"].get() == 1024


def test_link_parsing_without_optionals() -> None:
    """
    Tests an edge case where the link line has no optional fields.
    """
    line = "L\tseg5\t+\tseg6\t+\t0M"
    link = Link(line)

    assert link.src == "seg5"
    assert link.dest == "seg6"
    assert link.overlap == "0M"
    assert not link.optionals  # No optional fields should be present


def test_link_with_insufficient_fields_raises_error() -> None:
    """
    Tests that a link line with fewer than the required 6 fields raises
    a MissingFieldsError.
    """
    # This line is missing the overlap and orientation fields
    line = "L\tseg7\tseg8"
    with pytest.raises(MissingFieldsError):
        Link(line)


# Tests for the Containment class
def test_containment_parsing_with_overlap_and_optionals() -> None:
    """
    Tests parsing a standard containment line with an overlap and optional fields.
    """
    # GFA line format: C <Container> <ContOr> <Contained> <ContdOr> <Pos> <Overlap> [Tags..]
    line = "C\tseg1\t+\tseg2\t-\t100\t8M4I\tNM:i:4"
    containment = Containment(line)

    assert containment.container == "seg1"
    assert containment.container_orient == "+"
    assert containment.contained == "seg2"
    assert containment.contained_orient == "-"
    assert containment.pos == 100
    assert containment.overlap == "8M4I"
    assert "NM" in containment.optionals
    assert containment.optionals["NM"].get() == 4


def test_containment_parsing_with_dummy_overlap() -> None:
    """
    Tests parsing a containment line with a '*' dummy overlap.
    The overlap attribute should be None.
    """
    line = "C\tseg3\t-\tseg4\t+\t250\t*"
    containment = Containment(line)

    assert containment.container == "seg3"
    assert containment.contained == "seg4"
    assert containment.pos == 250
    assert containment.overlap is None
    assert not containment.optionals


def test_containment_with_invalid_pos_raises_error() -> None:
    """
    Tests that a containment line with a non-integer position raises a ValueError.
    """
    line = "C\tseg5\t+\tseg6\t+\tABC\t*"
    with pytest.raises(ValueError):
        Containment(line)


def test_containment_with_insufficient_fields_raises_error() -> None:
    """
    Tests that a containment line with fewer than the required 7 fields
    raises a MissingFieldsError.
    """
    line = "C\tseg7\t+\tseg8\t-"
    with pytest.raises(MissingFieldsError):
        Containment(line)


# Tests for the Path class
def test_path_parsing_with_overlaps() -> None:
    """
    Tests parsing a standard path line with segment names and overlaps.
    """
    # GFA line format: P <PathName> <SegmentNames> <Overlaps>
    line = "P\tpath1\tseg1+,seg2-,seg3+\t10M,5M"
    path = Path(line)

    assert path.name == "path1"
    assert path.segment_names == [("seg1", "+"), ("seg2", "-"), ("seg3", "+")]
    assert path.overlaps == ["10M", "5M"]


def test_path_parsing_with_dummy_overlaps() -> None:
    """
    Tests parsing a path line with a '*' dummy overlap field.
    The overlaps list should be empty.
    """
    line = "P\tpath2\tseg4+,seg5+\t*"
    path = Path(line)

    assert path.name == "path2"
    assert path.segment_names == [("seg4", "+"), ("seg5", "+")]
    assert not path.overlaps


def test_path_with_single_segment() -> None:
    """
    Tests a path with only one segment. The overlaps list should be empty
    as per the GFA specification (usually indicated by '*').
    """
    line = "P\tpath3\tseg6-\t*"
    path = Path(line)

    assert path.name == "path3"
    assert path.segment_names == [("seg6", "-")]
    assert not path.overlaps


def test_path_with_insufficient_fields_raises_error() -> None:
    """
    Tests that a path line with fewer than the required 4 fields
    raises a MissingFieldsError.
    """
    line = "P\tpath4"
    with pytest.raises(MissingFieldsError):
        Path(line)
