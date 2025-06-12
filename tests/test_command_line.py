import os
import subprocess
from tempfile import TemporaryDirectory

import pytest

GFA_IN = "test_data/MT.gfa"
FASTA_IN = "test_data/MT.fasta"


def test_update_segments():
    with TemporaryDirectory() as tmpdir:
        gfa_out = os.path.join(tmpdir, "test_out.gfa")
        subprocess.run(
            [
                "gfa_utils",
                "update-segments",
                GFA_IN,
                FASTA_IN,
                gfa_out,
            ],
            check=True,
        )
        assert os.path.exists(gfa_out)


def test_isolates():
    with TemporaryDirectory() as tmpdir:
        csv_out = os.path.join(tmpdir, "test_out.csv")
        subprocess.run(
            ["gfa_utils", "isolates", GFA_IN, csv_out],
            check=True,
        )
        assert os.path.exists(csv_out)


def test_dump_segments():
    with TemporaryDirectory() as tmpdir:
        fasta_out = os.path.join(tmpdir, "test_out.fasta")
        subprocess.run(
            [
                "gfa_utils",
                "dump-segments",
                GFA_IN,
                fasta_out,
            ],
            check=True,
        )
        assert os.path.exists(fasta_out)


@pytest.mark.parametrize(
    "format", ["graphml", "gml", "edgelist"]
)
def test_convert(format):
    with TemporaryDirectory() as tmpdir:
        graph_out = os.path.join(tmpdir, f"test_out.{format}")
        subprocess.run(
            [
                "gfa_utils",
                "convert",
                "--format",
                format,
                GFA_IN,
                graph_out,
            ],
            check=True,
        )
        assert os.path.exists(graph_out)