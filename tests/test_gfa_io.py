import os
from Bio import SeqIO
from tempfile import TemporaryDirectory

import pytest
from gfa_io import GFA, find_isolated_segments, update_segments

GFA_IN_1 = "test_data/MT.gfa"
GFA_IN_2 = "test_data/sample.gfa"
FASTA_IN = "test_data/MT.fasta"


def test_gfa_read():
    gfa = GFA(GFA_IN_1)
    assert len(gfa.segments) == 8
    assert len(gfa.links) == 11
    assert len(gfa.paths) == 0
    assert len(gfa.containments) == 0


def test_gfa_read_no_paths():
    gfa = GFA(GFA_IN_2)
    assert len(gfa.segments) == 6
    assert len(gfa.links) == 4
    assert len(gfa.paths) == 0
    assert len(gfa.containments) == 0


def test_to_networkx():
    gfa = GFA(GFA_IN_1)
    graph = gfa.to_networkx()
    assert graph.number_of_nodes() == len(gfa.segments)
    assert graph.number_of_edges() == len(gfa.links)


def test_to_fasta():
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


def test_find_isolated_segments():
    isolated = find_isolated_segments(GFA_IN_1)
    assert isinstance(isolated, list)


def test_update_segments():
    with TemporaryDirectory() as tmpdir:
        gfa_out = os.path.join(tmpdir, "test_out.gfa")
        with open(gfa_out, "w") as f:
            update_segments(GFA_IN_1, FASTA_IN, f)
        assert os.path.exists(gfa_out)