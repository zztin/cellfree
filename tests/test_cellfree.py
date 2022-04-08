from click.testing import CliRunner

import cellfree.cmd.bam as bam
import cellfree.cmd.motifs as motifs


def test_bam():
    runner = CliRunner()
    result = runner.invoke(bam.main)
    assert result.output == "bam command is called\n"


def test_motifs():
    runner = CliRunner()
    result = runner.invoke(motifs.main)
    assert result.output == "motif command is called\n"
