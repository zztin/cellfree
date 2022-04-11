from click.testing import CliRunner

import cellfree.cmd.prepare as prepare
import cellfree.cmd.end_motif as end_motif


def test_bam():
    runner = CliRunner()
    result = runner.invoke(prepare.prepare, ['--help'])
    assert result.output == "prepare subcommand is called, yay!n"


def test_motifs():
    runner = CliRunner()
    result = runner.invoke(end_motif.end_motif)
    assert result.output == "end motif."
