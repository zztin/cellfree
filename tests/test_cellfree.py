import inspect
import os

import pytest
from click.testing import CliRunner

import cellfree.cmd.end_motif as end_motif
import cellfree.cmd.prepare as prepare
from cellfree.cmd.cellfree import main


@pytest.mark.skip(reason="still work in progress")
def test_bamtagmultiome_benchmark():
    """FIXME: Only works with " tests$ pytest test_cellfree.py::test_bamtagmultiome_benchmark" but not "pytest"
    it is very likely the same root cause of the argparser workaround of cellfree.algorithm.bam.run_tagging
    """
    expected = """Started writing to /tmp/tagged.bam
Params: {'query_name_flagger': None, 'molecule_class': <class 'singlecellmultiomics.molecule.nlaIII.NlaIIIMolecule'>, 'fragment_class': <class 'singlecellmultiomics.fragment.nlaIII.NlaIIIFragment'>, 'molecule_class_args': {'umi_hamming_distance': 1, 'reference': None}, 'fragment_class_args': {'read_group_format': 0}, 'yield_invalid': True, 'yield_overflow': True, 'start': None, 'end': None, 'contig': None, 'every_fragment_as_molecule': False, 'skip_contigs': set(), 'progress_callback_function': None, 'pooling_method': 1, 'perform_allele_clustering': False}
"""
    # remove existing files generated by the previous test run if any
    for item in ["/tmp/tagged.bam", "/tmp/tagged.bam.bai", "/tmp/tagged.status.txt"]:
        if os.path.exists(item):
            os.remove(item)
    import tests

    mini_bam = os.path.join(
        os.path.dirname(inspect.getfile(tests)), "data", "mini_nla_test.bam"
    )
    # TODO: I was aware mini_nla_test.bam.bai differs from the upstream
    runner = CliRunner()
    result = runner.invoke(main, ["tag", "--bamfile", mini_bam])
    assert result.output == expected


@pytest.mark.skip(reason="still work in progress")
def test_bam():
    runner = CliRunner()
    result = runner.invoke(prepare.prepare, ["--help"])
    assert result.output == "prepare subcommand is called, yay!n"


@pytest.mark.skip(reason="still work in progress")
def test_end_motif():
    runner = CliRunner()
    result = runner.invoke(end_motif.end_motif)
    assert result.output == "end motif."
