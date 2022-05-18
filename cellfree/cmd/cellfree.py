import click

from cellfree.algorithm import bam
from cellfree.cmd import all, end_motif, frag_length, prepare, split, tsv


@click.group()
@click.option(
    "--verbose",
    type=click.Choice(["debug", "info", "warning", "error", "critical"]),
    default="info",
    help="Verbose level corresponding to logging level",
)
@click.option("--conf", help="Configuration file")
def main(verbose, conf):
    pass


@click.command()
@click.argument("bam_file")
def tag(bam_file):
    bam.run_tagging(bam_file)


main.add_command(prepare.prepare)
main.add_command(tag)
main.add_command(end_motif.end_motif)
main.add_command(frag_length.frag_length)
main.add_command(tsv.tsv)
main.add_command(all.all)
main.add_command(split.separate_barcoded)
