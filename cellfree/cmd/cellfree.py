import click

from cellfree.algorithm import bam
from cellfree.cmd import prepare, end_motif, tsv, all


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
@click.option("--bamfile", help="Single bamfile")
def tag(bamfile):
    bam.run_tagging(bamfile)


main.add_command(prepare.prepare)
main.add_command(tag)
main.add_command(end_motif.end_motif)
main.add_command(tsv.tsv)
main.add_command(all.all)

