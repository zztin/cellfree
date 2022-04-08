import click

from cellfree.cmd import bam, motifs


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


main.add_command(bam.main)
main.add_command(motifs.main)
