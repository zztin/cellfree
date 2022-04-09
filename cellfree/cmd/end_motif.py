import click
from cellfree.algorithm import end_motif


@click.command()
@click.option("--prepared-bamfile", help="Prepared file.")
@click.option("--png", help="Path to output statistics files")
def end_motif(prepared_bamfile):
    end_motif.run(prepared_bamfile)
