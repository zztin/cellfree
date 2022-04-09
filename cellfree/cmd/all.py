import click
from cellfree.algorithm import end_motif


@click.command()
@click.option("--prepared-bamfile", help="Single bamfile")
def all(prepared_bamfile):
    end_motif.run(prepared_bamfile)
