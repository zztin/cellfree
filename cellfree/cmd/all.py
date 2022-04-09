import click
from cellfree.algorithm import end_motif


@click.command()
@click.option("--prepared-bamfile", help="Single bamfile")
def main(prepared_bamfile):
    end_motif.run(prepared_bamfile)
