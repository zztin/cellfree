import click
from cellfree.algorithm import end_motif


@click.command()
@click.option("--bamfile", help="Single bamfile")
def main(bamfile):
    end_motif.run(bamfile)
