import click


@click.command()
@click.option("--prepared-bamfile", help="Prepared file.")
@click.option("--png", help="Path to output statistics files")
def end_motif(prepared_bamfile):
    print("end motif.")
