import click


@click.command()
@click.option("--prepared-bamfile", help="Single bamfile")
def prepare(prepared_bamfile):
    print("bam command is called")
