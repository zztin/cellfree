import click


@click.command()
@click.option("--bamfile", help="Single bamfile")
def main(bamfile):
    print("bam command is called")
