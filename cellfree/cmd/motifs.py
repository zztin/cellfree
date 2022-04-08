import click


@click.command()
@click.option("--bamfile", help="Single bamfile")
def main(bamfile):
    print("motif command is called")
