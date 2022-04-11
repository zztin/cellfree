import click


@click.command()
@click.argument("tagged-bamfile")
@click.option("--end-motif", is_flag=True, show_default=True, default=False,
              help="Insert result to the end-motif column in the tsv")
@click.option("--frag", is_flag=True, show_default=True, default=False,
              help="Insert result to the frag column in the tsv")
def tsv(tagged_bamfile, end_motif, frag):
    print(f"Tagged Bamfile: {tagged_bamfile}")
    if end_motif:
        print(f"Enabled:        end_motif")
    if frag:
        print(f"Enabled:        frag")