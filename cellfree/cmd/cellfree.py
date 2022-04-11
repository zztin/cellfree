import click

from cellfree.cmd import prepare, end_motif, all, frag_length


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


main.add_command(prepare.prepare)
main.add_command(end_motif.end_motif)
main.add_command(frag_length.frag_length)
main.add_command(all.all)
