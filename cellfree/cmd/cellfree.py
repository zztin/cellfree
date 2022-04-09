import click

from cellfree.cmd import prepare, end_motif, all


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


main.add_command(prepare.main)
main.add_command(end_motif.main)
main.add_command(all.main)
