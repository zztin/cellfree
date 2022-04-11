import click
from cellfree.prepare.prepare_bam import prepare_bam_single_thread


@click.command()
@click.option("--out-bam-path", help="Path to output bam file")
@click.argument("input_bam_path")
def prepare(input_bam_path,
            out_bam_path,
            ):


    prepare_bam_single_thread(input_bam_path,
                              out_bam_path,
                              consensus_model=None,
                              consensus_model_args={},
                              ignore_bam_issues=False,
                              head=None,
                              no_source_reads=False
                              )

    print("prepare subcommand is called, yay!")
