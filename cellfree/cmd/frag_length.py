import click
from cellfree.algorithm import cal_frag_length
from cellfree.plot.plot_lengths import plot_forward_length, plot_reverse_length

@click.command()
@click.argument("prepared_bamfile")
@click.option("--png", help="Path to output statistics files")
@click.option("--sample-name", default="TEST")
def frag_length(prepared_bamfile, png, sample_name):
    counters = cal_frag_length.run(prepared_bamfile)
    plot_forward_length(counters[0], png, sample_name=sample_name)
    plot_reverse_length(counters[1], png, sample_name=sample_name)
