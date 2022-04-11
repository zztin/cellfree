import click
from cellfree.algorithm import cal_frag_length
from cellfree.plot.plot_lengths import plot_forward_length

@click.command()
@click.argument("prepared_bamfile")
@click.option("--png", help="Path to output statistics files")
def frag_length(prepared_bamfile, png):
    counters = cal_frag_length.run(prepared_bamfile)
    print("fwd_read_lengths", counters[0])
    plot_forward_length(counters[0], png)
