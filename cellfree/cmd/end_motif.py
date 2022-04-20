import click
from pysam import FastaFile

from cellfree.algorithm import get_end_motif
from cellfree.plot.plot_end_motifs import plot_kmer_dist
from cellfree.utils.pyutils import examine_path


@click.command()
@click.argument("prepared_bamfile")
@click.option("--out-path", default="~/00_projects/cellfree-working/")
@click.option("--sample", default="TEST")
@click.option("--refpath", default="~/00")
def end_motif(prepared_bamfile, out_path, sample, refpath):
    reference = FastaFile(refpath)
    table = get_end_motif.run(prepared_bamfile, reference, table=None)
    # check if path exist, if not, create folder
    examine_path(f"{out_path}/feature_table/")
    examine_path(f"{out_path}/end_motifs/")
    table.to_csv(f"{out_path}/feature_table/{sample}.tsv", sep="\t")
    print("table successfully stored.")
    # plot_kmer_dist(f"{out_path}/end_motifs/{sample}.png",
    #                table['arr_5_end'],
    #                '5 prime end',
    #                ['red'],
    #                )
