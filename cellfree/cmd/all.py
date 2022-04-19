import click
from cellfree.algorithm import get_end_motif
from cellfree.utils.pyutils import examine_path
from pysam import FastaFile

@click.command()
@click.argument("prepared_bamfile")
@click.option("--out-path", default="~/00_projects/cellfree-working/")
@click.option("--sample", default="TEST")
@click.option("--refpath", default="~/00")
def all(prepared_bamfile, out_path, sample, refpath):
    reference = FastaFile(refpath)
    table = get_end_motif.run(prepared_bamfile, reference, table=None)
    # check if path exist, if not, create folder
    examine_path(f"{out_path}/feature_table/")
    table.to_csv(f"{out_path}/feature_table/{sample}.tsv", sep='\t')
    print("table successfully stored.")
