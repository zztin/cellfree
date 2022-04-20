import click
from pysam import FastaFile

from cellfree.algorithm import get_end_motif
from cellfree.utils.pyutils import examine_path


@click.command()
@click.argument("prepared_bamfile")
@click.option("--out-path", default="~/00_projects/cellfree-working/")
@click.option("--sample", default="TEST")
@click.option("--refpath", default="~/00_projects/cellfree-working/")
@click.option("--add-label", default=None, type=str)
def all(prepared_bamfile, out_path, sample, refpath, add_label):
    reference = FastaFile(refpath)
    table = get_end_motif.run(prepared_bamfile, reference, table=None)
    table["label"] = add_label
    print("Molecule count:", table.shape[0])
    # check if path exist, if not, create folder
    examine_path(f"{out_path}/feature_table/")
    table.to_csv(f"{out_path}/feature_table/{sample}.tsv", sep="\t")
    table.to_pickle(f"{out_path}/feature_table/{sample}.pickle.gz")
    print("table successfully stored.")
