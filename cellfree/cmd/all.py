import sys

import click
from pysam import FastaFile

from cellfree.algorithm import get_end_motif
from cellfree.utils.pyutils import examine_path


@click.command()
@click.argument("prepared_bamfile")
@click.option("--out-path")
@click.option("--sample", default="TEST")
@click.option("--refpath", help="Reference file is necessary to get end motif k-mer.")
@click.option("--add-label", default=None, type=str)
def all(prepared_bamfile, out_path, sample, refpath, add_label):
    if refpath is None:
        print(
            "[Error] Please add a path to reference genome file by supplying --refpath and resubmit."
        )
        sys.exit(1)
    else:
        reference = FastaFile(refpath)

    table = get_end_motif.run(prepared_bamfile, reference, table=None)
    table["label"] = add_label
    print("Molecule count:", table.shape[0])
    # check if path exist, if not, create folder
    examine_path(f"{out_path}/feature_table/")
    table.to_pickle(f"{out_path}/feature_table/{sample}.pickle.gz")
    table.to_csv(f"{out_path}/feature_table/{sample}.tsv", index=False, sep="\t")

    print("table successfully stored.")
