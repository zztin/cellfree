import glob

import click
import pysam
from singlecellmultiomics.bamProcessing.bamFunctions import (
    MapabilityReader,
    add_blacklisted_region,
    get_reference_from_pysam_alignmentFile,
    sorted_bam_file,
    verify_and_fix_bam,
    write_program_tag,
)

from cellfree.prepare.prepare_bam import prepare_bam_single_thread


@click.command()
@click.option("--out-bam-path", help="Path to output bam file")
@click.option("--bam-tag", default="BX", help="Path to output bam file")
@click.argument("input_bam_path")
def separate_barcoded(
    input_bam_path,
    out_bam_path,
    bam_tag,
    ignore_bam_issues=False,
):
    reads_no_bb = []
    reads_w_bb = []
    reads_weird = []
    input_bam = pysam.AlignmentFile(
        input_bam_path, "rb", ignore_truncation=ignore_bam_issues, threads=4
    )
    input_header = input_bam.header.as_dict()
    for x in [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]:
        for i, read in enumerate(input_bam.fetch(x, 0, 2140000000)):
            # This can be breaking if new consensus reads format is created.
            if not read.query_name.startswith("consensus"):
                pass
            elif read.get_tag(bam_tag) == "NotFound":
                reads_no_bb.append(read)
            elif read.get_tag(bam_tag) != "XXXX":
                reads_w_bb.append(read)
            else:
                reads_weird.append(read)

    read_groups = dict()  # Store unique read groups in this dict
    with sorted_bam_file(
        out_bam_path + "_no_bb.sorted.bam", header=input_header, read_groups=read_groups
    ) as out:
        for read in reads_no_bb:
            out.write(read)

    with sorted_bam_file(
        out_bam_path + "_with_bb.sorted.bam",
        header=input_header,
        read_groups=read_groups,
    ) as out:
        for read in reads_w_bb:
            out.write(read)

    with sorted_bam_file(
        out_bam_path + "_imcomplete_bb.sorted.bam",
        header=input_header,
        read_groups=read_groups,
    ) as out:
        for read in reads_weird:
            out.write(read)

    print("Files are separated.")
    c_no_b = pysam.view("-c", out_bam_path + "_no_bb.sorted.bam")
    c_w_b = pysam.view("-c", out_bam_path + "_with_bb.sorted.bam")
    c_i = pysam.view("-c", out_bam_path + "_imcomplete_bb.sorted.bam")

    print(
        f"read count of {out_bam_path.split('/')[-1]}: \n\t_no_bb: {c_no_b}\t_with_bb:{c_w_b}\timcomplete_bb:{c_i}"
    )
