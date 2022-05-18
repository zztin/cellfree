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
@click.argument("input_bam_path")
def separate_barcoded(
    input_bam_path,
    out_bam_path,
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
            if read.get_tag("BX") == "NotFound":
                reads_no_bb.append(read)
            elif read.get_tag("BX") != "XXXX":
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
