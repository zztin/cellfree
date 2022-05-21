#!/usr/bin/env python3
import argparse
import collections
import json

import pysam
from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file

# Supress [E::idx_find_and_load] Could not retrieve index file for, see https://github.com/pysam-developers/pysam/issues/939
pysam.set_verbosity(0)


def load_metadata(json_path):
    # load input:
    with open(json_path, "r") as f:
        data_raw = json.load(f)
    metadata = collections.defaultdict(dict)
    for entry in data_raw:
        # TODO: Change in Cycas to make metadata id same as bam file query name.
        metadata[entry["id"].rsplit("_", 3)[0]] = entry
    return metadata


def load_bam_file(bam_path):
    bam = pysam.AlignmentFile(bam_path)
    return bam


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="Path to bam file")
    parser.add_argument(
        "json", type=str, help="Path to json file containing per read information"
    )
    parser.add_argument("sample_name", type=str, help="Sample name for SM tag.")
    parser.add_argument(
        "out_path", type=str, default="./", help="Path to out tagged bam file"
    )
    parser.add_argument(
        "-li",
        "--library",
        type=str,
        default="CyclomicsSeq",
        help="Type of library prep",
    )
    parser.add_argument(
        "-pl", "--platform", type=str, default="Nanopore", help="Sequencing platform"
    )
    parser.add_argument(
        "-v", "--version", type=int, help="consensus version may influence json format."
    )
    args = parser.parse_args()
    # load metadata
    metadata = load_metadata(args.json)
    # load bam
    bam = load_bam_file(args.bam)
    input_header = bam.header.as_dict()
    # Get input filename without suffix.
    name_suffix = args.bam.rsplit("/", 1)[-1].split(".bam")[0]
    # collect reads of different classes
    no_backbone = []
    with_backbone = []
    backbone_summary = collections.Counter()
    for i, read in enumerate(bam):
        read_name = read.query_name
        # Add SM tag
        read.set_tag("SM", args.sample_name)
        # Add Barcode tag Y_tags
        # Tag YB contains barcode
        for Y_tag, value in metadata[read_name]["Y_tags"].items():
            read.set_tag(Y_tag, value)
        if read.get_tag("YB") == "NNNN":
            no_backbone.append(read)
        else:
            with_backbone.append(read)
            backbone_summary[read.get_tag("YB")] += 1
        # TODO: Add other tags - (a) nanopore sequencing info, (b) start and end of each unit
    # Write file no backbone
    read_groups = dict()
    with sorted_bam_file(
        f"{args.out_path}/{name_suffix}_no_bb.tagged.sorted.bam",
        header=input_header,
        read_groups=read_groups,
    ) as out:
        for read in no_backbone:
            out.write(read)
    # Write file with backbone
    read_groups = dict()
    with sorted_bam_file(
        f"{args.out_path}/{name_suffix}_w_bb.tagged.sorted.bam",
        header=input_header,
        read_groups=read_groups,
    ) as out:
        for read in with_backbone:
            out.write(read)
