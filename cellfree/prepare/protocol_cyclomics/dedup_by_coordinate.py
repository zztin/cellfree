import argparse
import os
import time

import pysam
import singlecellmultiomics
from singlecellmultiomics.bamProcessing import sorted_bam_file
from singlecellmultiomics.bamProcessing.bamFunctions import (
    get_reference_from_pysam_alignmentFile,
)
from singlecellmultiomics.fastaProcessing import CachedFastaNoHandle
from singlecellmultiomics.utils.prefetch import UnitialisedClass

from cellfree.molecule.cyclomics import CHICMolecule

# from singlecellmultiomics.molecule.chic import CHICMolecule
from cellfree.molecule.iterator import MoleculeIterator

# from singlecellmultiomics.fragment import CHICFragment
from cellfree.read_unit.chic import CHICFragment

if __name__ == "__main__":
    pid = os.getpid()
    parser = argparse.ArgumentParser(
        description="Deduplicate base on molecular location."
    )
    parser.add_argument("--read_bam", type=str, help="bam file path")
    parser.add_argument("--out_path", type=str, help="out path")
    parser.add_argument("--prefix", type=str, help="prefix of the bam file name")
    parser.add_argument("--SM", type=str, help="value to put in SM tag")
    parser.add_argument(
        "--ref", type=str, help="reference which bam file is mapped to. (str)"
    )
    parser.add_argument(
        "--merge",
        type=str,
        help="merge reads from different nanopore reads but covering the same start,"
        "end sites within max 2bp range. Direction of the read is set by the first added read.",
    )

    args = parser.parse_args()
    SM_bam = f"{args.out_path}/{args.prefix}_{pid}_test.SMtagged.sorted.bam"
    t_bam = f"{args.out_path}/{args.prefix}_{pid}_coordinate_dedup.tagged.bam"

    # autodetect reference:
    reference = None
    if args.ref is None:
        args.ref = get_reference_from_pysam_alignmentFile(args.read_bam)

    if args.ref is not None:
        try:
            reference = UnitialisedClass(CachedFastaNoHandle, args.ref)
            print(f"Loaded reference from {args.ref}")
        except Exception as e:
            print(
                "Error when loading the reference file, continuing without a reference"
            )
            reference = None

    if args.SM is None:
        with pysam.AlignmentFile(args.read_bam, "r", ignore_truncation=True) as g:
            with sorted_bam_file(SM_bam, origin_bam=g) as f:
                for i, read in enumerate(g):
                    read.set_tag("SM", args.SM)
                    f.write(read)
    else:
        SM_bam = args.read_bam
    timeA = time.time()
    print("SM tag written.", timeA)

    with pysam.AlignmentFile(SM_bam) as f:

        with sorted_bam_file(
            t_bam,
            origin_bam=f,
        ) as target_bam:
            for i, m in enumerate(
                MoleculeIterator(
                    alignments=f,
                    moleculeClass=CHICMolecule,
                    fragmentClass=CHICFragment,
                    every_fragment_as_molecule=False,
                    perform_qflag=False,
                    molecule_class_args={
                        "reference": reference,
                        "max_associated_fragments": 100,
                    },
                    fragment_class_args={"assignment_radius": 20},
                    max_buffer_size=1000000,
                    yield_overflow=False,
                )
            ):
                read_name = f"consensus_{m.get_a_reference_id()}_{i}"
                # write tags to all fragments associated with the molecule
                m.write_tags()

                m.write_pysam(
                    target_bam,
                    consensus=True,
                    consensus_name=read_name,
                    no_source_reads=True,
                )

    pysam.index(t_bam, f"{t_bam}.bai")
    if args.SM is not None:
        os.remove(SM_bam)
        os.remove(f"{SM_bam}.bai")

print((time.time() - timeA) / 60, "min")
