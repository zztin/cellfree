import sys
from datetime import datetime
from itertools import chain

import pysam
from singlecellmultiomics.bamProcessing.bamFunctions import (
    MapabilityReader,
    add_blacklisted_region,
    get_reference_from_pysam_alignmentFile,
    sorted_bam_file,
    verify_and_fix_bam,
    write_program_tag,
)
from singlecellmultiomics.molecule import MoleculeIterator, ReadIterator
from singlecellmultiomics.molecule.consensus import calculate_consensus
from singlecellmultiomics.universalBamTagger.tagging import (
    generate_tasks,
    prefetch,
    run_tagging_tasks,
    write_job_gen_to_bed,
)

import cellfree
from cellfree.molecule.cyclomics import CHICMolecule
from cellfree.read_unit.chic import CHICFragment


def write_status(output_path, message):
    status_path = output_path.replace(".bam", ".status.txt")
    with open(status_path, "w") as o:
        o.write(message + "\n")


def add_read_group_SM_tag(read_bam, bam_with_rg, SM_tag):
    with pysam.AlignmentFile(read_bam, "r", ignore_truncation=True) as g:
        with sorted_bam_file(bam_with_rg, origin_bam=g) as f:
            for i, read in enumerate(g):
                read.set_tag("SM", SM_tag)
                f.write(read)
    return bam_with_rg


def prepare_bam_single_thread(
    input_bam_path,
    out_bam_path,
    consensus_model=None,
    consensus_model_args={},  # Clearly the consensus model class and arguments should be part of molecule
    ignore_bam_issues=False,
    head=None,
    no_source_reads=False,
):

    input_bam = pysam.AlignmentFile(
        input_bam_path, "rb", ignore_truncation=ignore_bam_issues, threads=4
    )

    input_header = input_bam.header.as_dict()
    reference = get_reference_from_pysam_alignmentFile(input_bam_path)

    # Write provenance information to BAM header
    write_program_tag(
        input_header,
        program_name="cellfree",
        command_line=" ".join(sys.argv),
        version=cellfree.__version__,
        description=f'Cellfree molecule processing, executed at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}',
    )

    print(f"Started writing to {out_bam_path}")

    read_groups = dict()  # Store unique read groups in this dict

    with sorted_bam_file(
        out_bam_path, header=input_header, read_groups=read_groups
    ) as out:

        molecule_iterator_exec = chain(
            MoleculeIterator(
                alignments=input_bam,
                moleculeClass=CHICMolecule,
                fragmentClass=CHICFragment,
                every_fragment_as_molecule=False,
                perform_qflag=False,
                molecule_class_args={
                    "reference": reference,
                    "max_associated_fragments": 200,
                },
                fragment_class_args={"assignment_radius": 4},
                max_buffer_size=1000000,
                yield_overflow=False,
            )
        )
        try:
            for i, molecule in enumerate(molecule_iterator_exec):

                # Stop when enough molecules are processed
                if head is not None and (i - 1) >= head:
                    break

                # set unique molecule identifier
                molecule.set_meta("mi", f"{molecule.get_a_reference_id()}_{i}")

                # Write tag values
                molecule.write_tags()

                """
                if unphased_allele_resolver is not None:  # write unphased allele tag:
                    molecule.write_allele_phasing_information_tag(
                        unphased_allele_resolver, 'ua')
                """

                # Update read groups
                for fragment in molecule:
                    rgid = fragment.get_read_group()
                    if not rgid in read_groups:
                        read_groups[rgid] = fragment.get_read_group(True)[1]

                # Calculate molecule consensus
                if consensus_model is not None:
                    calculate_consensus(
                        molecule, consensus_model, i, out, **consensus_model_args
                    )

                # Write the reads to the output file
                if not no_source_reads:
                    molecule.write_pysam(out)
        except Exception as e:
            write_status(out_bam_path, "FAIL, The file is not complete")
            raise e

        # Reached the end of the generator
        write_status(out_bam_path, "Reached end. All ok!")


def read_prepared_bam_to_molecules(
    input_prepared_bam,
    reference=None,
    head=None,
    ignore_bam_issues=False,
):
    """
    No deduplication needed. Every fragment as a molecule.
    """
    input_bam = pysam.AlignmentFile(
        input_prepared_bam, "rb", ignore_truncation=ignore_bam_issues, threads=4
    )
    molecules = []
    molecule_iterator_exec = chain(
        MoleculeIterator(
            alignments=input_bam,
            moleculeClass=CHICMolecule,
            fragmentClass=CHICFragment,
            every_fragment_as_molecule=True,
            perform_qflag=False,
            max_buffer_size=1000000,
            yield_overflow=False,
        )
    )

    try:
        for i, molecule in enumerate(molecule_iterator_exec):

            # Stop when enough molecules are processed
            if head is not None and (i - 1) >= head:
                break
            molecules.append(molecule)
    except Exception as e:
        print("FAIL, The file is not complete")
        raise e

    return molecules
