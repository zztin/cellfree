import sys
import pysam
from singlecellmultiomics.molecule import MoleculeIterator, ReadIterator
from singlecellmultiomics.bamProcessing.bamFunctions import sorted_bam_file, get_reference_from_pysam_alignmentFile, write_program_tag, MapabilityReader, verify_and_fix_bam,add_blacklisted_region
from singlecellmultiomics.universalBamTagger.tagging import generate_tasks, prefetch, run_tagging_tasks, write_job_gen_to_bed
import cellfree
from datetime import datetime
from itertools import chain
from singlecellmultiomics.molecule.consensus import calculate_consensus




def write_status(output_path, message):
    status_path = output_path.replace('.bam','.status.txt')
    with open(status_path,'w') as o:
        o.write(message+'\n')


def prepare_bam_single_thread(
        input_bam_path,
        out_bam_path,
        molecule_iterator = None,
        molecule_iterator_args = None,
        consensus_model = None,
        consensus_model_args={}, # Clearly the consensus model class and arguments should be part of molecule
        ignore_bam_issues=False,
        head=None,
        no_source_reads=False
        ):

    input_bam = pysam.AlignmentFile(input_bam_path, "rb", ignore_truncation=ignore_bam_issues, threads=4)
    input_header = input_bam.header.as_dict()


    # Write provenance information to BAM header
    write_program_tag(
        input_header,
        program_name='cellfree',
        command_line=" ".join(
            sys.argv),
        version=cellfree.__version__,
        description=f'SingleCellMultiOmics molecule processing, executed at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

    print(f'Started writing to {out_bam_path}')


    molecule_iterator_args = prefetch(molecule_iterator_args['contig'],
                                    molecule_iterator_args['start'],
                                    molecule_iterator_args['end'],
                                    molecule_iterator_args['start'],
                                    molecule_iterator_args['end'],
                                    molecule_iterator_args)


    molecule_iterator_args_wo_alignment = {k:v for k, v in molecule_iterator_args.items() if k != 'alignments'}
    molecule_iterator_args_wo_alignment_unmapped = molecule_iterator_args_wo_alignment.copy()
    molecule_iterator_args_wo_alignment_unmapped['contig'] = '*'

    molecule_iterator_exec = chain(
                molecule_iterator(input_bam, **molecule_iterator_args_wo_alignment_unmapped), molecule_iterator(input_bam, **molecule_iterator_args_wo_alignment),
                )


    print('Params:',molecule_iterator_args)
    read_groups = dict()  # Store unique read groups in this dict

    with sorted_bam_file(out_bam_path, header=input_header, read_groups=read_groups) as out:
        try:
            for i, molecule in enumerate(molecule_iterator_exec):

                # Stop when enough molecules are processed
                if head is not None and (i - 1) >= head:
                    break

                # set unique molecule identifier
                molecule.set_meta('mi', f'{molecule.get_a_reference_id()}_{i}')

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
                    calculate_consensus(molecule,
                                        consensus_model,
                                        i,
                                        out,
                                        **consensus_model_args)

                # Write the reads to the output file
                if not no_source_reads:
                    molecule.write_pysam(out)
        except Exception as e:
            write_status(out_bam_path,'FAIL, The file is not complete')
            raise e

        # Reached the end of the generator
        write_status(out_bam_path,'Reached end. All ok!')

