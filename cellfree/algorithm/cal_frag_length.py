from collections import Counter
from dataclasses import dataclass

from cellfree.prepare.prepare_bam import read_prepared_bam_to_molecules


def run(prepared_bam):
    """
    strand (bool) : False for Forward, True for reverse

    """
    forward_strand_length = []
    reverse_strand_length = []
    unmapped_length = []
    molecules = read_prepared_bam_to_molecules(prepared_bam)
    for molecule in molecules:
        # TODO: molecule.strand is True if molecule is on reverse strand (should have only 1 value,
        # TODO: if contains flip, combine into forward strand?)
        if not molecule.strand:
            forward_strand_length.append(abs(molecule.spanEnd - molecule.spanStart))
        else:
            reverse_strand_length.append(abs(molecule.spanEnd - molecule.spanStart))

    fwd_len_counter = Counter(forward_strand_length)
    rev_len_counter = Counter(reverse_strand_length)
    # unmapped_len_counter = Counter(unmapped_length)
    # all_mapped_len_counter = Counter(forward_strand_length) + Counter(reverse_strand_length)
    # all_len_counter = Counter(forward_strand_length) + Counter(reverse_strand_length) + Counter(unmapped_length)

    return (fwd_len_counter, rev_len_counter)

    # Alternatively samtools stats @RL
    print("motif command is called")


def export_length(molecule):
    # tsv = None
    pass
