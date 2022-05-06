import collections
from dataclasses import dataclass

import pandas as pd

from cellfree.prepare.prepare_bam import read_prepared_bam_to_molecules

# def get_k_mers(molecule):
#     molecule.strand
#     molecule.start
#     molecule.end
#     fetch_sequence(molecule.sequence)
#     reference.fetch_coordinate


def get_fragment_neighbouring_sequence(molecule, extra_bp=4, reference=None):
    """Obtain the sequence between the start and end of the molecule

    Args:
        reference(pysam.FastaFile) : reference  to use.
            If not specified `self.reference` is used

    Returns:
        sequence (str)
    """
    if reference is None:
        if molecule.reference is None:
            raise ValueError("Please supply a reference (PySAM.FastaFile)")
    try:
        span_plus_extra = reference.fetch(
            molecule.chromosome,
            molecule.spanStart - extra_bp,
            molecule.spanEnd + extra_bp,
        ).upper()
        return span_plus_extra

    except ValueError as e:
        print(
            f"Error: {e} at position {molecule.chromosome}:{molecule.spanStart}-{molecule.spanEnd}. XXXX used as replacement."
        )
        return "X" * 1000


class Motifs:
    def __init__(self, molecules, bp_count=4, table=None, reference=None):
        self.bp_count = bp_count
        if table is None:
            self.table = pd.DataFrame()
        else:
            self.table = table
        self.molecules = molecules
        self.reference = reference

    def get_feature_table(self):
        # TODO: check all molecules are mapped
        # use arrays
        # arr_read_name = ['readname']
        # arr_dir = ['dir_true_is_rev']
        # arr_read_length = ['length']
        # arr_5_end = ['end_5_prime']
        # arr_3_end = ['end_3_prime']
        # arr_5_end_upstream = ['end_5_upstream']
        # arr_3_end_downstream = ['end_3_upstream']
        arr_read_name = []
        arr_chromosome = []
        arr_start = []
        arr_end = []
        arr_dir = []
        arr_read_length = []
        arr_5_end = []
        arr_3_end = []
        arr_5_end_upstream = []
        arr_3_end_downstream = []
        arr_score = []

        for i, molecule in enumerate(self.molecules):
            arr_read_name.append(f"m{i}")
            assert type(molecule.chromosome) is str
            arr_chromosome.append(str(molecule.chromosome))
            arr_start.append(molecule.spanStart)
            arr_end.append(molecule.spanEnd)
            assert type(molecule.strand) is bool
            arr_dir += ["-" if molecule.strand else "+"]
            arr_read_length.append(molecule.spanEnd - molecule.spanStart)
            frag_seq = get_fragment_neighbouring_sequence(
                molecule, self.bp_count, self.reference
            )
            arr_5_end.append(frag_seq[self.bp_count + 1 : self.bp_count * 2 + 1])
            arr_5_end_upstream.append(frag_seq[: self.bp_count])
            arr_3_end.append(frag_seq[-self.bp_count * 2 : -self.bp_count])
            arr_3_end_downstream.append(frag_seq[-self.bp_count :])
            # self.table[molecule.readname][f"5_end_{self.bp_count}"] = molecule.get_fragment_span_sequence[:4]
            # self.table[molecule.readname][f"3_end_{self.bp_count}"] = molecule.get_fragment_span_sequence[-4:]
            # self.table[molecule.readname][f"5_end_upstream_{self.bp_count}"] =
            # self.table[molecule.readname][f"3_end_downstream_{self.bp_count}"] = molecule.get_fragment_neighbouring_sequence[-4:]
        self.table["chrom"] = arr_chromosome
        self.table["chromStart"] = arr_start
        self.table["chromEnd"] = arr_end
        self.table["name"] = arr_read_name
        self.table["score"] = 0
        self.table["strand"] = arr_dir
        self.table["arr_read_length"] = arr_read_length
        self.table["arr_5_end_upstream"] = arr_5_end_upstream
        self.table["arr_5_end"] = arr_5_end
        self.table["arr_3_end"] = arr_3_end
        self.table["arr_3_end_downstream"] = arr_3_end_downstream


def run(prepared_bam, reference, table=None, out_path="../test_out.tsv"):
    """
    out_path can be ‘.gz’, ‘.bz2’, ‘.zip’, ‘.xz’, or ‘.zst’ , inferred.
    """
    molecules = read_prepared_bam_to_molecules(prepared_bam)
    motifs = Motifs(
        molecules,
        bp_count=4,
        reference=reference,
        table=table,
    )
    motifs.get_feature_table()

    return motifs.table
