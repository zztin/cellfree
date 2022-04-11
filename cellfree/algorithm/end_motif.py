import collections
from dataclasses import dataclass
from cellfree.prepare.prepare_bam import read_prepared_bam_to_molecules
import pandas as pd


# def get_k_mers(molecule):
#     molecule.strand
#     molecule.start
#     molecule.end
#     fetch_sequence(molecule.sequence)
#     reference.fetch_coordinate




class Motifs:
    def __init__(self,
                 molecules,
                 bp_count=4,
                 table=None):
        if table is None:
            self.table = collections.defaultdict()
        else:
            self.table = table
        self.molecules = []


    def get_motifs(self):
        self.table[self.molecule.readname]['5_end']
        self.table[self.molecule.readname]['3_end']
        self.table[self.molecule.readname]['5_end_upstream']
        self.table[self.molecule.readname]['3_end_downstream']
        self.table[self.molecule.readname]['middle']


def run(prepared_bam):
    motifs = Motifs()

    for molecule in read_prepared_bam_to_molecules(prepared_bam):
        motifs.get_motifs(molecule)

    return

    print("motif command is called")
