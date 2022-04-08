# A molecule can contain multiple reads mapping at the same location. (PCR/ amplification Deduplication)

# A molecule that mapped to different locations are molecule that contains Structure Variants (SV). These
# contains reads with same read name but different mapping location.

# molecule are some times called as cell-free DNA fragment.
# molecule method contains consensus building from multiple reads mapping at same location
# molecule methods calls SV while building molecule

# Molecule contains read1 and read2. Start and end coordinate of the "fragment"

class Molecule:
    def __init__(self, width, height):

        raise NotImplementedError

class DnaMolecule(Molecule):
    def __init__(self, width, height):
        fragments = [] # list
        molecule_base_quality = []
        molecule_mapping_quality = []
        map_coordinate = None
        strand = None
        adapter_sequence = None # to check with soft clip bases
        softclips = None  # list of bases
        sequence = None
        genome_version = None # should be the same as the tracks. Should be checked a level up.
        circular = False # bool
        unmapped = False
        raise NotImplementedError

    def collect_fragments(self):
        raise NotImplementedError

    def get_strand(self):
        raise NotImplementedError

    def get_ends(self):
        self.end_3_prime = None
        self.end_5_prime = None
        self.end_3_prime_outside = None
        self.end_5_prime_outside = None

        raise NotImplementedError

    def get_molecule_length(self):
        raise NotImplementedError

    def get_genomics_tracks_overlay(self):
        raise NotImplementedError

    def is_circular(self):
        raise NotImplementedError



class DnaBisulfideMolecule(DnaMolecule):
    def __init__(self, width, height):
        raise NotImplementedError

    def get_methyl_sites(self):
        raise NotImplementedError


class DnaNanoporeMolecule(DnaMolecule):
    def __init__(self, width, height):
        raise NotImplementedError

    def get_methyl_sites(self):
        raise NotImplementedError


class RnaMolecule(Molecule):
    def __init__(self, width, height):
        raise NotImplementedError
