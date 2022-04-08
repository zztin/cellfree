import pysam


class RunSummary:
    def __init__(self):
        raise NotImplementedError

    def collect(self):
        raise NotImplementedError


class EndMotifCollect(RunSummary):
    def __init__(self):
        self.rev_end_5_prime = None
        self.fwd_end_5_prime = None
        self.rev_end_3_prime = None
        self.fwd_end_3_prime = None
        raise NotImplementedError


class LengthCollect(RunSummary):
    def __init__(self):
        self.simple_molecule_length = None
        self.sv_molecules_length = None
        self.circular_molecules_length = None
        self.unmapped_molecules_length = None

        raise NotImplementedError


class CircularDnaCollect(RunSummary):
    def __init__(self):
        self.circular = None
        self.rev_end_3_prime = None
        self.fwd_end_3_prime = None
        raise NotImplementedError


class Export:
    def __init__(self):
        self.data = None
        self.format = None  # tsv, pickle.gz, ml_input
