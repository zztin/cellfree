from pysam import FastaFile
from pysamiterators import CachedFasta


class CachedFastaNoHandle(CachedFasta):
    def __init__(self, path: str):
        handle = FastaFile(path)
        CachedFasta.__init__(self, handle)


class FastaNoHandle(FastaFile):
    def __init__(self, path: str):
        handle = FastaFile(path)
        FastaFile.__init__(self, handle)
