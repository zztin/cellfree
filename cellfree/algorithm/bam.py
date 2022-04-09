from singlecellmultiomics.universalBamTagger.bamtagmultiome import argparser
from singlecellmultiomics.universalBamTagger.bamtagmultiome import run_multiome_tagging


def run_tagging(bamfile):
    args = argparser.parse_args()
    args.bamin = bamfile
    args.method = "nla"
    args.o = "/tmp/tagged.bam"
    run_multiome_tagging(args)
