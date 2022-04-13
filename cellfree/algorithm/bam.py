from singlecellmultiomics.universalBamTagger.bamtagmultiome import (
    argparser,
    run_multiome_tagging,
)


def run_tagging(bamfile):
    # TODO: this will result in the commandline is parsed by argparser but not click. In the end we would like to
    #  use this argparser workaround for porting
    args = argparser.parse_args()
    args.bamin = bamfile
    args.method = "nla"
    args.o = "/tmp/tagged.bam"
    run_multiome_tagging(args)
