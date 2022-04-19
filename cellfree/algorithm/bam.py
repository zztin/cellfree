from singlecellmultiomics_ng.universalBamTagger.bamtagmultiome import (
    run_multiome_tagging_cmd,
)


def run_tagging(bamfile):
    # TODO: this will result in the commandline is parsed by argparser but not click. In the end we would like to
    #  use this argparser workaround for porting
    cmd = f"{bamfile} -method nla -o /tmp/tagged.bam".split(" ")
    run_multiome_tagging_cmd(cmd)
