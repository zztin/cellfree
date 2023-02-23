import argparse
import os
import re
import shlex
import subprocess
import time
import warnings

import cyvcf2
import numpy as np
import pysam

from cellfree.utils.pyutils import examine_path

with warnings.catch_warnings():
    warnings.simplefilter("ignore")


def check_vcf_format(in_bam):
    if in_bam.endswith("vcf"):
        print("VCF file is not indexed")
        print(
            f"bgzip -c {in_bam} > {in_bam}.gz",
        )
        print(f"tabix -p vcf {in_bam}.gz")
        exit(1)


def pileup_cleanup(pileup_record):
    """
    exceptions includes:
    ^ begin of a sequence
    $ after end of a sequence
    X-NNN (N amount depends on the following deletion)
    * deletion
    >>>pileup_cleanup("^ATT")
    "ATT"
    >>>pileup_cleanup("A-1N*T")
    "A*T"
    >>>pileup_cleanup("A-2N**")
    "A"
    >>>pileup_cleanup("*T-1N*")
    "T"
    >>>pileup_cleanup("T$")
    "T"
    >>>pileup_cleanup("G-3NNN***")
    "G"
    """
    pileup_record = re.sub("[0-9]", "", pileup_record)
    pileup_record = pileup_record.replace("-", "")
    pileup_record = pileup_record.replace("N", "")
    pileup_record = pileup_record.replace("^", "")
    pileup_record = pileup_record.replace("$", "")
    pileup_record = pileup_record.replace("*", "")

    return pileup_record


# Get only the reads that overlap with somatic sSNV. Can skip.
def get_overlapping_vcf_reads_bedtools(
    somatic_vcf, input_bam_path, output_overlap_path
):
    pass
    subprocess.call(f"touch {output_overlap_path}.touchtest")
    cmd = shlex.split(
        f"bedtools intersect -wa -a {input_bam_path}  -b {somatic_vcf} > {output_overlap_path}"
    )
    subprocess.call(cmd)


class Alignments:
    def __init__(self, bam, out_ref_bam, out_alt_bam, out_exception_bam):
        self.alignments = []
        # Use this to assign alignments to out bam is faster
        self.ordered_alignments_tags = []
        self.phased_tags = []  # -1, 0, 1
        self.template_bam = bam
        self.MUT_alignments = []
        self.WT_alignments = []
        self.exception_alignments = []
        self.not_phased_ref_alignments = []
        self.ref_out_bam = out_ref_bam
        self.alt_out_bam = out_alt_bam
        self.exception_bam = out_exception_bam
        self.associated_gsnv = []
        self.associated_gsnv_labels = []

    def __len__(self):
        return len(self.alignments)

    # DA: allele, DS: siteCoordinate, DB: allele for gSNV , DC: siteCoordinate for gSNV
    def write_records(self):
        # Write wild type alignments by looking at ordered alignment tags
        array = np.array(self.ordered_alignments_tags)
        indices = list(np.where(array == 0)[0])
        self.WT_alignments = map(self.alignments.__getitem__, indices)
        # Get MUT
        array = np.array(self.ordered_alignments_tags)
        indices = list(np.where(array > 0)[0])
        self.MUT_alignments = map(self.alignments.__getitem__, indices)

        array = np.array(self.ordered_alignments_tags)
        indices = list(np.where(array == -2)[0])
        self.exception_alignments = map(self.alignments.__getitem__, indices)

        for read in self.WT_alignments:
            self.ref_out_bam.write(read)
        for read in self.MUT_alignments:
            self.alt_out_bam.write(read)
        for read in self.exception_alignments:
            self.exception_bam.write(read)


class SomaticNucleotideVariants:
    def __init__(self, snv, criteria=None, phasing=False, cfdna_length_range=350):
        # inherit gsnv object from cyvcf2
        if criteria is None:
            criteria = {"FILTER": True, "QUAL": 30, "LEN": 1}
        self.snv = snv
        self.CHROM = snv.CHROM
        self.start = snv.start
        self.end = snv.end
        self.ALT = snv.ALT
        self.REF = snv.REF
        # Boolean. If we want to phase it or not.
        self.phasing = phasing
        # add customized features
        self.ssnv_length = snv.end - snv.start
        self.snv_type = None
        # A list containing only alignments that contains gsnv. Updated over time.
        self.pileup_alignments = []
        self.query_sequences = []
        self.ssnv_ref_alleles = []
        self.ssnv_alt_alleles = []
        self.criteria = criteria
        self.pass_filter = False
        self.phasing_info = None  # Get information if it should be ALT-REF or ALT-ALT for the tumor allele. From vcf?
        self.gsnv_ref_alleles = []
        self.gsnv_alt_alleles = []
        self.remaining_positions = None
        self.align_to_ref = []
        self.align_to_alt = {}
        self.related_gsnv = []
        self.filter()
        self.cfdna_length_range = cfdna_length_range
        self.related_gsnv = []

    def get_gsnvs(self):
        self.related_gsnv = []

    def snv_type(self):
        if self.ssnv_length > 1:
            self.snv_type = "indel"
        else:
            self.snv_type = "gsnv"

    def filter(self):
        if self.criteria["LEN"] == 0:
            if self.snv.FILTER == self.criteria["FILTER"]:
                self.pass_filter = True

        elif self.snv.FILTER == self.criteria["FILTER"]:
            if self.snv.QUAL >= self.criteria["QUAL"]:
                if self.ssnv_length == self.criteria["LEN"]:
                    self.pass_filter = True
        else:
            self.pass_filter = False


class GermlineSomaticNucleotideVariants(SomaticNucleotideVariants):
    def __init__(self, snv, criteria=None):
        # inherit gsnv object from cyvcf2
        super().__init__(snv, criteria)
        if criteria is None:
            criteria = {"FILTER": True, "QUAL": 30, "LEN": 0}
        self.ssnv = snv
        self.CHROM = snv.CHROM
        self.start = snv.start
        self.end = snv.end
        self.ALT = snv.ALT
        self.REF = snv.REF
        # add customized features
        self.ssnv_length = snv.end - snv.start
        self.snv_type = None
        # A list containing only alignments that contain gsnv. Updated over time.
        self.pileup_alignments = []
        self.query_sequences = []
        self.ssnv_ref_alleles = []
        self.ssnv_alt_alleles = []
        self.criteria = criteria
        self.pass_filter = False
        self.phasing_info = None  # Get information if it should be ALT-REF or ALT-ALT for the tumor allele. From vcf?
        self.gsnv_ref_alleles = []
        self.gsnv_alt_alleles = []
        self.remaining_positions = None
        self.align_to_ref = []
        self.align_to_alt = {}
        # All reads that are aligned to a gsnv. Write DA tag.
        self.aligned_to_snv = []
        self.filter()

    def snv_type(self):
        if self.ssnv_length > 1:
            self.snv_type = "indel"
        else:
            self.snv_type = "gsnv"

    def filter(self):
        if self.criteria["LEN"] == 0:
            if self.ssnv.FILTER == self.criteria["FILTER"]:
                self.pass_filter = True

        if self.ssnv.FILTER == self.criteria["FILTER"]:
            if self.ssnv.QUAL >= self.criteria["QUAL"]:
                if self.ssnv_length == self.criteria["LEN"]:
                    self.pass_filter = True
        else:
            self.pass_filter = False


# For gsnv. End: Get reads overlap with each allele
def get_neighbour_gsnv_info(gsnv, bam, out_alignment_records):
    # Getting allele of gsnv per alignment
    pileup_columns = bam.pileup(
        gsnv.CHROM, gsnv.start, gsnv.end, min_mapping_quality=0, min_base_quality=0
    )
    # all columns of reads that overlap with this particular position.
    gsnv.remaining_positions = gsnv.ssnv_length
    for pileup_column in pileup_columns:
        # Collect all bases from gsnv.start to gsnv.end, append to correct index
        if pileup_column.reference_pos == gsnv.start:
            gsnv.pileup_alignments = [
                pileup_column.pileups[i].alignment
                for i in range(len(pileup_column.pileups))
            ]

            query_sequences = pileup_column.get_query_sequences(add_indels=True)
            gsnv.query_sequences = [v.upper() for v in query_sequences]
            gsnv.remaining_positions -= 1
            pileup_column_at_start_pos = pileup_column
            if gsnv.remaining_positions == 0:
                break
                # to collect more positions
        elif (pileup_column.reference_pos < gsnv.end) & (
            pileup_column.reference_pos > gsnv.start
        ):
            # If there is a deletion at this position, the length of the query_sequences would be shorter.
            # It is not possible to know which base should stitch to which. (Same problem as Inez faced)
            query_sequences_next = pileup_column.get_query_sequences(add_indels=True)
            query_sequences_next = [v.upper() for v in query_sequences_next]
            # Edge case where the alignment does not cover the first base of the VCF record.
            if len(gsnv.query_sequences) == 0:
                # First position is not covered. The covered part cannot be used to determine which allele it is.
                return 0

            # if all reads covering the indel positions are not prematurely short, then process normally
            elif len(query_sequences_next) == len(gsnv.query_sequences):
                gsnv.query_sequences = [
                    a + b for (a, b) in zip(gsnv.query_sequences, query_sequences_next)
                ]
            # Edge cases where one of the bam record is much shorter, and did not cover the complete indel length
            else:
                # kick out the alignments that are shorter than the gsnv/indel length
                # because we cannot determine the allele there (enter group "reads_not_aligned")
                alignments = [
                    pileup_column.pileups[i].alignment
                    for i in range(len(pileup_column.pileups))
                ]
                # Find the alignment that is not in the first position. Any alignment that either start later,
                # or ends prematurely will not be in the final list
                alignments_to_kick_out = set(gsnv.pileup_alignments) - set(alignments)
                for alignment in alignments_to_kick_out:
                    index_to_kick_out = gsnv.pileup_alignments.index(alignment)
                    query_sequences.pop(index_to_kick_out)
                    gsnv.pileup_alignments.pop(index_to_kick_out)
                # Now the query_sequences only contains the bases of the remaining alignments
                gsnv.query_sequences = [
                    a + b for (a, b) in zip(query_sequences, query_sequences_next)
                ]
                gsnv.remaining_positions -= 1
                if gsnv.remaining_positions == 0:
                    break
    # phasing to gsnv
    if gsnv.remaining_positions == 0:
        for i, query_sequence in enumerate(gsnv.query_sequences):
            query_sequence = pileup_cleanup(query_sequence)
            if (query_sequence not in gsnv.ALT) and (query_sequence != gsnv.REF):
                gsnv_type = None
            elif query_sequence == gsnv.REF:
                gsnv_type = 0
            for j in range(len(gsnv.ALT)):
                allele_num = j + 1  # start counting by 1
                if query_sequence == gsnv.ALT[j]:
                    gsnv_type = allele_num

            out_alignment_records.associated_gsnv_labels.append(
                (
                    gsnv.pileup_alignments[i].qname,
                    gsnv.CHROM,
                    gsnv.start,
                    gsnv.end,
                    query_sequence,
                    gsnv_type,
                )
            )
            out_alignment_records.associated_gsnv.append(gsnv.pileup_alignments[i])
            # print("GSNV!", (gsnv.pileup_alignments[i].qname,
            #                 gsnv.CHROM,
            #                 gsnv.start,
            #                 gsnv.end,
            #                 query_sequence,
            #                 gsnv_type))


# For snv
def get_alignments_and_query_per_snv(snv, bam, out_alignment_records):

    pileup_columns = bam.pileup(
        snv.CHROM, snv.start, snv.end, min_mapping_quality=0, min_base_quality=0
    )
    # all columns of reads that overlap with this particular position.

    snv.remaining_positions = snv.ssnv_length
    # If no read pileup at this location, the for loop would finish.
    for pileup_column in pileup_columns:
        if pileup_column.reference_pos == snv.start:
            snv.pileup_alignments = [
                pileup_column.pileups[i].alignment
                for i in range(len(pileup_column.pileups))
            ]

            query_sequences = pileup_column.get_query_sequences(add_indels=True)
            snv.query_sequences = [v.upper() for v in query_sequences]
            snv.remaining_positions -= 1
            if snv.remaining_positions == 0:
                break
                # to collect more positions
        elif (pileup_column.reference_pos < snv.end) & (
            pileup_column.reference_pos > snv.start
        ):
            # If there is a deletion at this position, the length of the query_sequences would be shorter.
            # It is not possible to know which base should stitch to which. (Same problem as Inez faced)
            query_sequences_next = pileup_column.get_query_sequences(add_indels=True)
            query_sequences_next = [v.upper() for v in query_sequences_next]

            # Edge case where the alignment does not cover the first base of the VCF record.
            if len(snv.query_sequences) == 0:
                # First position is not covered. The covered part cannot be used to determine which allele it is.
                return 0
            # if all reads covering the indel positions are not prematurely short, then process normally
            elif len(query_sequences_next) == len(snv.query_sequences):
                snv.query_sequences = [
                    a + b for (a, b) in zip(snv.query_sequences, query_sequences_next)
                ]

            # Edge cases where one of the bam record is much shorter, and did not cover the complete indel length
            else:
                # kick out the alignments that are shorter than the gsnv/indel length
                # because we cannot determine the allele there (enter group "reads_not_aligned")
                alignments = [
                    pileup_column.pileups[i].alignment
                    for i in range(len(pileup_column.pileups))
                ]
                # Find the alignment that is not in the first position. Any alignment that either start later,
                # or ends prematurely will not be in the final list
                alignments_to_kick_out = set(snv.pileup_alignments) - set(alignments)
                for alignment in alignments_to_kick_out:
                    index_to_kick_out = snv.pileup_alignments.index(alignment)
                    snv.query_sequences.pop(index_to_kick_out)
                    snv.pileup_alignments.pop(index_to_kick_out)
                # Now the query_sequences only contains the bases of the remaining alignments
                snv.query_sequences = [
                    a + b for (a, b) in zip(snv.query_sequences, query_sequences_next)
                ]
            snv.remaining_positions -= 1
            if snv.remaining_positions == 0:
                break
        else:
            pass
    # If bam overlaps with SNV:
    if snv.remaining_positions == 0:
        ref_alt_id = np.full((len(snv.query_sequences)), -1)
        for i, query_sequence in enumerate(snv.query_sequences):
            query_sequence = pileup_cleanup(query_sequence)
            # print(snv.REF, snv.ALT, query_sequence)

            if (query_sequence not in snv.ALT) and (query_sequence != snv.REF):
                # Dep ref_alt_id
                ref_alt_id[i] = -2
                out_alignment_records.ordered_alignments_tags.append(-2)
                # Does not align to any known allele
                snv.pileup_alignments[i].set_tag(
                    "DP", f"{snv.CHROM}_{snv.start}_{snv.end}"
                )
                snv.pileup_alignments[i].set_tag("DS", None)
                snv.pileup_alignments[i].set_tag("DA", -1)

                # Updated alignment is added to a new list.
                out_alignment_records.alignments.append(snv.pileup_alignments[i])
                # For debugging:
                # print("NOT MATCHED\n", f"{snv.CHROM}:{snv.start}-{snv.end}\n",  snv.pileup_alignments[i])
                # Write all unmapped sequence to file. Very valuable.

            elif query_sequence == snv.REF:
                # Dep ref_alt_id
                ref_alt_id[i] = 0
                # Add to ordered list
                out_alignment_records.ordered_alignments_tags.append(0)
                # Update tag to a bam  alignmentRead
                snv.pileup_alignments[i].set_tag(
                    "DP", f"{snv.CHROM}_{snv.start}_{snv.end}"
                )
                snv.pileup_alignments[i].set_tag("DS", f"{snv.REF}")
                snv.pileup_alignments[i].set_tag("DA", 0)
                # Updated alignment is added to a new list.
                out_alignment_records.alignments.append(snv.pileup_alignments[i])
            for j in range(len(snv.ALT)):
                allele_num = j + 1  # start counting by 1
                if query_sequence == snv.ALT[j]:
                    out_alignment_records.ordered_alignments_tags.append(allele_num)
                    ref_alt_id[i] = allele_num
                    snv.pileup_alignments[i].set_tag(
                        "DP", f"{snv.CHROM}_{snv.start}_{snv.end}"
                    )
                    snv.pileup_alignments[i].set_tag("DS", f"{snv.ALT[j]}")
                    snv.pileup_alignments[i].set_tag("DA", allele_num)
                    out_alignment_records.alignments.append(snv.pileup_alignments[i])
                    # For debugging:
                    # print(f"ALT: {snv.CHROM}:{snv.start}-{snv.end}\n", snv.pileup_alignments[i])

        # print("Allele info:", ref_alt_id)

    else:
        pass
        # print("This SNV has no pileup bam read.")


def find_ssnv(
    bam_path,
    somatic_vcf_path,
    germline_vcf_path,
    out_ref_bam_path,
    out_alt_bam_path,
    out_ERR_bam_path,
    phasing=False,
):
    check_vcf_format(somatic_vcf_path)
    check_vcf_format(germline_vcf_path)
    somatic_vcf = cyvcf2.VCF(somatic_vcf_path)
    germline_vcf = cyvcf2.VCF(germline_vcf_path)
    bam = pysam.AlignmentFile(bam_path)
    pysam.index(bam_path)
    out_ref_bam = pysam.AlignmentFile(out_ref_bam_path, "wb", template=bam)
    out_alt_bam = pysam.AlignmentFile(out_alt_bam_path, "wb", template=bam)
    out_ERR_bam = pysam.AlignmentFile(out_ERR_bam_path, "wb", template=bam)
    count_filtered_ssnv = 0
    count_bam_read = 0
    allele_ref_counts = 0
    allele_alt_counts = 0
    allele_not_matched = 0
    allele_debug_should_be_zero = 0

    for ssnv_raw in somatic_vcf:
        out_alignment_records = Alignments(bam, out_ref_bam, out_alt_bam, out_ERR_bam)
        ssnv = SomaticNucleotideVariants(
            ssnv_raw, criteria={"FILTER": None, "QUAL": 30, "LEN": 0}, phasing=phasing
        )
        # if variant pass filter, action per variant.
        if ssnv.pass_filter:
            count_filtered_ssnv += 1
            get_alignments_and_query_per_snv(ssnv, bam, out_alignment_records)

            if ssnv.phasing:
                for gsnv_raw in germline_vcf(
                    f"{ssnv.CHROM}:{ssnv.start-ssnv.cfdna_length_range}-{ssnv.end + ssnv.cfdna_length_range}"
                ):
                    gsnv = GermlineSomaticNucleotideVariants(gsnv_raw)
                    get_neighbour_gsnv_info(gsnv, bam, out_alignment_records)

                    # Per snv positions, write records
            if len(out_alignment_records) > 0:
                out_alignment_records.write_records()

                count_bam_read += len(out_alignment_records)
                arr = np.array(out_alignment_records.ordered_alignments_tags)
                allele_ref_counts += np.count_nonzero(arr == 0)
                allele_alt_counts += np.count_nonzero(arr > 0)
                allele_not_matched += np.count_nonzero(arr == -2)
                allele_debug_should_be_zero += np.count_nonzero(arr == -1)

        else:
            # The SNV did not pass filtering criteria
            continue

    # Check each gsnv alignment record if they are in snv list, if so, add info to snv read.
    print(
        "filtered SNV:",
        count_filtered_ssnv,
        "overlapping reads:",
        count_bam_read,
        "Ref allele:",
        allele_ref_counts,
        "Alt allele:",
        allele_alt_counts,
    )
    print("allele not matched either", allele_not_matched)
    print(
        "Tumor fraction derived from matching molecules =Alt/(Ref+Alt) ratio:",
        round(allele_alt_counts / (allele_alt_counts + allele_ref_counts), 5),
    )
    print(
        "Error rate =Err/(Ref+Alt+err) ratio:",
        round(
            allele_not_matched
            / (allele_alt_counts + allele_ref_counts + allele_not_matched),
            5,
        ),
    )

    if phasing:
        print("associated gSNV:", len(out_alignment_records.associated_gsnv))
    bam.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b",
        "--input_bam_path",
        type=str,
        help="Input bam file of interest for molecule to tag.",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--somatic_vcf_path",
        type=str,
        help="Input path of somatic vcf file (gz + tabix).",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--germline_vcf_path",
        type=str,
        help="Input germline vcf file. If nothing is given, no phasing is done.",
        default=None,
    )
    parser.add_argument(
        "-o",
        "--out_bam_path",
        type=str,
        help="Folder where output bam files should be stored.",
    )
    parser.add_argument(
        "-p",
        "--phasing",
        type=bool,
        default=False,
        help="If phasing should be applied.",
    )
    parser.add_argument("-q", "--quiet", action="store_true", default=False)
    args = parser.parse_args()
    examine_path(args.out_bam_path)
    filename = os.path.basename(args.input_bam_path).split(".")[0]
    out_ref_bam_path = f"{args.out_bam_path}/{filename}_vcf_molecules_ref.bam"
    out_alt_bam_path = f"{args.out_bam_path}/{filename}_vcf_molecules_alt.bam"
    out_ERR_bam_path = f"{args.out_bam_path}/{filename}_vcf_molecules_not_matched.bam"
    # Execute
    start_time = time.time()
    find_ssnv(
        args.input_bam_path,
        args.somatic_vcf_path,
        args.germline_vcf_path,
        out_ref_bam_path,
        out_alt_bam_path,
        out_ERR_bam_path,
        phasing=False,
    )
    end_time = time.time()
    print(
        f"Program finishes successfully. Time spent: {round((end_time - start_time)/60,2)} mins."
    )
    print("input bam:", args.input_bam_path)
    print("output bam folder:", args.out_bam_path)
