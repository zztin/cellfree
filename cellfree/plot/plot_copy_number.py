# Create CNV/CNA copy number alteration plots for multiple samples. CNA and tumor fraction must be called beforehand.
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from singlecellmultiomics.utils import get_contig_list_from_fasta, is_main_chromosome


# Define chromsome order:
def sort_chromosome_names(l):
    chrom_values = []
    for chrom in l:
        chrom_value = None
        chrom = chrom.replace("chr", "").upper()
        if chrom == "X":
            chrom_value = 99
        elif chrom == "Y":
            chrom_value = 100
        elif chrom == "M" or chrom == "MT":
            chrom_value = 101
        elif chrom == "EBV":
            chrom_value = 102
        elif chrom == "MISC_ALT_CONTIGS_SCMO":
            chrom_value = 999
        else:
            try:
                chrom_value = int(chrom)
            except Exception as e:
                chrom_value = 999 + sum((ord(x) for x in chrom))
        chrom_values.append(chrom_value)

    indices = sorted(range(len(chrom_values)), key=lambda x: chrom_values[x])
    return [l[idx] for idx in indices]


def contigs_conversion_t2t(x):
    contigs_conversion_dict = {
        "NC_060925.1": "chr1",
        "NC_060926.1": "chr2",
        "NC_060927.1": "chr3",
        "NC_060928.1": "chr4",
        "NC_060929.1": "chr5",
        "NC_060930.1": "chr6",
        "NC_060931.1": "chr7",
        "NC_060932.1": "chr8",
        "NC_060933.1": "chr9",
        "NC_060934.1": "chr10",
        "NC_060935.1": "chr11",
        "NC_060936.1": "chr12",
        "NC_060937.1": "chr13",
        "NC_060938.1": "chr14",
        "NC_060939.1": "chr15",
        "NC_060940.1": "chr16",
        "NC_060941.1": "chr17",
        "NC_060942.1": "chr18",
        "NC_060943.1": "chr19",
        "NC_060944.1": "chr20",
        "NC_060945.1": "chr21",
        "NC_060946.1": "chr22",
        "NC_060947.1": "chrX",
        "NC_060948.1": "chrY",
    }
    return contigs_conversion_dict[x]


class GenomicPlot:
    def __init__(self, ref_path, contigs=None, ignore_contigs=None, cluster=False):
        """
        Initialise genomic plot

        ref_path(str or pysam.FastaFile) : Path or handle to reference

        """

        if contigs is None:
            self.contigs = sort_chromosome_names(
                list(
                    filter(
                        lambda x: is_main_chromosome(x)
                        and (ignore_contigs is None or x not in ignore_contigs),
                        get_contig_list_from_fasta(ref_path),
                    )
                )
            )
        else:
            self.contigs = contigs

        # Obtain the lengths:
        if type(ref_path) is str:
            with pysam.FastaFile(ref_path) as reference:
                self.lengths = {
                    r: l
                    for r, l in zip(reference.references, reference.lengths)
                    if r in self.contigs
                }
        else:
            self.lengths = {
                r: l
                for r, l in zip(ref_path.references, ref_path.lengths)
                if r in self.contigs
            }

        self.total_bp = sum(self.lengths.values())
        # Prune contigs with no length:
        self.contigs = [contig for contig in self.contigs if contig in self.lengths]
        # self.contigs = [contigs_conversion_t2t(contig) for contig in self.contigs]
        # self.lengths = {contigs_conversion_t2t(r): l for (r, l) in self.lengths.items()}
        self.contigs = sort_chromosome_names(self.contigs)
        self.cluster = cluster

    def cn_heatmap(
        self,
        df,
        cell_font_size=3,
        max_cn=4,
        method="ward",
        cmap="bwr",
        yticklabels=True,
        figsize=(15, 20),
        xlabel="Chromosomes",
        ylabel="Samples",
        vmin=0,
        xtickfontsize=8,
        cluster=False,
        **kwargs,
    ):
        """
        Create a heatmap from a copy number matrix

        df: triple indexed dataframe with as columns ('contig', start, end ), as rows cells/samples

        cell_font_size (int): font size of the cell labels

        max_cn (int) : dataframe will be clipped to this value. (Maximum copy number shown)

        method (str) : clustering metric

        cmap (str) : colormap used

        figsize(tuple) : Size of the figure

        xlabel (str) : Label for the x-axis, by default this is Contigs

        ylabel (str) : Label for the x-axis, by default this is Cells

        **kwargs : Arguments which will be passed to seaborn.clustermap

        """

        allelic_mode = len(df.columns[0]) == 4
        if allelic_mode:
            alleles = [
                allele
                for allele in df.columns.get_level_values(0).unique()
                if not pd.isna(allele)
            ]
            contigs_to_plot = [
                contig
                for contig in self.contigs
                if contig in set(df.columns.get_level_values(1))
            ]
            # Resample the dataframe, drop columns with no allele assigned:
            df = df.loc[:, df.columns.isin(contigs_to_plot, level=1)][
                alleles
            ].sort_index(1)

            def m(k):
                allele, contig, start, end = k
                return self.contigs.index(contig), alleles.index(allele), start

            desired_order = sorted(
                list(
                    df.loc[:, df.columns.isin(self.contigs, level=1)][alleles]
                    .sort_index(1)
                    .columns
                ),
                key=m,
            )
            df = df[desired_order]

        else:

            # Figure out what contigs are present in the dataframe:
            contigs_to_plot = [
                contig
                for contig in self.contigs
                if contig in set(df.columns.get_level_values(0))
            ]
            df = df.sort_index(1)[contigs_to_plot]
        if cluster:
            try:
                clmap = sns.clustermap(
                    df,
                    col_cluster=False,
                    method=method,
                    cmap=cmap,
                    vmax=max_cn,
                    vmin=0,
                    yticklabels=yticklabels,
                    figsize=figsize,
                    **kwargs,
                )
                ax_heatmap = clmap.ax_heatmap
            except Exception as e:
                print(e)
                print("Falling back on heatmap without clustering")

                fig, ax_heatmap = plt.subplots(figsize=figsize)
                clmap = sns.heatmap(
                    df,
                    cmap=cmap,
                    vmax=max_cn,
                    vmin=vmin,
                    yticklabels=True,
                    ax=ax_heatmap,
                    **kwargs,
                )
        else:
            fig, ax_heatmap = plt.subplots(figsize=figsize)
            clmap = sns.heatmap(
                df,
                cmap=cmap,
                vmax=max_cn,
                vmin=vmin,
                yticklabels=True,
                ax=ax_heatmap,
                **kwargs,
            )
        prev = None
        xtick_pos = []
        xtick_label = []
        last_idx = 0

        allele = None
        for idx, key in enumerate(df.columns):
            if allelic_mode:
                (allele, contig, start, end) = key
            else:
                (contig, start, end) = key

            # Clean up contig label:
            contig = contig.replace("chr", "")
            if allele is not None:
                contig = f"{contig}:{allele}"
            if prev is not None and prev != contig:
                ax_heatmap.axvline(idx - 0.5, c="k", lw=1.5, zorder=10)
                xtick_pos.append((idx + last_idx) / 2)
                xtick_label.append(prev)
                last_idx = idx
            prev = contig

        # Plot last tick..
        xtick_pos.append((idx + last_idx) / 2)
        xtick_label.append(contig)

        ax_heatmap.set_xticks(xtick_pos)
        ax_heatmap.set_xticklabels(xtick_label, rotation=0, fontsize=xtickfontsize)
        ax_heatmap.set_xlabel(xlabel, labelpad=20)
        ax_heatmap.set_ylabel(ylabel, labelpad=20)

        return clmap, (ylabel, xtick_pos, xtick_label)

    def get_relative_widths(self):
        return [self.lengths[contig] / self.total_bp for contig in self.contigs]

    def reset_axis(self, contig):
        ax = self[contig]
        ax.clear()
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlabel(contig.replace("chr", ""))
        ax.set_xlim(0, self.lengths[contig])

    def get_figure(self, figsize=(20, 1)):

        widths = self.get_relative_widths()

        gs_kw = dict(width_ratios=widths)
        figure = plt.figure(figsize=figsize)
        figure.subplots_adjust(bottom=0.25, top=0.75)

        self.gridspec = gridspec.GridSpec(
            1, len(widths), figure=figure, wspace=0.1, width_ratios=widths
        )
        self.axis = {}
        prev_ax = None
        for i, contig in enumerate(self.contigs):
            # i = i + 1 # grid spec indexes from 0

            ax = plt.subplot(self.gridspec[i], sharey=prev_ax)
            self.axis[contig] = ax
            self.reset_axis(contig)
            prev_ax = ax
        sns.despine(left=True)
        figure.canvas.draw()
        return figure

    def __getitem__(self, contig):
        return self.axis[contig]


def get_tumor_fraction_from_ichorCNA(df):

    df["Tf", "tf_1", "tf_2"] = [0.19, 0.09, 0.09, 0.08, 0.52, 0.52, 0.12, 0.65, 0.32]
    new_df = df.apply(lambda row: (row - 2) * row["Tf", "tf_1", "tf_2"] + 2, axis=1)
    new_df.drop(columns=["Tf"], inplace=True)
    return new_df


def remove_noisy_bins(df):
    df = df.fillna(2)
    pre = df.diff()
    post = df.diff(periods=-1)
    subtract = pre.subtract(post)
    df2 = df.copy()
    ## diff != 0  , subtract == 0
    for column in df.columns:
        indx1 = np.where(pre[column] != 0)
        indx2 = np.where(subtract[column] == 0)
        indx_inter = set(indx1[0]).intersection(set(indx2[0]))
        df2[column][list(indx_inter)] = (
            df2[column][list(indx_inter)] - pre[column][list(indx_inter)]
        )
    return df2


if __name__ == "__main__":
    input_dir = "/Users/liting/01_data/AMC_ascites/02_CNA_compare/CNA/"
    logr_copy_number = []
    ccnas = []
    files = [
        "HC01_MAR6510_b24_dedup.cna.seg",
        "HC02_MAR6510_b23_dedup.cna.seg",
        "HC03_MAR6510_b22_dedup.cna.seg",
        "HC04_DER4387_raw.cna.seg",
        "DER5087.cna.seg",
        "MAR6939.cna.seg",
        "MAR7069.cna.seg",
        "MAR7195.cna.seg",
        "OVCA06_pass_dedup.cna.seg",
    ]  # os.listdir(input_dir):
    names = [
        "HC01",
        "HC02",
        "HC03",
        "HC04",
        "OVCA01",
        "OVCA03",
        "OVCA04",
        "OVCA05",
        "OVCA06",
    ]
    for file, plot_label in zip(files, names):
        column_name = file.split(".")[0]

        path = input_dir + file
        cna = pd.read_csv(path, sep="\t")
        cna["contig"] = cna.apply(
            lambda row: (row["chr"], row["start"], row["end"]), axis=1
        )
        index = pd.MultiIndex.from_tuples(cna["contig"])
        x = cna.set_index(index)
        # subset = x[[f'{column_name}.logR_Copy_Number']]
        # subset.columns= [plot_label]
        ccna = x[[f"{column_name}.Corrected_Copy_Number"]]
        ccna.columns = [plot_label]
        ccnas.append(ccna)
    #     logr_copy_number.append(subset)
    # logr_copy_number_df = pd.concat(logr_copy_number, axis = 1).T
    ccnas_df = pd.concat(ccnas, axis=1)
    # Remove CNV change that only exist in 1 MB bin, not neighbouring bins. Fill with the value before the bin.
    # ccnas_df = remove_noisy_bins(ccnas_df)
    # Get tumor fraction
    new_ccnas_df = get_tumor_fraction_from_ichorCNA(ccnas_df.T)
    # ccnas_df.fillna(2.0) # If values are filled, the rows can be clustered. Otherwise, showed in the order of the input files.
    p = GenomicPlot(
        "/Users/liting/01_data/GENOMES/GRCh37_decoy/references_hs37d5_hs37d5.fa"
    )
    x = p.cn_heatmap(new_ccnas_df, max_cn=4, figsize=(15, 6), cluster=False)
    plt.savefig(
        "/Users/liting/Documents/02SurfDrive_Sync/PhD/00_Projects/02_GW_cyclomics/04_FIGURES/CNA_of_samples_adjusted_by_tf.png",
        dpi=300,
    )
