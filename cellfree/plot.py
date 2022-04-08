import argparse
import collections

import numpy as np
import pysam
from matplotlib import pyplot as plt
from scipy.signal import find_peaks


def plot_peaks_prominence_compare(
    rl_lists,
    filename,
    namelist,
    xmin=30,
    xmax=1000,
    prominence_param=20,
    colors=["green", "orange", "blue", "violet", "red"],
    density=True,
    tickstep=10,
    vertical_lines="top3",
    outpath="./",
):
    """
    rl_lists: rl_lists are generated from the plot_length_multi

    """
    plt.subplots(figsize=(24, 8))
    cont_xy_list = []
    sum_y_per_bam = []
    max_y_per_bam = []

    for k, rl_list in enumerate(rl_lists):
        # construct continuous x axis and corresponding y value
        rl_count = collections.Counter(rl_list)
        zip_rl = sorted(rl_count.items())
        x = [x for (x, y) in zip_rl if x >= 30]
        name = namelist[k]
        print(f"{name} max mapped length = {max(x)}")
        y = [y for (x, y) in zip_rl if x >= 30]
        cont_x = range(xmin, max(x))
        if xmax:
            cont_x = cont_x[: xmax - xmin + 1]
        cont_y = []
        i = 0
        for xi in cont_x:
            if xi not in x:
                cont_y.append(0)
            else:
                index = x.index(xi)
                cont_y.append(y[index])
                assert xi == x[index]
                i += 1
        cont_xy_list.append((cont_x, cont_y))
        sum_y_per_bam.append(sum(cont_y))
        max_y_per_bam.append(max(cont_y))
        # print(cont_x[-1])
    max_y = max(max_y_per_bam)
    max_sum_y = max(sum_y_per_bam)
    # plot actual

    # plot density

    for k, (cont_x, cont_y) in enumerate(cont_xy_list):
        name = namelist[k]
        color = colors[k]
        sum_y_per = sum_y_per_bam[k]
        max_y_per = max_y_per_bam[k]
        data = np.array(cont_y)
        peaks = find_peaks(data, prominence=prominence_param)
        # largest 3 peaks:
        #        max_y_index = data.index(max(peaks[0]))
        #        print(f"{name} Peaks:", peaks[0][:30])
        if density:
            density_data = data / sum_y_per  # *(max_sum_y/max_y )
            # print("density")
            # print(max(data))
        else:
            density_data = data
            # print("not density")
            # print(max(data))
        plt.plot(cont_x, density_data, label=name, color=color)

        plt.plot(
            peaks[0] + xmin,
            density_data[peaks[0]],
            "x",
            label=f"{name}: {'bp, '.join([str(x + xmin) for x in peaks[0][:10]])}bp",
            color=f"dark{color}",
        )
        if vertical_lines == "all":
            for peak in peaks[0] + xmin:
                plt.axvline(peak, color=f"dark{color}", alpha=0.1)
        elif vertical_lines == "top3":
            ind = np.argpartition(density_data[peaks[0]], -3)[-3:]
            print("top3", density_data[peaks[0]][ind])
            for peak in (peaks[0] + xmin)[ind]:
                plt.axvline(peak, color=f"dark{color}", alpha=0.1)

    plt.title(filename)
    plt.xticks(np.arange(xmin, xmax + tickstep, tickstep))
    plt.legend(fontsize=10)
    plt.savefig(
        f"{outpath}/{filename}_{prominence_param}_{xmin}_{xmax}_density{density}_prominence_peaks_xmax{xmax}.png",
        dpi=80,
    )
    # plt.show()

    return cont_xy_list


def plot_length_multi(
    bamfile_list,
    filename,
    outpath="../",
    density=False,
    len_range=[10, 1500],
    histtype="step",
    alpha=0.7,
):
    # parse length from cigar string
    rl_list = []
    histogram_list = []
    fig, ax = plt.subplots(figsize=(12, 4))
    for name, file in bamfile_list:
        with pysam.AlignmentFile(file) as f:
            read_lengths = []
            for read in f:
                rn = 0
                if read.cigartuples:
                    for i, j in read.cigartuples:
                        if i == 0:
                            rn += j
                read_lengths.append(rn)
        rl_list.append(read_lengths)

        histogram = plt.hist(
            read_lengths,
            range=len_range,
            label=name,
            bins=int((len_range[1] - len_range[0]) / 2),
            alpha=alpha,
            density=density,
            histtype=histtype,
        )
        histogram_list.append(histogram)
    if density:
        plt.ylabel("read ratio")
    else:
        plt.ylabel("read count")
    # Mapped length
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    plt.xlabel("fragment length (bp)")
    plt.xticks(range(len_range[0], len_range[1], 50))
    plt.legend()
    plt.title(filename)
    plt.savefig(f"{outpath}/{filename}_read_length_max{len_range[1]}.png", dpi=250)
    # plt.show()
    return rl_list, histogram_list


parser = argparse.ArgumentParser(description="Plot the fragment length comparison")
parser.add_argument("filepaths", nargs="+", type=str, help="profile file path to plot")
parser.add_argument("--out_path", type=str, default="./", help="output path")
parser.add_argument("--tickstep", type=int, default=100, help="tickstep")
parser.add_argument("--xmin", type=int, default=10, help="min plotting range")
parser.add_argument("--xmax", type=int, default=800, help="max plotting range")


args = parser.parse_args()
print(args.filepaths)
filenames = [x.split("/")[-1].split("_")[0] for x in args.filepaths]
plotnames = [x.split("/")[-1].split(".")[0] for x in args.filepaths]
bamfilelist = list(zip(filenames, args.filepaths))
print(bamfilelist)
# file2 = "/Users/lchen/current_projects/concall-working/output/DER4387/DER4387_tide_fl.sorted.bam"
# file4 ="/Users/lchen/current_projects/concall-working/results/bam_to_tag/MAR5937_tide_fl_no_bb.sorted.bam"
# file5 = "/Users/lchen/current_projects/concall-working/results/circular/DER5088_th_only_sing_tide_no_bb_double.sorted.bam"

# bamfilelist = [("Control", file2), ("EAC_OES_1033-Flongle",file4)]
rl_list2, histogram2 = plot_length_multi(
    bamfilelist,
    f"{''.join(plotnames)}",
    outpath=args.out_path,
    histtype="bar",
    density=True,
    alpha=0.2,
    len_range=[args.xmin, args.xmax],
)

cont_xy = plot_peaks_prominence_compare(
    rl_list2,
    f"{'_'.join(filenames)}",
    plotnames,
    xmin=args.xmin,
    xmax=args.xmax,
    tickstep=args.tickstep,
    outpath=args.out_path,
)
cont_xy = plot_peaks_prominence_compare(
    rl_list2,
    f"{'_'.join(filenames)}",
    plotnames,
    density=False,
    xmin=args.xmin,
    xmax=args.xmax,
    tickstep=args.tickstep,
    outpath=args.out_path,
)
