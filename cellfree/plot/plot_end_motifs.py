import matplotlib.pyplot as plt


def plot_kmer_dist(
    out_path,
    motifs,
    ind,
    colors,  # =['blue', "red", "teal", "yellow"],
    relative=True,
    figsize=(10, 4),
):
    """
    motifs: [[AATT, AATC,...]]
    ind: name of the input list
    >>> plot_kmer_dist([df['fwd5'], df['fwd3']], ["fwd 5' end", "fwd 3' end"], [ "red", "blue"], path, figsize = figsize)
    """
    plt.figure(figsize=figsize)
    for i, x in enumerate(motifs):
        if relative:
            sum_total = sum(x.values())
        else:
            sum_total = 1
        plt.bar(
            [k for k, v in sorted(x.items())],
            [v / sum_total for k, v in sorted(x.items())],
            label=f"{ind[i]}",
            alpha=0.3,
            color=colors[i],
        )
    plt.xticks(rotation="vertical")
    plt.legend()
    plt.savefig(f"{out_path}", dpi=100)
