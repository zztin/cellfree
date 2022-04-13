import collections

import numpy as np
import pysam
from matplotlib import pyplot as plt
from scipy.signal import find_peaks


def plot_forward_length(length_counter, png, sample_name, dpi=100, start=0, end=1000):
    zip_rl = sorted(length_counter.items())
    plt.plot(
        [x for (x, y) in zip_rl], [y for (x, y) in zip_rl], label=sample_name, color="b"
    )
    plt.xlim(start, end)
    plt.xticks(range(start, end + 10, 100))
    plt.legend()
    plt.title("forward strand cfdna length")
    plt.xlabel("cfdna length (bp)")
    plt.savefig(f"{png}.forward.png", dpi=dpi)
    plt.cla()


def plot_reverse_length(length_counter, png, sample_name, dpi=100, start=0, end=1000):
    zip_rl = sorted(length_counter.items())
    plt.cla()
    plt.plot(
        [x for (x, y) in zip_rl],
        [y for (x, y) in zip_rl],
        label=sample_name,
        color="black",
    )
    plt.xlim(start, end)
    plt.xticks(range(start, end + 10, 100))
    plt.legend()
    plt.title("reverse strand cfdna length")
    plt.xlabel("cfdna length (bp)")
    plt.savefig(f"{png}.reverse.png", dpi=dpi)
    plt.cla()
