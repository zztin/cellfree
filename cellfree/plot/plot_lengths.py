import collections
import numpy as np
import pysam
from matplotlib import pyplot as plt
from scipy.signal import find_peaks


def plot_forward_length(length_counter, png, dpi=100, start=30, end=1000):
    zip_rl = sorted(length_counter.items())
    plt.plot([x for x in range(len(zip_rl))], [y for y in zip_rl], label="illumina", color='b')



    plt.legend()
    plt.savefig(png, dpi=dpi)
    plt.show()
