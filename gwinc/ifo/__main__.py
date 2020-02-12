import argparse
import numpy as np
import matplotlib.pyplot as plt

from . import IFOS, PLOT_STYLE
from .. import load_budget


FLO = 3
FHI = 10000
NPOINTS = 3000
YLIM = (1e-25, 1e-20)


def main():
    parser = argparse.ArgumentParser(
        description="Reference IFO comparison plot",
    )
    parser.add_argument(
        '--save', '-s',
        help="save plot to file (.pdf/.png/.svg)")
    args = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    freq = np.logspace(np.log10(FLO), np.log10(FHI), NPOINTS)

    for ifo in IFOS:
        Budget = load_budget(ifo)
        data = Budget(freq).calc()
        label = Budget.name

        ax.loglog(freq, np.sqrt(data), label=label)

    ax.grid(
        True,
        which='both',
        lw=0.5,
        ls='-',
        alpha=0.5,
    )

    ax.legend(
        # ncol=2,
        fontsize='small',
    )

    ax.autoscale(enable=True, axis='y', tight=True)
    ax.set_ylim(*YLIM)
    ax.set_xlim(freq[0], freq[-1])

    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(PLOT_STYLE['ylabel'])
    ax.set_title("PyGWINC reference IFO strain comparison")

    if args.save:
        fig.savefig(args.save)
    else:
        plt.show()


if __name__ == '__main__':
    main()
