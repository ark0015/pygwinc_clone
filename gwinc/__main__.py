import signal
import argparse
import numpy as np
from matplotlib import pyplot as plt
from IPython.terminal.embed import InteractiveShellEmbed

import logging
logging.basicConfig(level=logging.INFO)

from . import load_ifo
from . import gwinc

##################################################

description = """Plot GWINC noise budget for specified IFO
"""

IFO = 'aLIGO'
FLO = 5
FHI = 6000
NPOINTS = 3000

parser = argparse.ArgumentParser(prog='gwinc',
                                 description=description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--flo', '-fl', default=FLO,
                    help="lower frequency bound in Hz [{}]".format(FLO))
parser.add_argument('--fhi', '--fh', default=FHI,
                    help="upper frequency bound in Hz [{}]".format(FHI))
parser.add_argument('--title', '-t',
                    help="plot title")
group = parser.add_mutually_exclusive_group()
group.add_argument('--dump', '-d', dest='dump', action='store_true',
                   help="print IFO parameters to stdout and exit")
group.add_argument('--save', '-s',
                   help="save figure to file")
group.add_argument('--interactive', '-i', action='store_true',
                   help="open interactive shell when plotting")
parser.add_argument('IFO', nargs='?', default=IFO,
                    help="IFO name or description file path (.yaml or .mat)")


def main():
    args = parser.parse_args()

    ifo = load_ifo(args.IFO)

    if args.dump:
        for k,v in sorted(ifo.walk()):
            print('{:50} {}'.format(k,v))
        return

    if args.interactive:
        ipshell = InteractiveShellEmbed(
            user_ns={'gwinc': gwinc,
                     'flo': args.flo,
                     'fhi': args.fhi,
                     'ifo': ifo,
            },
            banner1='''
PYGWINC interactive plotter

You may interact with plot using "plt." methods, e.g.:

>>> plt.title("foo")
>>> plt.savefig("foo.pdf")
''')
        ipshell.enable_pylab()
        ipshell.run_code("output = gwinc(flo, fhi, ifo, fig=True)")
        if args.title:
            ipshell.run_code("plt.title('{}')".format(args.title))
        ipshell()
    else:
        gwinc(args.flo, args.fhi, ifo, fig=True)
        if args.title:
            plt.title(args.title)
        if args.save:
            plt.savefig(args.save)
        else:
            plt.show()


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    main()
