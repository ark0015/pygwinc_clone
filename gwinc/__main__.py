from __future__ import print_function
import os
import signal
import argparse
import numpy as np
from matplotlib import pyplot as plt
from IPython.terminal.embed import InteractiveShellEmbed

import logging
logging.basicConfig(format='%(message)s',
                    level=os.getenv('LOG_LEVEL', logging.INFO))

from .ifo import available_ifos, load_ifo
from .precomp import precompIFO
from .gwinc import gwinc as pygwinc
from . import gwinc_matlab
from . import plot_noise
from . import util

##################################################

description = """Plot GWINC noise budget for specified IFO.

Available included IFOs: {}

""".format(', '.join(["'{}'".format(ifo) for ifo in available_ifos()]))
# for ifo in available_ifos():
#     description += "  '{}'\n".format(ifo)
description += """
If the inspiral_range package is installed, various figures of merit
can be calculated for the resultant spectrum with the --fom argument,
e.g.:

  gwinc --fom range ...
  gwinc --fom range:m1=20,m2=20 ...

See documentation for inspiral_range package for details.

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
parser.add_argument('--npoints', '-n', default=NPOINTS,
                    help="number of frequency points [{}]".format(NPOINTS))
parser.add_argument('--title', '-t',
                    help="plot title")
parser.add_argument('--matlab', '-m', action='store_true',
                    help="use MATLAB gwinc engine to calculate noises")
parser.add_argument('--fom',
                    help="calculate inspiral range for resultant spectrum ('func[:param=val,param=val]')")
parser.add_argument('--displacement', '-D', action='store_true', default=False, dest='displacement',
                   help="supress adding displacement sensitivity axis")
group = parser.add_mutually_exclusive_group()
group.add_argument('--dump', '-d', dest='dump', action='store_true',
                   help="print IFO parameters to stdout and exit")
group.add_argument('--interactive', '-i', action='store_true',
                   help="open interactive shell when plotting")
group.add_argument('--save', '-s',
                   help="save noise (hdf5) or figure (pdf/png/etc) to file")
group.add_argument('--no-plot', '-np', action='store_false', dest='plot',
                   help="supress plotting")
parser.add_argument('IFO', default=IFO,
                    help="IFO name or description file path (.yaml, .mat, .m), or full HDF5 data file")


def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    args = parser.parse_args()

    ##########
    # initial arg processing

    if os.path.splitext(args.IFO)[1] == '.hdf5':
        title, ifo, noises = util.load_hdf5(args.IFO)

    else:
        ifo = load_ifo(args.IFO)
        noises = None
        if args.title:
            title = args.title
        else:
            title = '{} GWINC Noise Budget'.format(args.IFO)

    if args.dump:
        ifo = precompIFO(ifo)
        print(ifo.to_txt(), end='')
        return

    if args.matlab:
        gwinc = gwinc_matlab.gwinc_matlab
    else:
        gwinc = pygwinc

    if args.fom:
        import inspiral_range
        try:
            range_func, fargs = args.fom.split(':')
        except ValueError:
            range_func = args.fom
            fargs = ''
        range_params = {}
        for param in fargs.split(','):
            if not param:
                continue
            p,v = param.split('=')
            if not v:
                raise ValueError('missing parameter value "{}"'.format(p))
            if p != 'approximant':
                v = float(v)
            range_params[p] = v

    ##########
    # main calculations

    if noises:
        freq = noises['Freq']

    else:
        logging.info("calculating noises...")
        freq = np.logspace(np.log10(args.flo), np.log10(args.fhi), args.npoints)
        score, noises, ifo = gwinc(freq, ifo)

    logging.info('recycling factor: {: >0.3f}'.format(ifo.gwinc.prfactor))
    logging.info('BS power:         {: >0.3f} W'.format(ifo.gwinc.pbs))
    logging.info('arm finesse:      {: >0.3f}'.format(ifo.gwinc.finesse))
    logging.info('arm power:        {: >0.3f} kW'.format(ifo.gwinc.parm/1000))

    if args.fom:
        logging.info("calculating inspiral {}...".format(range_func))
        H = inspiral_range.CBCWaveform(freq, **range_params)
        logging.debug("params: {}".format(H.params))
        fom = eval('inspiral_range.{}'.format(range_func))(freq, noises['Total'], H=H)
        logging.info("{}({}) = {:.2f} Mpc".format(range_func, fargs, fom))
        fom_title = '{func} {m1}/{m2} Msol: {fom:.2f} Mpc'.format(
            func=range_func,
            m1=H.params['m1'],
            m2=H.params['m2'],
            fom=fom,
            )
        title += '\n{}'.format(fom_title)

    ##########
    # output

    # save noise traces to HDF5 file
    if args.save and os.path.splitext(args.save)[1] == '.hdf5':
        util.save_hdf5(
            title,
            ifo,
            noises,
            args.save,
        )

    # interactive shell plotting
    elif args.interactive:
        ipshell = InteractiveShellEmbed(
            user_ns={'freq': freq,
                     'noises': noises,
                     'ifo': ifo,
                     'plot_noise': plot_noise,
            },
            banner1='''
PYGWINC interactive plotter

You may interact with plot using "plt." methods, e.g.:

>>> plt.title("foo")
>>> plt.savefig("foo.pdf")
''')
        ipshell.enable_pylab()
        ipshell.run_code("plot_noise(ifo, noises)")
        ipshell.run_code("plt.title('{}')".format(title))
        ipshell()

    # standard plotting
    elif args.plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plot_noise(
            ifo,
            noises,
            ax=ax,
            displacement=args.displacement,
        )
        ax.set_title(title)
        fig.tight_layout()
        if args.save:
            fig.savefig(
                args.save,
            )
        else:
            plt.show()


if __name__ == '__main__':
    main()
