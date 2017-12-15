import signal
import argparse
import numpy as np
from matplotlib import pyplot as plt
from IPython.terminal.embed import InteractiveShellEmbed

import logging
logging.basicConfig(level=logging.INFO)

from . import load_ifo
from . import plot_noise

##################################################

description = """Plot GWINC noise budget for specified IFO

If the inspiral_range package is installed, various figures of merit
can be calculated for the resultant spectrum with the --fom argument,
e.g.:

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
                    help="calculate inspiral range for resultant spectrum ('func:param=val,param=val')")
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

    freq = np.logspace(np.log10(args.flo), np.log10(args.fhi), args.npoints)

    if args.matlab:
        from .gwinc_matlab import gwinc_matlab as gwinc
    # FIXME: why weird import issue requires doing this here instead
    # of above?
    else:
        from . import gwinc

    if args.fom:
        import inspiral_range
        logging.info("calculating inspiral range...")
        try:
            ffunc, fargs = args.fom.split(':')
        except ValueError:
            ffunc = args.fom
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
        range_params = inspiral_range.waveform._get_waveform_params(**range_params)
        range_func = eval('inspiral_range.{}'.format(ffunc))

    score, noises, ifo = gwinc(freq, ifo)

    if args.fom:
        psd = (freq, noises['Total'])
        fom = range_func(psd, **range_params)
        fom_title = '{func} {m1}/{m2}: {fom:.3f} Mpc'.format(
            func=range_func.__name__,
            m1=range_params['m1'],
            m2=range_params['m2'],
            fom=fom,
            )
        if args.title:
            args.title += '\n{}'.format(fom_title)
        else:
            args.title = fom_title

    if args.interactive:
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
        ipshell.run_code("plot_noise(noises)")
        if args.title:
            ipshell.run_code("plt.title('{}')".format(args.title))
        ipshell()
    else:
        plot_noise(noises)
        if args.title:
            plt.title(args.title)
        if args.save:
            plt.savefig(args.save)
        else:
            plt.show()


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    main()
