from __future__ import print_function
import os
import signal
import argparse
import numpy as np

import logging
logging.basicConfig(
    format='%(message)s',
    level=os.getenv('LOG_LEVEL', logging.WARNING),
)

from . import IFOS, load_budget, plot_noise

##################################################

description = """Plot GWINC noise budget for specified IFO.

Available included IFOs: {}

""".format(', '.join(["'{}'".format(ifo) for ifo in IFOS]))
# for ifo in available_ifos():
#     description += "  '{}'\n".format(ifo)
description += """

By default a GWINC noise budget of the specified IFO will calculated,
and plotted with an interactive plotter.  If the --save option is
specified the plot will be saved directly to a file (without display)
(various formats are supported, indicated by file extension).  If the
requested extension is 'hdf5' or 'h5' then the noise traces and IFO
parameters will be saved to an HDF5 file (without plotting).  The
input file (IFO) can be an HDF5 file saved from a previous call, in
which case all noise traces and IFO parameters will be loaded from
that file.

Individual IFO parameters can be overriden with the --ifo option:

  gwinc --ifo Optics.SRM.Tunephase=3.14 ...

If the inspiral_range package is installed, various figures of merit
can be calculated for the resultant spectrum with the --fom option,
e.g.:

  gwinc --fom horizon ...
  gwinc --fom range:m1=20,m2=20 ...

See documentation for inspiral_range package for details.

"""

IFO = 'aLIGO'
FREQ = '5:3000:6000'
DATA_SAVE_FORMATS = ['.hdf5', '.h5']

parser = argparse.ArgumentParser(
    prog='gwinc',
    description=description,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    '--freq', '-f', metavar='FSPEC',
    help="frequency array specification in Hz, either as 'flo:fhi' or 'flo:npoints:fhi' [{}]".format(FREQ))
parser.add_argument(
    '--ifo', '-o',
    #nargs='+', action='extend',
    action='append',
    help="override budget IFO parameter (may be specified multiple times)")
parser.add_argument(
    '--title', '-t',
    help="plot title")
parser.add_argument(
    '--fom',
    help="calculate inspiral range for resultant spectrum ('func[:param=val,param=val]')")
group = parser.add_mutually_exclusive_group()
group.add_argument(
    '--interactive', '-i', action='store_true',
    help="interactive plot with interactive shell")
group.add_argument(
    '--save', '-s', action='append',
    help="save plot (.png/.pdf/.svg) or budget traces (.hdf5/.h5) to file (may be specified multiple times)")
group.add_argument(
    '--yaml', '-y', action='store_true',
    help="print IFO as yaml to stdout and exit (budget not calculated)")
group.add_argument(
    '--text', '-x', action='store_true',
    help="print IFO as text table to stdout and exit (budget not calculated)")
group.add_argument(
    '--diff', '-d', metavar='IFO',
    help="show differences table between another IFO description and exit (budget not calculated)")
group.add_argument(
    '--no-plot', '-np', action='store_false', dest='plot',
    help="supress plotting")
parser.add_argument(
    'IFO',
    help="IFO name, or path to budget module (.py), description file (.yaml/.mat/.m), or HDF5 data file (.hdf5/.h5)")


def freq_from_spec(spec):
    fspec = spec.split(':')
    if len(fspec) == 2:
        fspec = fspec[0], FREQ.split(':')[1], fspec[1]
    return np.logspace(
        np.log10(float(fspec[0])),
        np.log10(float(fspec[2])),
        int(fspec[1]),
    )


def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    args = parser.parse_args()

    ##########
    # initial arg processing

    if os.path.splitext(os.path.basename(args.IFO))[1] in DATA_SAVE_FORMATS:
        from .io import load_hdf5
        Budget = None
        freq, traces, attrs = load_hdf5(args.IFO)
        ifo = getattr(attrs, 'IFO', None)
        plot_style = attrs
        if args.freq:
            logging.warning("ignoring frequency specification for frequencies defined in HDF5...")

    else:
        Budget = load_budget(args.IFO)
        ifo = Budget.ifo
        if args.freq:
            try:
                freq = freq_from_spec(args.freq)
            except IndexError:
                parser.error("improper frequency specification '{}'".format(args.freq))
        else:
            freq = getattr(Budget, 'freq', freq_from_spec(FREQ))
        plot_style = getattr(Budget, 'plot_style', {})
        traces = None

    out_data_files = set()
    out_plot_files = set()
    if args.save:
        for path in args.save:
            if os.path.splitext(path)[1] in DATA_SAVE_FORMATS:
                out_data_files.add(path)
        out_plot_files = set(args.save) - out_data_files

    if args.ifo:
        for paramval in args.ifo:
            param, val = paramval.split('=', 1)
            ifo[param] = float(val)

    if args.yaml:
        if not ifo:
            parser.exit(2, "no IFO structure available.")
        print(ifo.to_yaml(), end='')
        return
    if args.text:
        if not ifo:
            parser.exit(2, "no IFO structure available.")
        print(ifo.to_txt(), end='')
        return
    if args.diff:
        if not ifo:
            parser.exit(2, "no IFO structure available.")
        fmt = '{:30} {:>20} {:>20}'
        Budget = load_budget(args.diff)
        diffs = ifo.diff(Budget.ifo)
        if diffs:
            print(fmt.format('', args.IFO, args.diff))
            for p in diffs:
                k = str(p[0])
                v = repr(p[1])
                ov = repr(p[2])
                print(fmt.format(k, v, ov))
        return

    if args.title:
        plot_style['title'] = args.title
    elif Budget:
        plot_style['title'] = "GWINC Noise Budget: {}".format(Budget.name)
    else:
        plot_style['title'] = "GWINC Noise Budget: {}".format(args.IFO)

    if args.plot:
        if args.save:
            # FIXME: this silliness seems to be the only way to have
            # matplotlib usable on systems without a display.  There must
            # be a better way.  'AGG' is a backend that works without
            # displays.  but it has to be set before any other matplotlib
            # stuff is imported.  and we *don't* want it set if we do want
            # to show an interactive plot.  there doesn't seem a way to
            # set this opportunistically.
            import matplotlib
            matplotlib.use('AGG')
        try:
            from matplotlib import pyplot as plt
        except RuntimeError:
            logging.warning("no display, plotting disabled.")
            args.plot = False

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

    if not traces:
        logging.info("calculating budget...")
        traces = Budget(freq=freq, ifo=ifo).run()

    # logging.info('recycling factor: {: >0.3f}'.format(ifo.gwinc.prfactor))
    # logging.info('BS power:         {: >0.3f} W'.format(ifo.gwinc.pbs))
    # logging.info('arm finesse:      {: >0.3f}'.format(ifo.gwinc.finesse))
    # logging.info('arm power:        {: >0.3f} kW'.format(ifo.gwinc.parm/1000))

    if args.fom:
        logging.info("calculating inspiral {}...".format(range_func))
        H = inspiral_range.CBCWaveform(freq, **range_params)
        logging.debug("params: {}".format(H.params))
        fom = eval('inspiral_range.{}'.format(range_func))(freq, traces['Total'][0], H=H)
        logging.info("{}({}) = {:.2f} Mpc".format(range_func, fargs, fom))
        fom_title = '{func} {m1}/{m2} Msol: {fom:.2f} Mpc'.format(
            func=range_func,
            m1=H.params['m1'],
            m2=H.params['m2'],
            fom=fom,
            )
        plot_style['title'] += '\n{}'.format(fom_title)

    ##########
    # interactive

    # interactive shell plotting
    if args.interactive:
        from IPython.terminal.embed import InteractiveShellEmbed
        ipshell = InteractiveShellEmbed(
            banner1='''
GWINC interactive plotter

You may interact with plot using "plt." methods, e.g.:

>>> plt.title("foo")
>>> plt.savefig("foo.pdf")
''',
            user_ns={
                'freq': freq,
                'traces': traces,
                'ifo': ifo,
                'plot_style': plot_style,
                'plot_noise': plot_noise,
            },
        )
        ipshell.enable_pylab()
        ipshell.ex("fig = plot_noise(freq, traces, **plot_style)")
        ipshell.ex("plt.title('{}')".format(plot_style['title']))
        ipshell()

    ##########
    # output

    # save noise traces to HDF5 file
    if out_data_files:
        from .io import save_hdf5
        attrs = dict(plot_style)
        attrs['IFO'] = ifo.to_yaml()
        for path in out_data_files:
            logging.info("saving budget traces {}...".format(path))
            save_hdf5(
                path=path,
                freq=freq,
                traces=traces,
                **attrs,
            )

    # standard plotting
    if args.plot or out_plot_files:
        logging.info("plotting noises...")
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plot_noise(
            freq,
            traces,
            ax=ax,
            **plot_style
        )
        fig.tight_layout()
        if out_plot_files:
            for path in out_plot_files:
                fig.savefig(path)
        else:
            plt.show()


if __name__ == '__main__':
    main()
