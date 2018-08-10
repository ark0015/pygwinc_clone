from numpy import sqrt
import matplotlib.pyplot as plt

import logging

PRIORITY_MAP = {
    'Quantum Vacuum': 0,
    'Seismic': 1,
    'Newtonian Gravity': 2,
    'Suspension Thermal': 3,
    'Coating Brownian': 4,
    'Coating Thermo-Optic': 5,
    'ITM Thermo-Refractive': 6,
    'ITM Carrier Density': 7,
    'Substrate Brownian': 8,
    'Substrate Thermo-Elastic': 9,
    'Excess Gas': 10,
}

STYLE_MAP = {
    'Quantum Vacuum': dict(
        color = 'xkcd:vibrant purple',
    ),
    'Seismic': dict(
        color = 'xkcd:brown',
    ),
    'Newtonian Gravity': dict(
        color = 'xkcd:green',
    ),
    'Suspension Thermal': dict(
        color = 'xkcd:deep sky blue',
    ),
    'Coating Brownian': dict(
        color = 'xkcd:fire engine red',
    ),
    'Coating Thermo-Optic': dict(
        color = 'xkcd:bright sky blue',
        ls = '--',
    ),
    'ITM Thermo-Refractive': dict(
        color = 'xkcd:dark sky blue',
        ls = '--',
    ),
    'ITM Carrier Density': dict(
        color = 'xkcd:grey',
        ls = '--',
    ),
    'Substrate Brownian': dict(
        color = 'xkcd:pumpkin orange',
        ls = '--',
    ),
    'Substrate Thermo-Elastic': dict(
        color = 'xkcd:golden',
        ls = '--',
    ),
    'Excess Gas': dict(
        color = 'xkcd:baby shit brown',
        ls = '--',
    ),
}


def plot_noise(
        ifo,
        noises,
        ax = None,
        displacement = True,
):
    f = noises['Freq']

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    if displacement:
        ax_d = ax.twinx()
        ax_d.set_yscale('log')

        def convert_ax_h_to_d(ax):
            """
            Update second axis according with first axis.
            """
            y1, y2 = ax.get_ylim()
            ax_d.set_ylim(y1 * ifo.Infrastructure.Length, y2 * ifo.Infrastructure.Length)
            ax_d.figure.canvas.draw()
        ax.callbacks.connect("ylim_changed", convert_ax_h_to_d)
        ax_d.set_ylabel(u"Displacement [m/\u221AHz]")

    def plot_dict(noises):
        #use sorted to force a consistent ordering
        #The key lambda in sorted gets the (name, noise) pair and so nn[0] returns the name
        #the return tuple causes sorting by priority, with a fallback to lexical sort on the name
        for name, noise in sorted(noises.items(), key = lambda nn : (PRIORITY_MAP.get(nn[0], 100), nn[0])):
            if name in ['Freq', 'Total']:
                continue
            if isinstance(noise, dict):
                plot_dict(noise)
            else:
                logging.info("plotting '{}'...".format(name))
                noise = noises[name]
                stylekw = dict(
                    color = (0, 0, 0),
                    label = name,
                    lw = 3,
                )
                try:
                    style = STYLE_MAP[name]
                    stylekw.update(style)
                except KeyError:
                    pass
                ax.loglog(f, sqrt(noise), **stylekw)
    plot_dict(noises)

    ax.loglog(f, sqrt(noises['Total']), color='xkcd:black', alpha=0.6, label='Total', lw=4)

    ax.grid(
        True,
        which = 'both',
        lw    = 0.5,
        ls    = '-',
        alpha = 0.5,
    )

    ax.legend(
        ncol=2,
        fontsize = 'small',
    )

    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(u"Strain [1/\u221AHz]")
    ax.set_xlim([min(f), max(f)])
    ax.set_ylim([3e-25, 1e-21])
