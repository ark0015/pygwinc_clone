from numpy import sqrt

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
        # color = 'xkcd:vibrant purple',
        color = '#ad03de',
    ),
    'Seismic': dict(
        # color = 'xkcd:brown',
        color = '#653700',
    ),
    'Newtonian Gravity': dict(
        # color = 'xkcd:green',
        color = '#15b01a',
    ),
    'Atmospheric Infrasound': dict(
        # color = 'xkcd:neon green',
        color = '#0cff0c',
        ls = '--',
    ),
    'Suspension Thermal': dict(
        # color = 'xkcd:deep sky blue',
        color = '#0d75f8',
    ),
    'Coating Brownian': dict(
        # color = 'xkcd:fire engine red',
        color = '#fe0002',
    ),
    'Coating Thermo-Optic': dict(
        # color = 'xkcd:bright sky blue',
        color = '#02ccfe',
        ls = '--',
    ),
    'ITM Thermo-Refractive': dict(
        # color = 'xkcd:dark sky blue',
        color = '#448ee4',
        ls = '--',
    ),
    'ITM Carrier Density': dict(
        # color = 'xkcd:grey',
        color = '#929591',
        ls = '--',
    ),
    'Substrate Brownian': dict(
        # color = 'xkcd:pumpkin orange',
        color = '#fb7d07',
        ls = '--',
    ),
    'Substrate Thermo-Elastic': dict(
        # color = 'xkcd:golden',
        color = '#f5bf03',
        ls = '--',
    ),
    'Excess Gas': dict(
        # color = 'xkcd:baby shit brown',
        color = '#ad900d',
        ls = '--',
    ),
}


def plot_noise(
        ifo,
        noises,
        ax=None,
        displacement=True,
):
    """Plot a GWINC noise budget from calculated noises

    If an axis handle is provided it will be used for the plot.
    `displacement` may be set to False to supress the right had
    displacement axis.

    Returns the figure handle used.

    """
    f = noises['Freq']

    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
    else:
        fig = ax.figure

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

    ax.loglog(f, sqrt(noises['Total']), color='#000000', alpha=0.6, label='Total', lw=4)

    ax.grid(
        True,
        which='both',
        lw=0.5,
        ls='-',
        alpha=0.5,
    )

    ax.legend(
        ncol=2,
        fontsize='small',
    )

    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(u"Strain [1/\u221AHz]")
    ax.set_xlim([min(f), max(f)])
    ax.set_ylim([3e-25, 1e-21])

    if displacement:
        ax_d = ax.twinx()
        ax_d.set_yscale('log')
        y1, y2 = ax.get_ylim()
        ax_d.set_ylim(y1*ifo.Infrastructure.Length, y2*ifo.Infrastructure.Length)
        ax_d.set_ylabel(u"Displacement [m/\u221AHz]")

    return fig
