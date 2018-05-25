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
        color = (0.41568627450980394, 0.23921568627450981, 0.60392156862745094),
    ),
    'Seismic': dict(
        color = (0.69411764705882351, 0.34901960784313724, 0.15686274509803921),
    ),
    'Newtonian Gravity': dict(
        color = (0.20000000000000001, 0.62745098039215685, 0.17254901960784313),
    ),
    'Suspension Thermal': dict(
        color = (0.12156862745098039, 0.47058823529411764, 0.70588235294117652),
    ),
    'Coating Brownian': dict(
        color = (0.8901960784313725, 0.10196078431372549, 0.10980392156862745),
    ),
    'Coating Thermo-Optic': dict(
        color = (0.98431372549019602, 0.60392156862745094, 0.59999999999999998),
        ls = '--',
    ),
    'ITM Thermo-Refractive': dict(
        color = (1.0, 0.49803921568627452, 0.0),
        ls = '--',
    ),
    'ITM Carrier Density': dict(
        color = (0.99215686274509807, 0.74901960784313726, 0.43529411764705883),
        ls = '--',
    ),
    'Substrate Brownian': dict(
        color = (0.69803921568627447, 0.87450980392156863, 0.54117647058823526),
        ls = '--',
    ),
    'Substrate Thermo-Elastic': dict(
        color = (0.65098039215686276, 0.80784313725490198, 0.8901960784313725),
        ls = '--',
    ),
    'Excess Gas': dict(
        color = (0.792156862745098, 0.69803921568627447, 0.83921568627450982),
        ls = '--',
    ),
}


def plot_noise(noises,):
    f = noises['Freq']

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
                plt.loglog(f, sqrt(noise), **stylekw)
    plot_dict(noises)

    plt.loglog(f, sqrt(noises['Total']), color='black', label='Total', lw=4)

    plt.grid(
        True,
        which='both',
        lw = .5,
        ls = '-',
        alpha = .5,
    )

    plt.legend(
        ncol=2,
        fontsize = 'small',
    )

    plt.xlabel('Frequency [Hz]')
    plt.ylabel(u"Strain [1/\u221AHz]")
    plt.xlim([min(f), max(f)])
    plt.ylim([3e-25, 1e-21])
