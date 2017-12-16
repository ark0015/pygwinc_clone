from numpy import sqrt
import matplotlib.pyplot as plt

import logging

COLOR_MAP = {
    'Quantum Vacuum': (0.41568627450980394, 0.23921568627450981, 0.60392156862745094),
    'Seismic': (0.69411764705882351, 0.34901960784313724, 0.15686274509803921),
    'Newtonian Gravity': (0.20000000000000001, 0.62745098039215685, 0.17254901960784313),
    'Suspension Thermal': (0.12156862745098039, 0.47058823529411764, 0.70588235294117652),
    'Coating Brownian': (0.8901960784313725, 0.10196078431372549, 0.10980392156862745),
    'Coating Thermo-Optic': (0.98431372549019602, 0.60392156862745094, 0.59999999999999998),
    'ITM Thermo-Refractive': (1.0, 0.49803921568627452, 0.0),
    'ITM Carrier Density': (0.99215686274509807, 0.74901960784313726, 0.43529411764705883),
    'Substrate Brownian': (0.69803921568627447, 0.87450980392156863, 0.54117647058823526),
    'Substrate Thermo-Elastic': (0.65098039215686276, 0.80784313725490198, 0.8901960784313725),
    'Excess Gas': (0.792156862745098, 0.69803921568627447, 0.83921568627450982),
    }


def plot_noise(noises):
    f = noises['Freq']

    def plot_dict(noises):
        for name, noise in noises.items():
            if name in ['Freq', 'Total']:
                continue
            if isinstance(noise, dict):
                plot_dict(noise)
            else:
                logging.info("plotting '{}'...".format(name))
                noise = noises[name]
                try:
                    color = COLOR_MAP[name]
                except KeyError:
                    color = (0, 0, 0)
                plt.loglog(f, sqrt(noise), color=color, label=name, lw=3)
    plot_dict(noises)

    plt.loglog(f, sqrt(noises['Total']), color='black', label='Total', lw=4)

    plt.grid(True, which='both')
    plt.legend(ncol=2)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel(u"Strain [1/\u221AHz]")
    plt.xlim([min(f), max(f)])
    plt.ylim([3e-25, 1e-21])
