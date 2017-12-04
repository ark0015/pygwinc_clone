from numpy import pi, sqrt
import scipy.constants

def computeFCParams(ifo, fcParams):
    """Compute ideal filter cavity Tin, detuning [Hz] and bandwidth [Hz]

    """
    # FC parameters
    c = scipy.constants.c
    fsrFC = c / (2 * fcParams.L)
    lossFC = fcParams.Lrt + fcParams.Te
  
    # detuning and cavity bandwidth (D&D paper P1400018 and/or PRD)
    eps = 4 / (2 + sqrt(2 + 2 * sqrt(1 + (4 * pi * ifo.gwinc.fSQL / (fsrFC * lossFC))**4)))
    s1eps = sqrt(1 - eps)

    # cavity bandwidth [Hz]
    gammaFC = ifo.gwinc.fSQL / sqrt(s1eps + s1eps**3)
    # cavity detuning [Hz]
    detuneFC = s1eps * gammaFC
  
    # input mirror transmission
    TinFC = 4 * pi * gammaFC / fsrFC - lossFC
    if TinFC < lossFC:
        raise RuntimeError('IFC: Losses are too high! %.1f ppm max.' % 1e6 * gammaFC / fsrFC)
  
    # Add to fcParams structure
    fcParams.Ti = TinFC
    fcParams.fdetune = -detuneFC
    fcParams.gammaFC = gammaFC
    fcParams.fsrFC = fsrFC

    return fcParams
