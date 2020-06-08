from ..unit_def import unit
import numpy as np

def etminan(co2, ch4, n2o):
    """Calculate the radiative forcing from CO2, CH4 and N2O.

    This function uses the updated formulas of Etminan et al. (2016),
    (10.1002/2016GL071930) including the overlaps between CO2, methane and
    nitrous oxide.

    Args:
        co2: instance of CO2 object
        ch4: instance of CH4 object
        n2o: instance of N2O object
    """

    cbar = 0.5*(co2.reference_concentration + co2.concentrations)
    mbar = 0.5*(ch4.reference_concentration + ch4.concentrations)
    nbar = 0.5*(n2o.reference_concentration + n2o.concentrations)
    
    # units really come into their own here
    co2.effective_radiative_forcing = (
        (-2.4e-7/unit.ppm**2*(co2.concentrations - co2.reference_concentration)**2
         + 7.2e-4/unit.ppm*abs(co2.concentrations-co2.reference_concentration) 
         - 2.1e-4/unit.ppb * nbar + 5.36) 
        * np.log(co2.concentrations/co2.reference_concentration)
    ) * unit.watt / unit.m**2
    
    ch4.effective_radiative_forcing = (
        (-1.3e-6/unit.ppb*mbar - 8.2e-6/unit.ppb*nbar + 0.043) 
        * (np.sqrt(ch4.concentrations/unit.ppb) - np.sqrt(ch4.reference_concentration/unit.ppb))
    ) * unit.watt / unit.m**2
    
    n2o.effective_radiative_forcing = (
        (-8.0e-6/unit.ppm*cbar + 4.2e-6/unit.ppb*nbar - 4.9e-6/unit.ppb*mbar + 0.117)
        * (np.sqrt(n2o.concentrations/unit.ppb) - np.sqrt(n2o.reference_concentration/unit.ppb))
    ) * unit.watt / unit.m**2


def MN(M, N):
    return 0.47 * np.log(1 + 2.01e-5*(M*N)**(0.75) + 5.31e-15*M*(M*N)**(1.52))


def co2_log(C, Cpi, F2x=3.71):
    return F2x/np.log(2) * np.log(C/Cpi)


def myhre(C, Cpi, F2x=3.71):
    """Calculate the radiative forcing from CO2, CH4 and N2O.

    This uses the Myhre et al. (1998) relationships including the band
    overlaps between CH4 and N2O. It is also used in AR5.

    Reference: Myhre et al, 1998, JGR, doi: 10.1029/98GL01908

    Inputs:
        C: [CO2, CH4, N2O] concentrations, [ppm, ppb, ppb]
        Cpi: pre-industrial [CO2, CH4, N2O] concentrations

    Keywords:
        F2x: radiative forcing from a doubling of CO2.

    Returns:
        3-element array of radiative forcing: [F_CO2, F_CH4, F_N2O]
    """

    F = np.zeros(3)

    F[0] = co2_log(C[0], Cpi[0], F2x)
    F[1] = 0.036 * (np.sqrt(C[1]) - np.sqrt(Cpi[1])) - (
      MN(C[1],Cpi[2]) - MN(Cpi[1],Cpi[2]))
    F[2] = 0.12 * (np.sqrt(C[2]) - np.sqrt(Cpi[2])) - (
      MN(Cpi[1],C[2]) - MN(Cpi[1],Cpi[2])) 

    return F
