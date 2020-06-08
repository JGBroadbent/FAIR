import numpy as np

def _forward(c0, em, ts, lt, vm):
    """Calculates concentrations.

    Args:
        c0: Quantity of concentrations at time_bound t
        em: Quantity of emissions at time step t
        ts: Quantity of timestep length
        lt: Quantity of greenhouse gas atmospheric lifetime
        vm: Quantity of kilogram to volume mixing ratio

    Returns:
        c1: Quantity of concentrations at time bound t+1
    """
    c1 = c0*np.exp(-ts/lt) + lt*em*vm * (1.-np.exp(-ts/lt))
    return c1

        
