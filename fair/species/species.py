"""Species definition classes."""

import numpy as np

from ..unit_def import unit
from ..forcing.ghg import etminan
from ..constants import molwt
from ..constants.general import M_ATMOS
from ._forward import _forward


class Species:
    """Container class for all FaIR species."""
    pass
    

class GreenhouseGas(Species):
    """Greenhouse gas object (subclass of Species).

    Attributes:
        radiative_efficiency: Quantity of radiative efficiency of greenhouse
            gas (radiative forcing per unit concentration)
        reference_concentration: Quantity of the reference concentrations for
            radiative forcing calculations.
        concentration_t0: Quantity of the concentrations of greenhouse gas at
            timestep zero.
        molecular_weight: Quantity of molecular weight of the greenhouse gas.
        name: String describing the greenhouse gas name.
        efficacy: Optional float describing the temperature response per unit
            forcing or ``efficacy'' (Hansen et al., 2005:
            https://doi.org/10.1029/2005JD005776). Default 1.
    """

    def __init__(self, radiative_efficiency, reference_concentration, 
        concentration_t0, molecular_weight, lifetime, name, efficacy=1):
        """Initiator for GreenhouseGas class."""
        self.name = name
        self.efficacy = efficacy
        self.radiative_efficiency    = radiative_efficiency
        self.reference_concentration = reference_concentration
        self.concentration_t0        = concentration_t0
        self.molecular_weight        = molecular_weight
        self.lifetime                = lifetime
        self.kg_to_vmr               = (
            molwt.AIR/self.molecular_weight / M_ATMOS * 1.0 * unit.vmr)
        #if name in library:   # pandas datatable of stats or a csv?
            # load up all the stuff above
            # self.lifetime = library[name].lifetime
            # self.radiative_efficiency = library[name].radiative_efficiency and so on
    
    def set_emissions(self, anthropogenic, natural=0, timestep=1*unit.year, time_points=None):
        """Sets emissions of greenhouse gas.

        Args:
            anthropogenic: Quantity array of anthropogenic emissions.
            natural: optional Quantity array of natural emissions. Default 0.
            timestep: Time Quantity describing length of timestep.
            time_points: optional time Quantity array describing the time
                points corresponding to the emissions arrays.
        """
        self.anthropogenic_emissions = anthropogenic
        self.natural_emissions = natural
        self.emissions = anthropogenic + natural
        self.timestep = timestep
        if time_points is None:
            if hasattr(self.emissions, '__len__'):
                nt = len(self.emissions)
                time_points = (np.arange(nt) + 0.5) * timestep
# ERROR CHECKING: ensure sensible inputs
        self.time_points = time_points
    
    def set_concentrations(self, concentrations, timestep=1*unit.year, time_bounds=None):
        """Sets concentrations of greenhouse gases.

        Args:
            concentrations: Quantity array of concentrations.
            timestep: Time Quantity describing length of timestep.
            time_bounds: optional time Quantity array describing the time
                bounds corresponding to the concentration arrays.
        """

        self.concentrations = concentrations
        self.timestep = timestep
        if time_bounds is None:
            if hasattr(self.concentrations, '__len__'):
                nt = len(self.concentrations)
                time_bounds = (np.arange(nt) + 0.5) * timestep
# ERROR CHECKING: ensure sensible inputs
        self.time_bounds = time_bounds

    def set_effective_radiative_forcing(self, effective_radiative_forcing):
        """Sets effective radiative forcing of greenhouse gas.

        Args:
            effective_radiative_forcing: Quantity array of effective radiative
                forcing.
        """
        self.effective_radiative_forcing = effective_radiative_forcing
        
    def set_lifetime(self, lifetime):
        """Sets atmospheric e-folding lifetime of greenhouse gas.

        Args:
            lifetime: Quantity of lifetime
        """
        self.lifetime = lifetime

    def set_radiative_efficiency(self, radiative_efficiency):
        """Sets radiative efficiency of greenhouse gas.

        Args:
            lifetime: Quantity of radiative efficiency
        """
        self.radiative_efficiency = radiative_efficiency

    def calculate_concentrations(self):
        """Calculates greenhouse gas concentrations from emissions."""
        # for GHGs with no temperature feedback or inter-species dependence
        # only makes sense when emissions are defined
        if hasattr(self.emissions, '__len__'):
            nt = len(self.emissions) 
            # ensures that either emissions or natural_emissions are array - 
            # still need to check it is 1D
            emissions = self.emissions
        elif np.isscalar(self.emissions) and hasattr(self.time_points, '__len__'):
            nt = len(self.time_points)
            emissions = self.emissions * np.ones(nt)
        concentrations = np.ones(nt+1) * np.nan * unit.vmr
        time_bounds = self.time_points - self.timestep/2
        time_bounds = np.append(time_bounds, self.time_points[-1] + self.timestep/2)
        
        concentrations[0] = self.concentration_t0
        
        for t in range(nt):
            concentrations[t+1] = _forward(concentrations[t], emissions[t], 
                                         self.timestep, self.lifetime, self.kg_to_vmr)

        self.concentrations = concentrations
        self.time_bounds = time_bounds

    def calculate_forcing(self, scale=1):
        """Calculates effective radiative forcing from concentrations."""
        # need checks to determine that scale is scalar or array with same length as concentations
        self.effective_radiative_forcing = (
            (self.concentrations - self.reference_concentration) 
            * self.radiative_efficiency * scale
        )
        self.effective_radiative_forcing.ito(unit.watt / unit.m**2)


# move to a module
def _calculate_alpha(cumulative_emissions, airborne_emissions, temperature, r0,
        rC, rT, g0, g1, iirf_max = 97.0*unit.year):
    """Calculate scaling value for CO2 time constants.

    Args:
        cumulative_emissions: Quantity of total CO2 emissions this timestep
        airborne_emissions: Quantity of airborne CO2 emissions this timestep
        temperature: Quantity of global mean near-surface air temperature
        r0: Quantity of pre-industrial time-integrated airborne fraction
        rC: Quantity of sensitivity in time-integrated airborne fraction with
            total carbon in land and ocean sinks
        rT: Quantity of sensitivity in time-integrated airborne fraction with
            global mean temperature anomaly

    Returns:
        alpha: scaling factor for CO2 time constants.
    """

    iirf = r0 + rC * (cumulative_emissions-airborne_fraction) + rT * temperature
    iirf = np.abs(iirf)  # not sure I understand this - it should never be negative
    iirf = (iirf>iirf_max) * iirf_max + iirf * (iirf<iirf_max)
    alpha = g0 * np.sinh(iirf / g1)
    return alpha


def _step_concentration(carbon_boxes0, emissions, alpha, a, tau, timestep,
        concentration_t0, kg_to_vmr):
    """Calculate concentrations of CO2 from emissions

    Args:
        carbon_boxes0: Quantity array of atmospheric CO2 boxes at previous
            timestep.
        emissions: Quantity of emissions of CO2 this timestep.
        alpha: time constant scaling factor.
        a: carbon boxes partitioning fraction.
        tau: Quantity array of carbon boxes atmospheric time constants.
        timestep: Quantity of time step
        reference_concentration: reference concentration (pre-industrial)
            of CO2.
        kg_to_vmr: Quantity of unit conversion from kg CO2 to atmospheric
            concentration.

    Returns:
        concentration: Quantity of CO2 concentration at end of timestep.
        carbon_boxes1: Quantity array of atmospheric CO2 boxes at end of
            timestep.
        airborne_emissions: Quantity of atmospheric carbon store.
    """
    carbon_boxes1 = emissions * kg_to_vmr * a * alpha * (tau/timestep) * (1 - np.exp(-timestep/(alpha*tau))) + carbon_boxes0 * np.exp(-timestep/(alpha*tau))
    concentration = concentration_t0 + (carbon_boxes1 + carbon_boxes0) / 2
    airborne_emissions = np.sum(carbon_boxes0) / kg_to_vmr
    return concentration, carbon_boxes1, airborne_emissions


class CO2(GreenhouseGas):

    def __init__(self, radiative_efficiency, reference_concentration,
            concentration_t0, molecular_weight, lifetime, name, a, tau, r0, rC,
            rT, efficacy=1, iirf_horizon=100*unit.year, iirf_max=97*unit.year):
        super().__init__(radiative_efficiency, reference_concentration,
            concentration_t0, molecular_weight, lifetime, name, efficacy=1)
        self.a = a
        self.tau = tau
        self.r0 = r0
        self.rC = rC
        self.rT = rT
        self.iirf_horizon = iirf_horizon
        self.iirf_max = iirf_max

    def calculate_concentrations(self):
        """Calculates CO2 gas concentrations from emissions.

        Args:
            temperature: Quantity array of scenario temperature anomaly.
        """

        if hasattr(self.emissions, '__len__'):
            nt = len(self.emissions) 
            # ensures that either emissions or natural_emissions are array - 
            # still need to check it is 1D
            emissions = self.emissions
        elif np.isscalar(self.emissions) and hasattr(self.time_points, '__len__'):
            nt = len(self.time_points)
            emissions = self.emissions * np.ones(nt)

        # nothing here vectorised
        g1 = np.sum(self.a * self.tau * (1 - (1 + self.iirf_horizon/self.tau)
            * np.exp(-self.iirf_horizon/self.tau)), axis=-1)
        g0 = 1/(np.sinh(np.sum(self.a * self.tau * (1 - 
            np.exp(-self.iirf_horizon/self.tau)), axis=-1)/g1))
        alpha1 = _calculate_alpha(cumulative_emissions)

        concentrations = np.ones(nt+1) * np.nan * unit.vmr
        #time_bounds = self.time_points - self.timestep/2
        #time_bounds = np.append(time_bounds, self.time_points[-1] + self.timestep/2)
        airborne_emissions = np.zeros(nt+1) * unit.Gt
        alpha = np.ones(nt) * np.nan
        concentrations[0] = self.concentration_t0
        cumulative_emissions = np.cumsum(emissions)
        carbon_boxes = np.zeros_like(self.a)

        for t in range(nt):
            alpha[t] = _calculate_alpha(cumulative_emissions[t], 
                airborne_emissions[t], self.temperature[t], self.r0, self.rC, 
                self.rT, g0, g1, iirf_max=self.iirf_max)        
            concentrations[t+1], carbon_boxes, airborne_emissions[t+1] = (
                _step_concentration(carbon_boxes, emissions[t], alpha[t], 
                    self.a, self.tau, self.timestep,
                    self.concentration_t0, self.kg_to_vmr)
            )

        self.alpha = alpha
        self.concentrations = concentrations
        self.time_bounds = time_bounds
        self.cumulative_emissions = cumulative_emissions
        self.airborne_emissions = airborne_emissions
        self.airborne_fraction = airborne_emissions/cumulative_emissions

    def calculate_forcing(self, ch4, n2o):
        etminan(self, ch4, n2o)


class CH4(GreenhouseGas):
    def calculate_forcing(self, co2, n2o):
        etminan(co2, self, n2o)


class N2O(GreenhouseGas):
    def calculate_forcing(self, co2, ch4):
        etminan(co2, ch4, self)
