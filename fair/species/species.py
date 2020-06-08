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
        concentration_t0, molecular_weight, name, efficacy=1):
        """Initiator for GreenhouseGas class."""
        self.name = name
        self.efficacy = efficacy
        self.radiative_efficiency    = radiative_efficiency
        self.reference_concentration = reference_concentration
        self.concentration_t0        = concentration_t0
        self.molecular_weight        = molecular_weight
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


class CO2(GreenhouseGas):
    def calculate_forcing(self, ch4, n2o):
        etminan(self, ch4, n2o)


class CH4(GreenhouseGas):
    def calculate_forcing(self, co2, n2o):
        etminan(co2, self, n2o)


class N2O(GreenhouseGas):
    def calculate_forcing(self, co2, ch4):
        etminan(co2, ch4, self)
