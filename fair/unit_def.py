import pint


class ScmUnitRegistry(pint.UnitRegistry):
    """
    Unit registry class for FaIR, adapted from OpenSCM.
    """
    
    def add_standards(self):
        """
        Add standard units.
        
        Has to be done separately because of pint's weird initializing.
        """

        self.define("a = 1 * year = annum = yr")
        self.define("h = hour")
        self.define("d = day")
        self.define("kt = 1000 * t")  # since kt is used for "knot" in the defaults
        self.define('volume_mixing_ratio = [vmr] = vmr')
        self.define('parts_per_million = 1e-6 * volume_mixing_ratio = ppm')
        self.define('parts_per_billion = 1e-9 * volume_mixing_ratio = ppb')
        self.define('parts_per_trillion = 1e-12 * volume_mixing_ratio = ppt')
        self.define('kilometric_ton = kt')
        

unit = ScmUnitRegistry()
unit.add_standards()
