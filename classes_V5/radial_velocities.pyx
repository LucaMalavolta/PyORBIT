from common import *
from abstract_model import *

"""
New changes: 
    mc.variables  is now called  mc.converter
    mv.var_list  is now called  mc.variable_index
"""


class RV_keplerian(AbstractModel):

    def __init__(self, model_name, common_ref, mc):
        self.model_class = 'rv_keplerian'
        self.model_name = model_name
        self.common_ref = common_ref
        self.common_model = mc.common_models[self.common_ref]

        self.list_pams_common = {
            'P': 'LU',  # Period, log-uniform prior
            'K': 'LU',  # RV semi-amplitude, log-uniform prior
            'f': 'U',  # RV vurve phase, log-uniform prior
            'e': 'U',  # eccentricity, uniform prior - to be fixed
            'o': 'U'}  # argument of pericenter
        self.list_pams_dataset = {}

        """ This is the list of planets that should use that model
        This list must be filled somehow
        """
        self.planet_list = {}

    def define_common_special_bounds(self, mc, var):
        if not(var == "e" or var == "o"):
            return False

        if 'e' in self.common_model.fix_list or \
           'o' in self.common_model.fix_list:
            return False

        if 'coso' in mc.variable_sampler[self.common_ref] or \
            'esino' in mc.variable_sampler[self.common_ref]:
            return False

        self.common_model.transformation['e'] = get_2var_e
        self.common_model.variable_index['e'] = [mc.ndim, mc.ndim + 1]
        self.common_model.transformation['o'] = get_2var_o
        self.common_model.variable_index['o'] = [mc.ndim, mc.ndim + 1]

        self.common_model.variable_sampler['ecoso'] = mc.ndim
        self.common_model.variable_sampler['esino'] = mc.ndim + 1
        mc.bounds_list.append(self.common_model.default_bounds['ecoso'])
        mc.bounds_list.append(self.common_model.default_bounds['esino'])
        mc.ndim += 2

        return True

    def define_common_special_starting_point(self, mc, var):
        '''eccentricity and argument of pericenter require a special treatment
         since they can be provided as fixed individual values or may need to be combined
         in ecosw and esinw if are both free variables'''

        if not(var == "e" or var == "o"):
            return False

        if 'ecoso' in self.common_model.variable_sampler and \
           'esino' in self.common_model.variable_sampler:

            if 'e' in self.common_model.starts and 'o' in self.common_model.starts:
                mc.starting_point[self.common_model.variable_sampler['ecoso']] = \
                    np.sqrt(self.common_model.starts['e']) * np.cos(self.common_model.starts['o'])
                mc.starting_point[self.common_model.variable_sampler['esino']] = \
                    np.sqrt(self.common_model.starts['e']) * np.sin(self.common_model.starts['o'])

            elif 'ecoso' in self.common_model.starts and 'esino' in self.common_model.starts:
                mc.starting_point[self.common_model.variable_sampler['ecoso']] = self.common_model.starts['ecoso']
                mc.starting_point[self.common_model.variable_sampler['ecoso']] = self.common_model.starts['esino']

        return True
