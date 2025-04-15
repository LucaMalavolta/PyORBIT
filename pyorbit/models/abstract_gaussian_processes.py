from pyorbit.subroutines.common import np
from pyorbit.keywords_definitions import *

class AbstractGaussianProcesses(object):

    def __init__(self, *args, **kwargs):
        self.use_stellar_rotation_period = False
        self.use_stellar_activity_decay = False
        self.model_class = 'gaussian_process'


    def _prepare_hyperparameter_conditions(self, mc, **kwargs):

        if kwargs.get('hyperparameters_condition', False):
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

        if kwargs.get('halfrotation_decay_condition', False):
            self.halfrotdec_condition = self._hypercond_03
        else:
            self.halfrotdec_condition = self._hypercond_00

        self.rotdec_factor_condition = self._hypercond_00   
        if kwargs.get('decay_rotation_factor', False):
            self.decay_rotation_factor = kwargs.get('decay_rotation_factor', 0)
            self.rotdec_factor_condition = self._hypercond_04

        if kwargs.get('rotation_decay_factor', False):
            self.decay_rotation_factor = kwargs.get('rotation_decay_factor', 0)
            self.rotdec_factor_condition = self._hypercond_04

    def _prepare_rotation_replacement(self, mc, parameter_name ='Prot', common_pam=True, **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                for keyword in keywords_stellar_rotation:

                    self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], keyword, self.use_stellar_rotation_period)
                break


        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)


        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            if common_pam:
                self.list_pams_common.discard(parameter_name)
            else:
                self.list_pams_dataset.discard(parameter_name)

    def _prepare_decay_replacement(self, mc, parameter_name ='Pdec' , **kwargs):

        for common_ref in self.common_ref:
            if mc.common_models[common_ref].model_class == 'activity':
                self.use_stellar_activity_decay = getattr(mc.common_models[common_ref], 'use_stellar_activity_decay', False)
                break

        for keyword in keywords_stellar_activity_decay:
            self.use_stellar_activity_decay = kwargs.get(keyword, self.use_stellar_activity_decay)

        if self.use_stellar_activity_decay:
            self.list_pams_common.update(['activity_decay'])
            self.list_pams_common.discard(parameter_name)

    def _set_derivative_option(self, mc, dataset, return_flag=False, **kwargs):

        if 'derivative'in kwargs:
            use_derivative = kwargs['derivative'].get(dataset.name_ref, False)
        elif dataset.name_ref in kwargs:
            use_derivative = kwargs[dataset.name_ref].get('derivative', False)
        else:
            if dataset.kind == 'H-alpha' or \
                dataset.kind == 'S_index' or \
                dataset.kind == 'Ca_HK' or \
                dataset.kind == 'FWHM':
                    use_derivative = False
            else:
                use_derivative = True

        """ instead of taking an action on the parameter, the flag is returned"""
        if return_flag:
            return use_derivative

        if not use_derivative:
            self.fix_list[dataset.name_ref] = {'rot_amp': [0., 0.]}

    def update_parameter_values(self,
                                parameter_values,
                                prepend='',
                                replace_rotation='Prot',
                                replace_decay='Pdec'):

        if self.use_stellar_rotation_period:
            parameter_values[replace_rotation] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values[replace_decay] = parameter_values['activity_decay']

    def check_hyperparameter_values(self, parameter_values):

        if not self.hyper_condition(parameter_values):
            return -np.inf
        if not self.rotdec_condition(parameter_values):
            return -np.inf
        if not self.halfrotdec_condition(parameter_values):
            return -np.inf
        if not self.rotdec_factor_condition(parameter_values):
            return -np.inf

        return True

    def _hypercond_04(self, parameter_values):
        #Condition on Rotation period and decay timescale
        return parameter_values['Pdec'] > self.decay_rotation_factor * parameter_values['Prot']

    @staticmethod
    def _hypercond_00(parameter_values):
        #Condition from Rajpaul 2017, Rajpaul+2021
        return True

    @staticmethod
    def _hypercond_01(parameter_values):
        # Condition from Rajpaul 2017, Rajpaul+2021
        # Taking into account that Pdec^2 = 2*lambda_2^2
        return parameter_values['Pdec']**2 > (3. / 2. / np.pi) * parameter_values['Oamp']**2 * parameter_values['Prot']**2

    @staticmethod
    def _hypercond_02(parameter_values):
        #Condition on Rotation period and decay timescale
        return parameter_values['Pdec'] > 2. * parameter_values['Prot']

    @staticmethod
    def _hypercond_03(parameter_values):
        #Condition on Rotation period and decay timescale
        return parameter_values['Pdec'] > 0.5 * parameter_values['Prot']

