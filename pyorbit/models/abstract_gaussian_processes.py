from pyorbit.subroutines.common import np
from pyorbit.keywords_definitions import *

class AbstractGaussianProcesses(object):

    def __init__(self, *args, **kwargs):
        self.use_stellar_rotation_period = False
        self.use_stellar_activity_decay = False

        """ Additional flags for those GP kernels that are not inherently connected to the stellar activity
            but can be used to model it
        """
        self.use_activity_Prot = False
        self.use_activity_Pdec = False
        self.use_shared_scale = False
        self.use_shared_decay = False
        self.use_shared_hyperparameters = False

        self.model_class = 'gaussian_process'


    def _prepare_hyperparameter_conditions(self, mc, **kwargs):


        if kwargs.get('hyperparameters_condition', False):
            print('GP model: {0:20s} hyperparameters_condition:    True'.format(self.model_name))
            self.hyper_condition = self._hypercond_01
        else:
            self.hyper_condition = self._hypercond_00

        if kwargs.get('rotation_decay_condition', False):
            print('GP model: {0:20s} rotation_decay_condition:     True'.format(self.model_name))
            self.rotdec_condition = self._hypercond_02
        else:
            self.rotdec_condition = self._hypercond_00

        if kwargs.get('halfrotation_decay_condition', False):
            print('GP model: {0:20s} halfrotation_decay_condition: True'.format(self.model_name))
            self.halfrotdec_condition = self._hypercond_03
        else:
            self.halfrotdec_condition = self._hypercond_00

        self.rotdec_factor_condition = self._hypercond_00
        if kwargs.get('decay_rotation_factor', False):
            print('GP model: {0:20s} decay_rotation_factor:        True'.format(self.model_name))
            self.decay_rotation_factor = kwargs.get('decay_rotation_factor', 0)
            self.rotdec_factor_condition = self._hypercond_04

        if kwargs.get('rotation_decay_factor', False):
            print('GP model: {0:20s} rotation_decay_factor:        True'.format(self.model_name))
            self.decay_rotation_factor = kwargs.get('rotation_decay_factor', 0)
            self.rotdec_factor_condition = self._hypercond_04

    def _prepare_shared_hyperparameters(self, pam_scale=None, pam_decay=None, **kwargs):

        for keyword in keywords_shared_hyperparameters:
            self.use_shared_hyperparameters =  kwargs.get(keyword, self.use_shared_hyperparameters)
        if self.use_shared_hyperparameters:
            self.use_shared_scale  = True
            self.use_shared_decay  = True
            pams_copy = self.list_pams_dataset.copy()
            for pam in pams_copy:
                self.list_pams_common.update([pam])
                self.list_pams_dataset.discard(pam)

        for keyword in keywords_shared_timescale:
            self.use_shared_scale =  kwargs.get(keyword, self.use_shared_scale)
        if self.use_shared_scale and pam_scale is not None:
            self.list_pams_common.update([pam_scale])
            self.list_pams_dataset.discard(pam_scale)

        for keyword in keywords_shared_decay:
            self.use_shared_decay =  kwargs.get(keyword, self.use_shared_decay)
        if self.use_shared_decay and pam_decay is not None:
            self.list_pams_common.update([pam_decay])
            self.list_pams_dataset.discard(pam_decay)

    def _check_extra_conditions(self, **kwargs):

        flag_check = [self.use_activity_Prot, self.use_activity_Pdec, self.use_stellar_rotation_period, self.use_stellar_activity_decay]
        print(flag_check )
        if sum(flag_check) > 1:
            raise ValueError('The rotation and decay flags are mutually exclusive. Please choose one of the following options: '
                                'use_activity_Prot, use_activity_Pdec, use_stellar_rotation_period, use_stellar_activity_decay')

    def _prepare_rotation_replacement(self, mc, parameter_name ='Prot', common_pam=True, check_common=True, **kwargs):

        if check_common:
            for common_ref in self.common_ref:
                if mc.common_models[common_ref].model_class == 'activity':
                    for keyword in keywords_stellar_rotation:

                        self.use_stellar_rotation_period = getattr(mc.common_models[common_ref], keyword, self.use_stellar_rotation_period)
                    break

        for keyword in keywords_activity_Prot:
            self.use_activity_Prot = kwargs.get(keyword, self.use_activity_Prot)

        if self.use_activity_Prot:
            self.list_pams_common.update(['Prot'])
            if common_pam:
                self.list_pams_common.discard(parameter_name)
            else:
                self.list_pams_dataset.discard(parameter_name)

        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        if self.use_stellar_rotation_period:
            self.list_pams_common.update(['rotation_period'])
            if common_pam:
                self.list_pams_common.discard(parameter_name)
            else:
                self.list_pams_dataset.discard(parameter_name)

        flag_check = [self.use_activity_Prot, self.use_stellar_rotation_period]
        if sum(flag_check) > 1:
            raise ValueError('The rotation and Prot flags are mutually exclusive. Please choose one of the following options: '
                                'use_activity_Prot, use_stellar_rotation_period')

    def _prepare_decay_replacement(self, mc, parameter_name ='Pdec', common_pam=True,  check_common=True, **kwargs):

        if check_common:
            for common_ref in self.common_ref:
                if mc.common_models[common_ref].model_class == 'activity':
                    for keyword in keywords_stellar_rotation:
                        self.use_stellar_activity_decay = getattr(mc.common_models[common_ref], keyword, self.use_stellar_activity_decay)
                    break

        for keyword in keywords_activity_Pdec:
            self.use_activity_Pdec = kwargs.get(keyword, self.use_activity_Pdec)

        if self.use_activity_Pdec:
            self.list_pams_common.update(['Pdec'])
            if common_pam:
                self.list_pams_common.discard(parameter_name)
            else:
                self.list_pams_dataset.discard(parameter_name)


        for keyword in keywords_stellar_activity_decay:
            self.use_stellar_activity_decay = kwargs.get(keyword, self.use_stellar_activity_decay)

        if self.use_stellar_activity_decay:
            self.list_pams_common.update(['activity_decay'])
            if common_pam:
                self.list_pams_common.discard(parameter_name)
            else:
                self.list_pams_dataset.discard(parameter_name)

        flag_check = [self.use_activity_Pdec, self.use_stellar_activity_decay]
        if sum(flag_check) > 1:
            raise ValueError('The rotation and Prot flags are mutually exclusive. Please choose one of the following options: '
                                'use_activity_Pdec, use_stellar_activity_decay')


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

        if self.use_activity_Prot:
            parameter_values[replace_rotation] = parameter_values['Prot']

        if self.use_activity_Pdec:
            parameter_values[replace_decay] = parameter_values['Pdec']

        if self.use_stellar_rotation_period:
            parameter_values[replace_rotation] = parameter_values['rotation_period']

        if self.use_stellar_activity_decay:
            parameter_values[replace_decay] = parameter_values['activity_decay']

    def check_hyperparameter_values(self, parameter_values, pam_scale=None, pam_decay=None):

        if pam_scale is not None:
            parameter_values['Prot'] = parameter_values[pam_scale]
        if pam_decay is not None:
            parameter_values['Pdec'] = parameter_values[pam_decay]

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

