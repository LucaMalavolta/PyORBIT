from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerEmcee(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        """ pyde/emcee variabless """
        self.pyde_dir_output = None
        self.emcee_dir_output = None

        self.emcee_parameters = {'nsave': 0,
                                 'npop_mult': 4,
                                 'thin': 1,
                                 'nburn': 0,
                                 'multirun': None,
                                 'multirun_iter': 20,
                                 'shutdown_jitter': False,
                                 'use_threading_pool': True,
                                 'starts_relative_dispersion': True,
                                 #'progress_bar': True
                                }

        self.pyde_parameters = {'ngen': 8000,
                                'npop_mult': 4,
                                'shutdown_jitter': False,
                                'use_threading_pool': True,
                                }
