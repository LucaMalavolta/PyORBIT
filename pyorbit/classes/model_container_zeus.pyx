from pyorbit.classes.model_container_abstract import ModelContainer


class ModelContainerZeus(ModelContainer):

    def __init__(self):
        super(self.__class__, self).__init__()

        """ pyde/zeus variabless """
        self.pyde_dir_output = None
        self.zeus_dir_output = None

        self.zeus_parameters = {'nsave': 0,
                                 'npop_mult': 4,
                                 'thin': 1,
                                 'nburn': 0,
                                 'multirun': None,
                                 'multirun_iter': 20,
                                 'shutdown_jitter': False,
                                 'use_threading_pool': True
                                }

        self.pyde_parameters = {'ngen': 8000,
                                'npop_mult': 4,
                                'shutdown_jitter': False,
                                'use_threading_pool': True
                                }
