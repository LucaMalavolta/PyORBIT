from common import *
from model_container_abstract import ModelContainer


class ModelContainerEmcee(ModelContainer):

    """ pyde/emcee variabless """
    pyde_dir_output = None
    emcee_dir_output = None

    emcee_parameters = {'nsave': 0,
                        'npop_mult': 4,
                        'thin': 1,
                        'nburn':0,
                        'multirun': None,
                        'multirun_iter': 20,
                        'version': '2.2.1',
                        'shutdown_jitter': False
                        }

    pyde_parameters = {'ngen': 8000,
                       'npop_mult': 4,
                       'shutdown_jitter': False
                       }

