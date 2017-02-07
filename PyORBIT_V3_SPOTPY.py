from PyORBIT_V3_Classes import *
import numpy as np
import spotpy
import cPickle as pickle
import os
import argparse

parser = argparse.ArgumentParser(prog='PyORBIT_V3_SPOTPY.py', description='PyDE+emcee runner')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('config_file', type=str, nargs=1, help='config file')


args = parser.parse_args()


class spotpy_setup(object):
    def __init__(self, args):
        file_conf = args.config_file[0]
        self.mc = ModelContainerSPOTPY()
        yaml_parser(file_conf, self.mc)
        dir_output = './' + self.mc.planet_name + '/'

        self.mc.create_bounds()
        print 'Dimensions = ', self.mc.ndim
        print '   '
        print 'Variable list:', self.mc.variable_list
        print
        print 'Variable bounds:', self.mc.bounds
        print

        population = pickle.load(open(dir_output + 'pyde_pops.pick', 'rb'))
        pyde_mean = pickle.load(open(dir_output + 'pyde_mean.pick', 'rb'))
        self.mc.recenter_bounds(pyde_mean, population)

        self.mc.results_resumen(pyde_mean)
        self.params = self.mc.spotpy_priors(pyde_mean)
        print pyde_mean
    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        theta=np.array(vector)
        return self.mc(theta)

    def evaluation(self):
        return self.mc.output_concatenated()

    def objectivefunction(self, simulation, evaluation):
        objectivefunction= -spotpy.objectivefunctions.rmse(simulation, evaluation)
        return objectivefunction

results = []
spotpy_setup = spotpy_setup(args)
rep = 1000

sampler = spotpy.algorithms.mc(spotpy_setup,    dbname='aaaMC',    dbformat='csv')
sampler.sample(rep)
results.append(sampler.getdata())


sampler=spotpy.algorithms.demcz(spotpy_setup, dbname='aaaDEMCz', dbformat='csv')
sampler.sample(rep, nChains=30)
results.append(sampler.getdata())


#spotpy.analyser.plot_parametertrace(results)

algorithms=['MC', 'DEMCz']
#spotpy.analyser.plot_heatmap_griewank(results,algorithms)

evaluation = spotpy_setup.evaluation()
spotpy.analyser.plot_likelihoodtraces(results, evaluation, algorithms)