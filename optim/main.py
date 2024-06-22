#!/usr/bin/python

import numpy as np
import pandas as pd
import pygmo as pg
from problem import SignalProcessingParams
import matplotlib.pyplot as plt

""" -------------------------- Options -------------------------- """
# Pandas display options
pd.options.display.max_columns = 139  #(26 + 32 + 29 + 25 + 27)
""" ----------------------------------------------------------------------- """

""" -------------------------- Misc Functions -------------------------- """
def convergence_plot(algos, figsize=(16, 50)) -> None:
  """
  Convergence plot from all islands

  :param algos:
  list of algorithms

  :return:
  """
  # Plot style
  plt.style.use('ggplot')
  # Extracting algorithms and logs
  # algos = [isl.get_algorithm() for isl in islands]
  udas = [alg.extract(algorithm_name) for alg in algos]
  logs = [uda.get_log() for uda in udas]
  # Creating figure
  fig = plt.figure(figsize=figsize)
  # Number of rows and columns
  cols = 2
  rows = (len(algos)//2) if not (len(algos) % 2) else len(algos)//2 + 1
  for num in range(1, len(algos)+1):
    ax = fig.add_subplot(rows, cols, num)
    ax.plot([entry[0] for entry in logs[num-1]],
            [entry[1] for entry in logs[num-1]],
            'k-',
            linewidth=2.0,
            label='Minimal capital cost')
    ax.plot([entry[0] for entry in logs[num-1]],
            [entry[2] for entry in logs[num-1]],
            'r--',
            linewidth=2.0,
            label='Average capital cost')
    ax.set_title('cm: {}, pm: {}'.format(udas[num-1].cm, udas[num-1].pm))
    ax.set_xlabel('Generation number')
    ax.set_ylabel('Capital cost (thousand USD)')
    ax.legend()
  plt.show()
  return


def get_best_individual_cm_pm(pops, algos):
  """
  Fitness najbolje jedinke svake populacije s pripadajućim crossover
  probability (cm) i mutation probability (pm)

  :param pops:
  List of populations

  :param algos:
  List of algorithms

  :return:
  """
  pops_best_f = [p.champion_f[0] for p in pops]
  udas = [alg.extract(algorithm_name) for alg in algos]
  probabs = [(uda.cm, uda.pm) for uda in udas]
  for i, best_f, probab in zip(range(len(probabs)), pops_best_f, probabs):
    print('Run {}'.format(i+1))
    print('\tMin Capital cost: {} × 1000€, cm: {}, pm: {}'.format(
        best_f, probab[0], probab[1]
    ))
  return


""" ------------------------------------------------------------------------ """

if __name__ == '__main__':
  """ -------------------------- Input Data -------------------------- """

  emg_grip_data = pd.read_csv(
    'data/2022-07-14-12-53-38_grip_emg_test_data_upsampled_aligned.csv',
    sep=',',
    decimal='.',
    dtype=np.float64,
    header=0,
    index_col=False
  )

  """ ------------------------ Filtering parameters ------------------------ """
  window_size_relative = 0.5

  """ ---------------------------  Algorithm Data -------------------------- """
  # Population size
  pop_size = 60
  # 20
  # Random number generation seed
  rnd_seed = None
  # Crossover probability
  cm = 0.3
  # Mutation probability
  pm = 0.5
  # Maximum generations
  max_gen = 2500
  # 1000
  # Phi termination criterion threshold, 1 >> eps > 0
  eps = 0.001
  # Logging level
  log_lvl = 1
  # Mutation probability list
  pm_lst = np.linspace(start=0.1, stop=0.7, num=4)  # 4 0.1 ÷ 0.7
  # Crossover probability list
  cm_lst = np.linspace(start=0.1, stop=0.7, num=5)  # 5 0.1 ÷ 0.7
  # Number of islands
  isl_num = pm_lst.shape[0] * cm_lst.shape[0]

  """ ------------------- Problem, Population, Algorithm ------------------- """
  # # Problem instance and variables
  prob_var = SignalProcessingParams(
    time_dat=emg_grip_data['measure_time'].to_numpy(),
    emg_dat=emg_grip_data['emg_ch_1'].to_numpy(),
    grip_dat=emg_grip_data['grip_force_inter'].to_numpy(),
    window_size_relative=window_size_relative
  )
  prob = pg.problem(prob_var)
  print(prob)
  # # Population instance with initial values
  # pop = pg.population(prob)
  # [pop.push_back(next(GenerateIndividual)) for _ in range(pop_size)]
  # # Algorithm instance
  # algo = pg.algorithm(GA_based_optimize(
  #   pm, cm, prob_var.left_ind, prob_var.mid_ind, prob_var.right_ind,
  #   prob_var.Nmax, gen=max_gen, eps=eps, seed=rnd_seed, log_lvl=log_lvl
  # ))
  # # Generating isl_num initial populations
  # pop_lst = [pg.population(prob) for _ in range(isl_num)]
  # [[p.push_back(next(GenerateIndividual)) for p in pop_lst]
  #   for _ in range(pop_size)]
  # # Generating pm_lst × cm_lst algorithms with different crossover and
  # # mutation probability
  # alg_lst = [pg.algorithm(GA_based_optimize(
  #   p, c, prob_var.left_ind, prob_var.mid_ind, prob_var.right_ind,
  #   prob_var.Nmax, gen=max_gen, eps=eps, seed=rnd_seed, log_lvl=log_lvl
  # )) for c in cm_lst for p in pm_lst]
  # # Creating islands
  # islands = [
  #   pg.island(algo=a, pop=p, udi=pg.mp_island()) for a, p in zip(alg_lst, pop_lst)
  #   ]