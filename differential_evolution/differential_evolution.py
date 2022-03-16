import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import os
import time

import utils
import cost

# Directory of differential evolution output
dir_output = "differential_evolution/output/"

array_param = "epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c" 
bounds = [(0.1, 0.99), (0.1, 0.99), (0.1, 0.99), (20, 60), (0.1, 10), (0.01, 2), (0.1, 2), (1, 15), (1, 2), (1, 2), (10, 25)] 

# Population size
pop_size = 30

# Mutation factor [0,2]
mutate = 0.7

# Recombination rate [0,1]
recombination = 0.5

# Max number of generations (iterations)
max_iter = 30

if __name__ == "__main__":

    # Creation of necessary directory
    if not os.path.isdir(dir_output):
        os.system("mkdir" + dir_output)

    output = open(dir_output + "report_DE.txt", "a")
    output.writelines("\n\n =========== Experiment ===========\n")
    output.writelines("\nPopulation size: " + str(pop_size))
    output.writelines("\nNumber of generations: " + str(max_iter) + "\n")
    output.writelines("\nParameters: " + array_param + "\n")
    output.close()

    pat_cont = 1
    
    patients = utils.reads_patients_names()

    for p in patients:
        output = open(dir_output + "report_DE.txt", "a")
        output.writelines("- " + p)

        # Gets patient experimental data
        exp_time, exp_viral_load = utils.reads_experimental_data(p)
        
        # Diferential evolution with time 
        t_ini = time.perf_counter()
        solution = differential_evolution(cost.model_cost, bounds, args=(exp_time, exp_viral_load, True), maxiter=max_iter, popsize=pop_size, mutation=mutate, recombination=recombination, tol=0.1)
        t_fin = time.perf_counter()
        
        output.writelines(f"Execution time: {t_fin - t_ini:0.4f} seconds")
        # output.writelines('\nCusto do melhor conjunto de parametros: ' + str(solution.fun) + '\n\n')
        output.writelines("\nSolution found: " + str(solution.x) + "\n\n")
        
        
        # utils.plot_experiment_patient(p)
        
        # # Experimental plot
        # plt.plot(t_exp[pat_cont-1], pat, 'ro')
        # plt.title(str(pat_name[pat_cont-1]))

        # plt.xlabel("dias")
        # plt.ylabel("Carga viral $log_{10}$")
        # plt.legend()
        # # # Plot da solucao com os melhores parametros
        # # cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont-1, t_exp)
        # # print('\nCusto do melhor conjunto de parametros: ' +
        # #       str(sol_pat.fun) + '\n\n')
        # # print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        # # plt.savefig("../figs/pat_"+str(pat_name[pat_cont-1])+"NGem_" +
        # #             str(maxiter)+"NPop_"+str(popsize)+".png")

        pat_cont += 1
        break
    
    output.close()