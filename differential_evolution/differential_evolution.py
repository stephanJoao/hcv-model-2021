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

    # Compiles model
    print("========== Compiling C++ model===========")
    directory = "hcv_model/"
    os.system("make clean -C " + directory)
    os.system("make -C " + directory)
    os.system("make run -C " + directory)
    print("=========================================")

    # Header of experiment report
    output = open(dir_output + "report_DE.txt", "a")
    output.writelines("\n\n=========== Experiment information ===========\n")
    output.writelines("\nPopulation size -------- " + str(pop_size))
    output.writelines("\nNumber of generations -- " + str(max_iter) + "\n")
    output.writelines("\nParameters -- " + array_param)
    output.writelines("\nBounds ------ " + str(bounds) + "\n\n\n")
    output.writelines("\n\n")
    output.close()

    # Reads patients names
    patients = utils.reads_patients_names()
    pat_cont = 1

    for i in range(0,len(patients)):
        
        if i > 16:
            output = open(dir_output + "report_DE.txt", "a")
            output.writelines("- " + patients[i])

            # Gets patient experimental data
            exp_time, exp_viral_load = utils.reads_experimental_data(patients[i])
            
            # Diferential evolution with time measurement
            t_ini = time.perf_counter()
            solution = differential_evolution(cost.model_cost, bounds, args=(exp_time, exp_viral_load), maxiter=max_iter, popsize=pop_size, mutation=mutate, recombination=recombination, tol=0.1, workers=1)
            t_fin = time.perf_counter()
            
            utils.writes_de_parameters(patients[i], solution.x)
            
            output.writelines(f"\n  - Execution time: {t_fin - t_ini:0.4f} seconds")
            # output.writelines('\nCusto do melhor conjunto de parametros: ' + str(solution.fun) + '\n\n')
            output.writelines("\n  - Solution found: " + str(solution.x) + "\n\n\n")
            output.close()
            
            utils.plot_experiment_patient(patients[i], solution.x)
            
            pat_cont += 1    
        

    # Footer of experiment report
    output = open(dir_output + "report_DE.txt", "a")
    output.writelines("\n\n==============================================\n")
    output.close()
