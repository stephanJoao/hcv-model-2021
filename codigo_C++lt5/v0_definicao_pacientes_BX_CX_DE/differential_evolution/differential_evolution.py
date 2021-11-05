# --- IMPORT DEPENDENCIES ------------------------------------------------------+

import random
import numpy as np
from depend.cost import viralmodelfit
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, NonlinearConstraint
import sys
import os
import time
# --- MAIN ---------------------------------------------------------------------+


# --- CONSTANTS ----------------------------------------------------------------+

cost_func = viralmodelfit                                  # Cost function
# Bounds   alpha*    r*     delta*    mu_c*  rho* , theta, sigma    c
bounds = [(20,60),(0.1,10),(0.01,2),(0.1,2),(1,15), (1,2), (1,2), (10,25)]
array_param = "alpha, r, delta, mu_c, rho, theta, sigma, c"
popsize = 10                                               # Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 10

# Vetor com todos os pacientes
patients = []
PATB07 = [4.1510]
PATB09 = [6.3541]
PATB08 = [5.7831]
PATB16 = [6.2943]
PATB17 = [6.1927]
PATB06 = [6.2584]
PATC05 = [6.208568]
PATC06 = [6.808956]
PATC09 = [6.296968]
PATC10 = [5.547272]

patients.append(PATB07)
patients.append(PATB09)
patients.append(PATB08)
patients.append(PATB16)
patients.append(PATB17)
patients.append(PATB06)
patients.append(PATC05)
patients.append(PATC06)
patients.append(PATC09)
patients.append(PATC10)

t_exp = np.zeros(10)

# --- RUN ----------------------------------------------------------------------+

if __name__ == "__main__":
    
    cwd = os.getcwd()
    saida = open(cwd+"/relatorio_DE.txt", "a")
    saida.writelines("\n\n-------NOVA TENTATIVA-------\n\n")
    saida.writelines('Population size: ' + str(popsize) +
                     '\nNumber of generations: ' + str(maxiter) + '\n')
    saida.writelines(array_param + '\n')
    saida.close()
    pat_cont = 1
    for pat in patients:
        saida = open(cwd+"/relatorio_DE.txt", "a")
        os.system("make clean")
        os.system("make")
        print(pat_cont, " Patient")
        saida.writelines(str(pat_cont) + " Patient\n\n")
        tic = time.perf_counter()
        sol_pat = differential_evolution(cost_func, bounds, args=(
            pat, 10**pat[0], pat_cont-1, t_exp), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination, tol=0.1)
        toc = time.perf_counter()
        # sol_pat.x => parametros--- sol_pat.fun => custo/ retorno da cost.py
        saida.writelines(f"Tempo de execução: {toc - tic:0.4f} segundos")
        saida.writelines(
            '\nCusto do melhor conjunto de parametros: ' + str(sol_pat.fun) + '\n\n')
        saida.writelines(
            '\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.clf()
        # Plot experimental
        plt.plot(t_exp[pat_cont-1], pat, 'ro')
        if pat_cont == 1:
            plt.title("PATB07")
        elif pat_cont == 2:
            plt.title("PATB09")
        elif pat_cont == 3:
            plt.title("PATB08")
        elif pat_cont == 4:
            plt.title("PATB16")
        elif pat_cont == 5:
            plt.title("PATB17")
        elif pat_cont == 6:
            plt.title("PATB06")
        elif pat_cont == 7:
            plt.title("PATC05")
        elif pat_cont == 8:
            plt.title("PATC06")
        elif pat_cont == 9:
            plt.title("PATC09")
        elif pat_cont == 10:
            plt.title("PATC10")

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        # Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont-1, t_exp)
        print('\nCusto do melhor conjunto de parametros: ' +
              str(sol_pat.fun) + '\n\n')
        print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.savefig("../figs/pat_"+str(pat_cont)+"NGem_" +
                    str(maxiter)+"NPop_"+str(popsize)+".png")

        pat_cont += 1
        saida.close()

    
# --- END ----------------------------------------------------------------------+
