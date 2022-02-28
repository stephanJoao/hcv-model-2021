import random
import numpy as np
from depend.cost import viralmodelfit
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, NonlinearConstraint
import os
import time



cost_func = viralmodelfit                                  #Cost function
bounds = [(0.1,0.99),(0.1,0.99),(0.1,0.99),(20,60),(0.1,10),(0.01,2),(0.1,2),(1,15), (1,2), (1,2), (10,25)] 
array_param = "epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c" 
popsize = 30                                               #Population size
# Mutation factor [0,2]
mutate = 0.7
# Recombination rate [0,1]
recombination = 0.5
# Max number of generations (maxiter)
maxiter = 30

t_exp = [np.array([0.0, 0.04, 0.08, 0.17, 0.34, 0.50, 1.00, 1.50, 3.01, 6.01, 7.00, 9.99, 14.00, 20.98, 30.01, 44.01, 55.01]) #PATB16
        ,np.array([0.0, 0.04, 0.08, 0.17, 0.33, 0.50, 1.00, 1.50, 3.00, 3.95, 6.98, 9.98, 15.97, 20.99, 28.02, 42.99,]) #PATB10
        ,np.array([0,0.04,0.08,0.17,0.34,0.50,1.00,1.50,2.97,5.99,7.97,10.06,14.06,21.02, 28.04, 44.07]) #PATC14
        ]
patients = []
PATB16 = [6.2943, 6.3541, 6.3212, 5.5915, 4.1949, 3.8517, 3.6651, 3.4814, 3.2529, 3.0120, 3.0302, 2.7528, 2.3838, 2.1818, 1.9243, 1.7160, 1]
PATB10 = [6.3840, 6.4587, 6.4587, 6.1718, 5.4932, 4.9461, 3.5412, 3.1535, 2.4669, 2.3945, 1.9294, 1.6990, 1, 1, 1, 1]
PATC14 = [5.9922,5.8673,5.8764,5.6997,5.0688,4.5903,4.0052,3.3286,2.6702,1.8865,1.5185, 1, 1, 1, 1, 1]

patients.append(PATB16)
patients.append(PATB10)
patients.append(PATC14)

pat_name = []
pat_name.append("PATB16")
pat_name.append("PATB10")
pat_name.append("PATC14")

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
        saida.writelines(str(pat_name[pat_cont-1]) + " Patient\n\n")
        tic = time.perf_counter()
        sol_pat = differential_evolution(cost_func, bounds, args=(pat, 10**pat[0], pat_cont-1, t_exp), maxiter=maxiter, popsize=popsize, mutation=mutate, recombination=recombination)#, tol=0.1)
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
        plt.title(str(pat_name[pat_cont-1]))

        plt.xlabel("dias")
        plt.ylabel("Carga viral $log_{10}$")
        plt.legend()
        # Plot da solucao com os melhores parametros
        cost_func(sol_pat.x, pat, (10**pat[0]), pat_cont-1, t_exp)
        print('\nCusto do melhor conjunto de parametros: ' +
              str(sol_pat.fun) + '\n\n')
        print('\nO melhor conjunto de parametros: '+str(sol_pat.x) + '\n\n')
        plt.savefig("../figs/pat_"+str(pat_name[pat_cont-1])+"NGem_" +
                    str(maxiter)+"NPop_"+str(popsize)+".png")

        pat_cont += 1
        saida.close()

    
