import numpy as np 
from scipy.spatial import distance
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def viralmodelfit(poi, exp, V0, pat_cont, t_exp):
    epsilon_r = poi[0]
    epsilon_alpha = poi[1]
    epsilon_s = poi[2]
    alpha = poi[3]
    r = poi[4]
    delta = poi[5]
    mu_c = poi[6]
    rho = poi[7]
    theta = poi[8]
    sigma = poi[9]
    c = poi[10]

    with open('../parametros_DE.txt', 'w') as filep:
        filep.write(str(V0)+","+str(epsilon_r)+","+str(epsilon_alpha)+","+str(epsilon_s)+
        ","+str(alpha)+","+str(r)+","+str(delta)+
        ","+str(mu_c)+","+str(rho)+","+str(theta)+
        ","+str(sigma)+","+str(c))
    #os.system("make clean")
    #os.system("make")
    os.system("make run")#Executa o modelo C++
    tempoPt = np.empty(0)
    V = np.empty(0)
    V_log = np.empty(0)
    with open("saida.txt", 'r') as f:
        lista = [line.split(',')  for line in f]
        for linha in lista:
            tempoPt = np.append(tempoPt, float(linha[0]))
            V = np.append(V, float(linha[1]))
    try:
      # Passa para a base log o resultado
      V_log = np.log10(V)
      V_pts = []
      for t in t_exp[pat_cont]:
        V_pts.append(V_log[int(t*100+1)])
      #ius = InterpolatedUnivariateSpline(t_exp[pat_cont], exp)
      #yi = ius(tempoPt)
      #plt.plot(tempoPt, yi, '--r', label='polinomio')
      
      # dst = distance.euclidean(V_log, yi) Fica muito ruim
      plt.plot(tempoPt, V_log, '-g', label='Modelo')
      # plt.plot(t_exp[pat_cont], V_pts, '^g') Prova de que esta pegando os pontos certos
      plt.plot(t_exp[pat_cont], exp, 'or', label='dados experimentais')
      # plt.show()
      dst = distance.euclidean(V_pts, exp)/len(V_pts)
    except:
      dst = 1000
    return dst