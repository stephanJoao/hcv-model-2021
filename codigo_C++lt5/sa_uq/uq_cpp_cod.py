import uncertainpy as un
import chaospy as cp
import numpy as np
import time
import os


def viralmodelfit(alpha, r, rho, epsilon_s, epsilon_alpha, epsilon_r, delta):
    
    #Escreve os parametros gerados pela un no arquivo a ser lido pelo modelo em C++
    with open('parametros.txt', 'w') as filep:
        filep.write(str(alpha)+","+str(r)+","+str(rho)+","+str(epsilon_s)+","+str(epsilon_alpha)+","+str(epsilon_r)+","+str(delta))
    os.system("make run")#Executa o modelo C++
    
    #Le a saida do c++ e retorna para a un
    tempo = np.empty(0)
    viral_load = np.empty(0)
    with open("saida.txt", 'r') as f:
        lista = [line.split(' ')  for line in f]
        for linha in lista:
            tempo = np.append(tempo, float(linha[0]))
            viral_load = np.append(viral_load, float(linha[1]))
            
    return tempo, viral_load

##############FIM MODELO################


model = un.Model(
    run = viralmodelfit,
    labels=["Tempo (dias)",
        "Carga viral (log10)"]
)

# create distributions

print('Criando distribuicoes normais')
#epsilon_r = cp.Normal(0.3,0.1)
epsilon_r = cp.Normal(0.3,0.001)
epsilon_s = cp.Normal(0.998,0.001)
#epsilon_alpha = cp.Normal(0.92,0.002)
epsilon_alpha = cp.Normal(0.92,0.001)
delta = cp.Normal(0.07,0.01)
alpha = cp.Normal(22,0.001)
r = cp.Normal(11.1,0.1)
rho = cp.Normal(12.0,1.0)

print('Distribuicoes criadas')



# define parameter dictionary
parameters = {"alpha": alpha,
        "r": r,
        "rho": rho,
        "epsilon_s": epsilon_s,
        "epsilon_alpha": epsilon_alpha,
        "epsilon_r": epsilon_r,
        "delta": delta
            }

# set up UQ
UQ = un.UncertaintyQuantification(
    model=model,
    parameters=parameters
)
#Garante que o modeloC++ está na ultima versao
os.system("make")

print('Inicio UQ --- ')
begin = time.perf_counter()
data = UQ.monte_carlo(nr_samples=10)
end = time.perf_counter()
print(' --- Fim UQ')
tempo = end-begin
with open('tempo.txt', 'w') as file:
    file.write("Tempo de execução: "+str(tempo)+" segundos")