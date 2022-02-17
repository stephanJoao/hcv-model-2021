import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

selected_patient = "PATB06"

# Directory of experimental data
exp_data = "differential_evolution/data/experimentalData.csv"

# Directory of diferential evolution parameters
DE_params = "differential_evolution/output/DE_parameters.csv"

# Reads experimental data
df_exp_data = pd.read_csv(exp_data)

# Gets patient's names
patients = []
for i in range(len(df_exp_data.columns) // 2):
      patients.append(df_exp_data.columns[i*2 + 1])

print("Experimental data on " + selected_patient + ": ")
print(df_exp_data[selected_patient].tolist())
print("\n")

# Reads diferential evolution parameters
df_DE_params = pd.read_csv(DE_params)

# Gets selected patient parameters
params = df_DE_params.loc[df_DE_params["Patient"] == selected_patient]
params = params.values.tolist()
params = params[0]

print("Parameters of " + selected_patient + ": ")
print(params)
print("\n")

# Writes parameters on file for model solver
pat_aux, V0, epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c = params


with open("hcv_model/input/parametros_DE.txt", 'w') as params_file:
      params_file.write(str(V0) + "," + str(epsilon_r) + "," + str(epsilon_alpha) + "," + str(epsilon_s) + "," + str(alpha) + "," + str(r) + "," + str(delta) +  "," + str(mu_c) + "," + str(rho) + "," + str(theta) + "," + str(sigma) + "," + str(c))

# Execute the model with the parameters
cwd = os.getcwd()
directory = cwd + "/hcv_model/"
os.system("make clean -C " + directory)
os.system("make -C " + directory)
os.system("make run -C " + directory)

# Gets data with the solved model
time = np.empty(0)
viral_load = np.empty(0)
with open(cwd + "/hcv_model/output/saida.txt", 'r') as solved_file:
    data_solved = [line.split(',') for line in solved_file]
    for i in data_solved:
        time = np.append(time, float(i[0]))
        viral_load = np.append(viral_load, float(i[1]))

viral_load = np.log10(viral_load)

# Gets data of selected patient
y = df_exp_data[selected_patient].tolist()
x = df_exp_data[selected_patient + "_time"].tolist()

# Plots data
plt.plot(x, y, 'or', label="experimental data")
plt.plot(time, viral_load, '-b', label="results model C++")
plt.title(selected_patient)
plt.ylabel("Viral load $log_{10}$")
plt.xlabel("Days")
plt.legend()
plt.savefig("differential_evolution/figures/" + selected_patient + ".png", dpi=300)
plt.show()