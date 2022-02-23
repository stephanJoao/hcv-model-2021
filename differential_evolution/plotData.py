import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

selected_patient = "PATB06"


def patients_names(dir):
      """
      Gets patients names in a dataframe.
      Arguments:
            dir: directory and name of .csv that has the experimental data
      Returns:
            l: list of patients names
      """
      df = pd.read_csv(dir)
      
      l = []
      for i in range(len(df.columns) // 2):
            l.append(df.columns[i*2 + 1])
      
      return l


# Directory of experimental data
dir_exp_data = "differential_evolution/data/experimentalData.csv"

# Directory of diferential evolution parameters
dir_de_params = "differential_evolution/output/DE_parameters.csv"

# Reads experimental data
df_exp_data = pd.read_csv(dir_exp_data)

# Gets patients' names
patients = []
for i in range(len(df_exp_data.columns) // 2):
      patients.append(df_exp_data.columns[i*2 + 1])

# Reads experimental data of selected patient
df = df_exp_data[selected_patient]
exp_data = df.to_numpy() # already in logaritmic scale
df = df_exp_data[selected_patient + "_time"]
exp_data_time = df.to_numpy()

# Prints experimental data of selected patient
print("Experimental data on " + selected_patient + ": ")
print(exp_data)
print("\n")

# Reads diferential evolution parameters
df_de_params = pd.read_csv(dir_de_params)

# Gets selected patient parameters
selected_patient_index = df_de_params.index[df_de_params["Patient"] == selected_patient]
de_params = df_de_params.to_numpy()
de_params = de_params[selected_patient_index].flatten()
de_params = np.delete(de_params, 0)

# Prints diferential evolution parameters of selected patient
print("Parameters of " + selected_patient + ": ")
print(de_params)
print("\n")

# Writes parameters on file for model solver
epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c = de_params

with open("hcv_model/input/parametros_DE.txt", 'w') as params_file:
      params_file.write(str(10**exp_data[0]) + "," + str(epsilon_r) + "," + str(epsilon_alpha) + "," + str(epsilon_s) + "," + str(alpha) + "," + str(r) + "," + str(delta) +  "," + str(mu_c) + "," + str(rho) + "," + str(theta) + "," + str(sigma) + "," + str(c))

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

# Plots data
plt.plot(exp_data_time, exp_data, 'or', label="experimental data")
plt.plot(time, viral_load, '-b', label="results model C++")
plt.title(selected_patient)
plt.ylabel("Viral load $log_{10}$")
plt.xlabel("Days")
plt.legend()
plt.savefig("differential_evolution/figures/" + selected_patient + ".png", dpi=300)
plt.show()