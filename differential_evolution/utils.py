import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# Directory of experimental data
dir_exp_data = "differential_evolution/data/experimentalData.csv"

# Directory of differential evolution parameters
dir_de_params = "differential_evolution/output/DEs_parameters.csv"

# Directory of model's input
dir_model_input = "hcv_model/input/DE_parameters.txt"

# Directory of model's output
dir_model_output = "hcv_model/output/solution.txt"

# Directory where images will be saved
dir_images = "differential_evolution/output/images/"

def reads_patients_names():
	"""
	Gets patients names in a dataframe.
	Returns:
		l: a list that contains patients names
	"""

	df = pd.read_csv(dir_exp_data)

	n = []
	for i in range(len(df.columns) // 2):
		n.append(df.columns[i*2 + 1])

	return n

def reads_experimental_data(patient):
	"""
	Reads experimental data of selected patient.
	Arguments:
		patient: a string that contains the name of the selected patient
	Returns:
		v: a list that contains the viral load of specified patient
		t: a list that contains the time corresponding a certain viral load
	"""

	df = pd.read_csv(dir_exp_data)

	df_aux = df[patient + "_time"]
	t = df_aux.to_numpy()
	df_aux = df[patient]
	v = df_aux.to_numpy() # already in logaritmic scale

	# # Prints experimental data of selected patient
	# print("Experimental data on " + patient + ": ")
	# print(v)
	# print("\n")

	return t, v

def reads_de_parameters(patient):
	"""
	Reads the differential evolution parameters of selected patient.
	Arguments:
		dir: a string that contains the directory and name of .csv that has the parameters
		patient: a string that contains the name of the selected patient
	Returns:
		p: a list that contains the parameters
	"""

	df = pd.read_csv(dir_de_params)

	patient_index = df.index[df["Patient"] == patient]
	p = df.to_numpy()
	p = p[patient_index].flatten()
	p = np.delete(p, 0)

	# # Prints differential evolution parameters of selected patient
	# print("Parameters of " + patient + ": ")
	# print(p)
	# print("\n")

	return p

def reads_model_data():
	"""
	Reads model data of previous selected patient.
	Returns:
		v: a list that contains the viral load of specified patient in the model
		t: a list that contains the time corresponding a certain viral load in the model
	"""

	df = pd.read_csv(dir_model_output)

	df_aux = df["time"]
	t = df_aux.to_numpy()
	df_aux = df["viral_load"]
	v = df_aux.to_numpy() # not in logaritmic scale
	v = np.log10(v)

	# # Prints experimental data of selected patient
	# print("Experimental data on " + patient + ": ")
	# print(v)
	# print("\n")

	return t, v

def plot_experiment_patient(patient, de_parameters = None):
	"""
	Plots the solved model with the experimental data.
	Arguments:
		patient: a string that contains the name of the selected patient
	"""

	# Reads all necessary data
	exp_time, exp_viral_load = reads_experimental_data(patient)
	if de_parameters == None:
		de_params = reads_de_parameters(patient)
	else:
		de_params = de_parameters

	# Writes parameters on file for model solver
	epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c = de_params

	with open(dir_model_input, 'w') as params_file:
		params_file.write(str(10**exp_viral_load[0]) + "," + str(epsilon_r) + "," + str(epsilon_alpha) + "," + str(epsilon_s) + "," + str(alpha) + "," + str(r) + "," + str(delta) + "," + str(mu_c) + "," + str(rho) + "," + str(theta) + "," + str(sigma) + "," + str(c))

	# Execute the model with the parameters
	print("=========== Running C++ model ===========")
	directory = "hcv_model/"
	os.system("make clean -C " + directory)
	os.system("make -C " + directory)
	os.system("make run -C " + directory)
	print("=========================================")

	# Gets data with the solved model
	model_time, model_viral_load = reads_model_data()

	# Plots data
	plt.plot(exp_time, exp_viral_load, 'or', label="Experimental data")
	plt.plot(model_time, model_viral_load, '-b', label="Results C++ model")
	plt.title(patient)
	plt.ylabel("Viral load $log_{10}$")
	plt.xlabel("Days")
	plt.legend()
	plt.savefig(dir_images + patient + ".png", dpi=300)
	# plt.show()
	plt.clf()

if __name__ == "__main__":      
	
	# Creation of necessary directory
	if not os.path.isdir(dir_images):
		os.system("mkdir" + dir_images)	

	# #####* TESTE
	# patients = patients_names(dir_exp_data)
	# exp_viral_load, exp_time = experimental_viral_load(dir_exp_data, patients[2])
	# print(viralmodelfit(de_parameters(dir_de_params, patients[2]), exp_viral_load, exp_time))
	
	# Program
	patients = reads_patients_names()

	for p in patients:
		plot_experiment_patient(p)