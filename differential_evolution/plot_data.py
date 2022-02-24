import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def patients_names(dir):
	"""
	Gets patients names in a dataframe.
	Arguments:
		dir: a string that contains the directory and name of .csv that has the experimental data
	Returns:
		l: a list that contains patients names
	"""

	df = pd.read_csv(dir)

	n = []
	for i in range(len(df.columns) // 2):
		n.append(df.columns[i*2 + 1])

	return n

def experimental_data(dir, patient):
	"""
	Reads experimental data of selected patient.
	Arguments:
		dir: a string that contains the directory and name of .csv that has the experimental data
		patient: a string that contains the name of the selected patient
	Returns:
		v: a list that contains the viral load of specified patient
		t: a list that contains the time corresponding a certain viral load
	"""

	df = pd.read_csv(dir)

	df_aux = df[patient]
	v = df_aux.to_numpy() # already in logaritmic scale
	df_aux = df[patient + "_time"]
	t = df_aux.to_numpy()

	# # Prints experimental data of selected patient
	# print("Experimental data on " + patient + ": ")
	# print(v)
	# print("\n")

	return v, t

def de_parameters(dir, patient):
	"""
	Reads the differential evolution parameters of selected patient.
	Arguments:
		dir: a string that contains the directory and name of .csv that has the parameters
		patient: a string that contains the name of the selected patient
	Returns:
		p: a list that contains the parameters
	"""

	df = pd.read_csv(dir)

	patient_index = df.index[df["Patient"] == patient]
	p = df.to_numpy()
	p = p[patient_index].flatten()
	p = np.delete(p, 0)

	# # Prints differential evolution parameters of selected patient
	# print("Parameters of " + patient + ": ")
	# print(p)
	# print("\n")

	return p

def plot_data_patient(dir_exp_data, dir_de_params, dir_model_input, dir_model_output, dir_images, patient):
	"""
	Plots the solved model with the experimental data.
	Arguments:
		dir_exp_data: a string that contains the directory and name of .csv that has the experimental data
		dir_de_params: a string that contains the directory and name of .csv that has the parameters
		dir_model_input: a string that contains the directory and name of .txt that that will be written in the model's input
		dir_model_input: a string that contains the directory and name of .txt that that will be read in the model's input
		dir_images: a string that contains the directory in which the plots will be saved
		patient: a string that contains the name of the selected patient
	"""

	# Reads all necessary data
	exp_data, exp_data_time = experimental_data(dir_exp_data, patient)
	de_params = de_parameters(dir_de_params, patient)

	# Writes parameters on file for model solver
	epsilon_r, epsilon_alpha, epsilon_s, alpha, r, delta, mu_c, rho, theta, sigma, c = de_params

	with open(dir_model_input, 'w') as params_file:
		params_file.write(str(10**exp_data[0]) + "," + str(epsilon_r) + "," + str(epsilon_alpha) + "," + str(epsilon_s) + "," + str(alpha) + "," + str(r) + "," + str(delta) + "," + str(mu_c) + "," + str(rho) + "," + str(theta) + "," + str(sigma) + "," + str(c))

	# Execute the model with the parameters
	print("=========== Running C++ model ===========")
	directory = "hcv_model/"
	os.system("make clean -C " + directory)
	os.system("make -C " + directory)
	os.system("make run -C " + directory)
	print("=========================================")

	# Gets data with the solved model
	time = np.empty(0)
	viral_load = np.empty(0)
	with open(dir_model_output, 'r') as solved_file:
		data_solved = [line.split(',') for line in solved_file]
	for i in data_solved:
		time = np.append(time, float(i[0]))
		viral_load = np.append(viral_load, float(i[1]))
	viral_load = np.log10(viral_load)

	# Plots data
	plt.plot(exp_data_time, exp_data, 'or', label="Experimental data")
	plt.plot(time, viral_load, '-b', label="Results C++ model")
	plt.title(patient)
	plt.ylabel("Viral load $log_{10}$")
	plt.xlabel("Days")
	plt.legend()
	plt.savefig(dir_images + patient + ".png", dpi=300)
	# plt.show()
	plt.clf()

from cost import viralmodelfit

def main():
      
	# Creation of necessary directories
	if not os.path.isdir("differential_evolution/output"): #TODO Esse vai pra differential evolution quando eu fizer
		os.system("mkdir differential_evolution/output")
	if not os.path.isdir("differential_evolution/output/images"):
		os.system("mkdir differential_evolution/output/images")

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

	# #####* TESTE
	# patients = patients_names(dir_exp_data)
	# exp_data, exp_data_time = experimental_data(dir_exp_data, patients[2])
	# print(viralmodelfit(de_parameters(dir_de_params, patients[2]), exp_data, exp_data_time))
	
	# Program
	patients = patients_names(dir_exp_data)

	for p in patients:
		plot_data_patient(dir_exp_data, dir_de_params, dir_model_input, dir_model_output, dir_images, p)


main()




