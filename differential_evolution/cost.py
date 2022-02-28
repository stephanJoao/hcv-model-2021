import numpy as np 
import os
import matplotlib.pyplot as plt
from scipy.spatial import distance

import utils

#recebe como parametro os parametros estocasticos, individuos, e os valores experimentais
def model_cost(de_params, exp_time, exp_viral_load):
	"""
	Calculates the cost of a model according to certain parameters.
	Arguments:
		de_params: the parameters for executing the model
		exp_time: the experimental viral load related time
		exp_viral_load: the experimental viral load
	Returns:
		dist: the summed distance of experimental points to the model representing the cost of the model
	"""

	# Gets parameters of the model (not the most efficient/elegant way, but its more understandable)
	epsilon_r, epsilon_alpha, epsilon_s, alpha,	r, delta, mu_c,	rho, theta, sigma, c = de_params

	with open(utils.dir_model_input, 'w') as params_file:
		params_file.write(str(10**exp_viral_load[0]) + "," + str(epsilon_r) + "," + str(epsilon_alpha) + "," + str(epsilon_s) + "," + str(alpha) + "," + str(r) + "," + str(delta) + "," + str(mu_c) + "," + str(rho) + "," + str(theta) + "," + str(sigma) + "," + str(c))
	
	# Execute the model with the parameters
	print("=========== Running C++ model ===========")
	directory = "hcv_model/"
	os.system("make clean -C " + directory)
	os.system("make -C " + directory)
	os.system("make run -C " + directory)
	print("=========================================")
	
	# Gets data with the solved model
	model_time, model_viral_load = utils.reads_model_data()
		
	# Selects those points in the model with nearest time to those in the experimental data
	model_selected = np.empty(0)
	for i in exp_time:
		try:
			model_selected = np.append(model_selected, model_viral_load[int(i * (10**2))])
		except:
			model_selected = np.append(model_selected, model_viral_load[len(model_viral_load) - 1]) # Gets last in the model

	# Plots data and model with selected points
	plt.plot(exp_time, model_selected, 'og', label='Model points selected')
	plt.plot(exp_time, exp_viral_load, 'or', label="Experimental data")
	plt.plot(model_time, model_viral_load, '-b', label="Results C++ model")
	plt.title("Visualization of cost function")
	plt.legend()
	plt.savefig(utils.dir_images + "test.png", dpi=300)
	# plt.show()
	plt.clf()

	# Calculates the cost of the model
	try:
		dist = 0
		for i in range(0, len(exp_time)):
			dist = dist + distance.euclidean([exp_time[i], model_selected[i]], [exp_time[i], exp_viral_load[i]])
	except:
		dist = 100

	return dist

if __name__ == "__main__":
	patients = utils.reads_patients_names()
	
	for p in patients:
			exp_time, exp_viral_load = utils.reads_experimental_data(p)
			de_params = utils.reads_de_parameters(p)

			cost = model_cost(de_params, exp_time, exp_viral_load) 

			print("Cost of patient " + p + " = " + str(cost) + "\n")