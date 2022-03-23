#include <iostream>
#include <unistd.h>

#include "HCV_Model.hpp"

using namespace std;

int main()
{  	
	char cwd[1024];	
	
	if(getcwd(cwd, sizeof(cwd)) != NULL ) 
	{
		string dir(cwd);
		string input = dir + "/input/DE_parameters.txt";
		string output = dir + "/output/solution.txt";

		HCV_Model* model = new HCV_Model(input, output);
		model->solve();
		delete model;
	}  	
  	
	return 0;
}