#include <iostream>
#include <unistd.h>

#include "HCV_Model.hpp"

using namespace std;

int main()
{  	
	char cwd[1024];	
  	getcwd(cwd, sizeof(cwd));
  	
	string dir(cwd);
	string input = dir + "/input/parametros_DE.txt";
	string output = dir + "/output/saida.txt";

	HCV_Model* model = new HCV_Model(input, output);
  	
	model->solve();
	return 0;
}